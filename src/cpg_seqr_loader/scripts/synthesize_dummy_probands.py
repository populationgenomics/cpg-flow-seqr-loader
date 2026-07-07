"""
Synthesize a 'dummy' male proband gVCF from a pair of parental gVCFs (RD-1168).

For cohorts containing only unaffected parental duos (no probands), we can enable
trio analysis in seqr we create a synthetic proband that inherits, at worst case, every variant
its parents carry:

    autosomes / chrX PAR : maternal allele = ALT iff mother carries any ALT
                           paternal allele = ALT iff father carries any ALT
                           e.g. 0/1 + 0/1 -> 1/1 ; 0/0 + 0/1 -> 0/1 ; 0/0 + 0/0 -> 0/0
    non-PAR chrX (male)  : hemizygous, maternal-only -> haploid 1 iff mother carries ALT
    chrY                 : copy the father's records verbatim (already hemizygous)
    chrM                 : copy the mother's records verbatim (maternal inheritance)

The output is a spec-valid GATK-style gVCF *with <NON_REF> reference blocks*, so that once it is
combined with the real samples it reads as confident hom-ref (0/0) - not missing - at sites that
are variant only in other families. The dummy's reference blocks are the intersection of the two
parents' hom-ref bands, split at the union of their variant sites.

This phase is gVCF synthesis only. Metamist registration (participant/sample/SG creation,
pedigree import, gvcf analysis linking) and the VDS/seqr workflow wiring are handled separately.
"""

import gzip
import shutil
import tempfile
from argparse import ArgumentParser
from pathlib import Path

import pysam
from cpg_utils import to_path
from loguru import logger

# symbolic non-reference alleles used by GATK gVCFs
NON_REF_ALLELES = {'<NON_REF>', '<*>'}

# GRCh38 pseudoautosomal regions on chrX (1-based, inclusive). Inside these the dummy is diploid;
# outside them (non-PAR chrX) the male dummy is hemizygous and inherits only from the mother.
PAR_X_REGIONS = [
    (10_001, 2_781_479),  # PAR1
    (155_701_383, 156_030_895),  # PAR2
]

# contig name groups (handle both 'chr'-prefixed and bare names)
CHRX_NAMES = {'chrX', 'X'}
CHRY_NAMES = {'chrY', 'Y'}
CHRM_NAMES = {'chrM', 'chrMT', 'M', 'MT'}

# FORMAT/INFO/ALT meta we rely on - added defensively if a parent header omits them
REQUIRED_FORMAT = {
    'GT': ('1', 'String', 'Genotype'),
    'GQ': ('1', 'Integer', 'Genotype Quality'),
    'DP': ('1', 'Integer', 'Approximate read depth'),
    'MIN_DP': ('1', 'Integer', 'Minimum DP observed within the GVCF block'),
    'AD': ('R', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed'),
    'PL': ('G', 'Integer', 'Normalized, Phred-scaled likelihoods for genotypes'),
}


def is_ref_block(rec: pysam.VariantRecord) -> bool:
    """A gVCF reference block has a single <NON_REF> ALT and an END-defined span."""
    return rec.alts is None or (len(rec.alts) == 1 and rec.alts[0] in NON_REF_ALLELES)


def is_par(pos: int) -> bool:
    """Is this 1-based chrX position within a pseudoautosomal region?"""
    return any(start <= pos <= end for start, end in PAR_X_REGIONS)


def carried_alt(rec: pysam.VariantRecord | None) -> str | None:
    """Return the (first) real ALT allele the single sample carries, or None if hom-ref/missing.

    A reference block (GT 0/0) returns None, which is exactly the 'contributes a ref allele'
    behaviour the worst-case rule needs.
    """
    if rec is None:
        return None
    gt = rec.samples[0].get('GT')
    if gt is None:
        return None
    for allele_index in gt:
        if allele_index in (0, None):
            continue
        allele = rec.alleles[allele_index]
        if allele in NON_REF_ALLELES:
            continue
        return allele
    return None


def _reanchor_alt(alt: str | None, parent_ref: str, common_ref: str) -> str | None:
    """Re-anchor an ALT allele defined against parent_ref onto a (longer) common_ref.

    When two parents have variants at the same locus with different REF lengths, all alleles must
    be expressed against one common REF (the longest). An ALT defined against a shorter REF needs
    the extra reference bases appended to it. parent_ref is always a prefix of common_ref (both are
    the reference sequence starting at the same position), so appending common_ref's suffix is the
    correct VCF allele-merge operation (the same thing bcftools norm / GATK do).

    Example: father C->CAAAAT against common REF CT becomes CAAAAT + 'T' = CAAAATT (i.e. CT->CAAAATT),
    which is the father's insertion correctly expressed without consuming the reference T.
    """
    if alt is None:
        return None
    if len(parent_ref) < len(common_ref):
        return alt + common_ref[len(parent_ref) :]
    return alt


def _int(value) -> int:
    """Coerce an optional FORMAT scalar (which may be None or a 1-tuple) to a non-negative int."""
    if value is None:
        return 0
    if isinstance(value, tuple):
        value = value[0] if value else 0
    return int(value or 0)


def gq_of(rec: pysam.VariantRecord) -> int:
    return _int(rec.samples[0].get('GQ'))


def dp_of(rec: pysam.VariantRecord) -> int:
    sample = rec.samples[0]
    return _int(sample.get('DP') or sample.get('MIN_DP'))


def min_dp_of(rec: pysam.VariantRecord) -> int:
    sample = rec.samples[0]
    return _int(sample.get('MIN_DP') or sample.get('DP'))


class PeekIter:
    """A single-item lookahead wrapper over a record iterator."""

    def __init__(self, iterator) -> None:
        self._iter = iter(iterator)
        self._head = next(self._iter, None)

    def peek(self) -> pysam.VariantRecord | None:
        return self._head

    def pop(self) -> pysam.VariantRecord | None:
        head = self._head
        self._head = next(self._iter, None)
        return head


def build_output_header(template: pysam.VariantHeader, sample_name: str) -> pysam.VariantHeader:
    """Clone a parent header (contigs + meta lines) into a single-sample output header."""
    header = pysam.VariantHeader()
    for record in template.records:
        header.add_record(record)

    # GATK gVCF headers already declare these; add them defensively in case a template omits one.
    # Adding a meta line that already exists is a no-op, so this is safe.
    if 'NON_REF' not in header.alts:
        header.add_meta('ALT', items=[('ID', 'NON_REF'), ('Description', 'Represents any possible alternative allele')])
    if 'END' not in header.info:
        header.add_meta(
            'INFO',
            items=[
                ('ID', 'END'),
                ('Number', '1'),
                ('Type', 'Integer'),
                ('Description', 'Stop position of the interval'),
            ],
        )
    for key, (number, vtype, desc) in REQUIRED_FORMAT.items():
        if key not in header.formats:
            header.add_meta('FORMAT', items=[('ID', key), ('Number', number), ('Type', vtype), ('Description', desc)])

    header.add_sample(sample_name)
    return header


def clone_record(
    header: pysam.VariantHeader,
    sample_name: str,
    src: pysam.VariantRecord,
) -> pysam.VariantRecord:
    """Copy a source record into the single-sample output header (used for chrY / chrM).

    Site-level INFO annotations are deliberately dropped: the Hail gVCF combiner does not need
    them (site annotations are recomputed downstream via VQSR/Hail), and copying them verbatim is
    fragile across headers (e.g. a multi-valued MLEAC/MLEAF from a multiallelic site fails against
    the mother-derived output header). END is preserved via stop=, matching the rest of the dummy
    gVCF, which also carries no INFO. Only FORMAT fields the output header declares are copied.
    """
    out = header.new_record(
        contig=src.contig,
        start=src.start,
        stop=src.stop,
        alleles=src.alleles,
        id=src.id,
        qual=src.qual,
    )
    out_sample = out.samples[sample_name]
    for key, value in src.samples[0].items():
        if key in header.formats:
            out_sample[key] = value
    return out


def emit_variant(
    out_vcf: pysam.VariantFile,
    header: pysam.VariantHeader,
    sample_name: str,
    contig: str,
    pos: int,
    ref: str,
    maternal_alt: str | None,
    paternal_alt: str | None,
    gq: int,
    dp: int,
    *,
    haploid: bool,
) -> None:
    """Write one dummy variant record applying the worst-case inheritance rule at this site."""
    alts: list[str] = []
    for alt in (maternal_alt, paternal_alt):
        if alt and alt not in alts:
            alts.append(alt)
    alleles = [ref, *alts, '<NON_REF>']

    def index_of(allele: str | None) -> int:
        return alleles.index(allele) if allele else 0

    maternal_index = index_of(maternal_alt)
    paternal_index = index_of(paternal_alt) if not haploid else 0
    genotype = (maternal_index,) if haploid else tuple(sorted((maternal_index, paternal_index)))

    out = header.new_record(contig=contig, start=pos - 1, stop=pos - 1 + len(ref), alleles=alleles)
    out.qual = None
    sample = out.samples[sample_name]
    sample['GT'] = genotype
    sample.phased = False
    sample['GQ'] = gq
    if dp:
        sample['DP'] = dp

    n_alleles = len(alleles)
    allelic_depths = [0] * n_alleles
    present = list(genotype)
    if any(index != 0 for index in present):
        per = dp // len(present) if dp else 0
        for index in present:
            allelic_depths[index] += per
        if dp:
            allelic_depths[present[0]] += dp - per * len(present)
    else:
        allelic_depths[0] = dp
    sample['AD'] = allelic_depths

    # PL consistent with the called genotype (0 at the call, GQ everywhere else)
    if haploid:
        likelihoods = [gq] * n_alleles
        likelihoods[genotype[0]] = 0
    else:
        n_genotypes = n_alleles * (n_alleles + 1) // 2
        low, high = sorted(genotype)
        called = high * (high + 1) // 2 + low
        likelihoods = [gq] * n_genotypes
        likelihoods[called] = 0
    sample['PL'] = likelihoods

    out_vcf.write(out)


def emit_ref_block(
    out_vcf: pysam.VariantFile,
    header: pysam.VariantHeader,
    sample_name: str,
    contig: str,
    start: int,
    end: int,
    ref_base: str,
    gq: int,
    min_dp: int,
    *,
    haploid: bool,
) -> None:
    """Write one dummy <NON_REF> reference block spanning [start, end] (1-based inclusive)."""
    # passing stop=end makes pysam emit INFO/END for the block (END is a reserved field set via stop)
    out = header.new_record(contig=contig, start=start - 1, stop=end, alleles=(ref_base, '<NON_REF>'))
    out.qual = None
    sample = out.samples[sample_name]
    sample['GT'] = (0,) if haploid else (0, 0)
    sample.phased = False
    sample['GQ'] = gq
    if min_dp:
        sample['MIN_DP'] = min_dp
    out_vcf.write(out)


def merge_contig(
    out_vcf: pysam.VariantFile,
    header: pysam.VariantHeader,
    sample_name: str,
    contig: str,
    mother_iter: 'PeekIter',
    father_iter: 'PeekIter',
    *,
    is_chrx: bool,
) -> None:
    """Two-pointer interval sweep over both parents' records for one autosome (or chrX).

    Streams directly from the sorted iterators (only the current record per parent is held in
    memory at any time) so whole-genome contigs with millions of records do not exhaust RAM.
    """
    while _on_contig(mother_iter, contig) and _on_contig(father_iter, contig):
        mother, father = mother_iter.peek(), father_iter.peek()
        m_start, m_end = mother.pos, mother.stop
        f_start, f_end = father.pos, father.stop

        seg_start = max(m_start, f_start)
        seg_end = min(m_end, f_end)

        # no overlap (a single-parent coverage gap) - drop the earlier interval, dummy is missing there
        if seg_start > seg_end:
            if m_end < f_start:
                mother_iter.pop()
            else:
                father_iter.pop()
            continue

        mother_var = not is_ref_block(mother)
        father_var = not is_ref_block(father)
        haploid = is_chrx and not is_par(seg_start)

        # On non-PAR chrX the son inherits only the maternal X - the father is ignored entirely.
        if haploid:
            if mother_var and mother.pos == seg_start:
                emit_variant(
                    out_vcf,
                    header,
                    sample_name,
                    contig,
                    seg_start,
                    mother.ref,
                    carried_alt(mother),
                    None,
                    gq_of(mother),
                    dp_of(mother),
                    haploid=True,
                )
            elif not mother_var:
                ref_base = mother.ref[0] if m_start == seg_start else father.ref[0]
                emit_ref_block(
                    out_vcf,
                    header,
                    sample_name,
                    contig,
                    seg_start,
                    seg_end,
                    ref_base,
                    gq_of(mother),
                    min_dp_of(mother),
                    haploid=True,
                )
            # else: indel-continuation sub-segment of a maternal variant - already emitted at its pos
        elif mother_var or father_var:
            # variant site - emit once, at the variant's own POS (dedupes indel-spanning segments)
            mother_starts_here = mother_var and m_start == seg_start
            father_starts_here = father_var and f_start == seg_start
            if mother_starts_here or father_starts_here:
                # A parent contributes its ALT only if its variant record STARTS here - a parent
                # whose indel merely spans this position was already emitted at its own start, so
                # including it here would create a spurious (often REF==ALT) allele.
                # The common REF is the longest contributing REF; each ALT is re-anchored to it
                # (see _reanchor_alt) so a multiallelic merge of differently-anchored indels stays
                # normalised and left-aligned for sparse_split_multi.
                contributing_refs = [
                    rec.ref for rec, starts in ((mother, mother_starts_here), (father, father_starts_here)) if starts
                ]
                common_ref = max(contributing_refs, key=len)
                maternal_alt = (
                    _reanchor_alt(carried_alt(mother), mother.ref, common_ref) if mother_starts_here else None
                )
                paternal_alt = (
                    _reanchor_alt(carried_alt(father), father.ref, common_ref) if father_starts_here else None
                )
                emit_variant(
                    out_vcf,
                    header,
                    sample_name,
                    contig,
                    seg_start,
                    common_ref,
                    maternal_alt,
                    paternal_alt,
                    min(gq_of(mother), gq_of(father)),
                    min(dp_of(mother), dp_of(father)),
                    haploid=False,
                )
        else:
            ref_base = mother.ref[0] if m_start == seg_start else father.ref[0]
            emit_ref_block(
                out_vcf,
                header,
                sample_name,
                contig,
                seg_start,
                seg_end,
                ref_base,
                min(gq_of(mother), gq_of(father)),
                min(min_dp_of(mother), min_dp_of(father)),
                haploid=False,
            )

        if m_end == seg_end:
            mother_iter.pop()
        if f_end == seg_end:
            father_iter.pop()

    # discard any single-parent tail left on this contig (dummy is missing there)
    _skip_contig(mother_iter, contig)
    _skip_contig(father_iter, contig)


def synthesize_dummy_gvcf(mother_gvcf: str, father_gvcf: str, sample_name: str, out_path: str) -> None:
    """Generate the dummy proband gVCF at out_path from the two parental gVCFs."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        local_mother = _localise(mother_gvcf, tmp)
        local_father = _localise(father_gvcf, tmp)
        # pysam forces INFO/END onto every record it writes (END is declared in the header), so we
        # write to a raw file first then strip END off variant lines in a text pass (see below).
        local_raw = str(tmp / 'raw.g.vcf.gz')
        local_out = str(tmp / 'dummy.g.vcf.gz')

        with pysam.VariantFile(local_mother) as mother_vcf, pysam.VariantFile(local_father) as father_vcf:
            header = build_output_header(mother_vcf.header, sample_name)
            contig_rank = {contig: rank for rank, contig in enumerate(header.contigs)}

            mother_iter = PeekIter(mother_vcf)
            father_iter = PeekIter(father_vcf)

            with pysam.VariantFile(local_raw, 'wz', header=header) as out_vcf:
                while mother_iter.peek() is not None or father_iter.peek() is not None:
                    contig = _next_contig(mother_iter, father_iter, contig_rank)
                    logger.info(f'Processing {contig}')

                    if contig in CHRY_NAMES:
                        # son inherits the entire paternal Y - copy father's records verbatim
                        _copy_contig(father_iter, contig, out_vcf, header, sample_name)
                        _skip_contig(mother_iter, contig)  # discard any anomalous maternal chrY
                    elif contig in CHRM_NAMES:
                        # mitochondria are maternally inherited - copy mother's records verbatim
                        _copy_contig(mother_iter, contig, out_vcf, header, sample_name)
                        _skip_contig(father_iter, contig)
                    else:
                        merge_contig(
                            out_vcf,
                            header,
                            sample_name,
                            contig,
                            mother_iter,
                            father_iter,
                            is_chrx=contig in CHRX_NAMES,
                        )

        # strip the spurious INFO/END off variant records (ref blocks keep it), then bgzip + index
        _strip_end_from_variants(local_raw, local_out, tmp)
        pysam.tabix_index(local_out, preset='vcf', force=True)
        _delocalise(local_out, out_path)
        _delocalise(local_out + '.tbi', out_path + '.tbi')
        logger.info(f'Wrote dummy proband gVCF for {sample_name} to {out_path}')


def _strip_end_from_variants(raw_path: str, out_path: str, tmp: Path) -> None:
    """Remove INFO/END from variant records, keeping it only on <NON_REF>-only reference blocks.

    pysam writes INFO/END on every record (it is declared in the header), but the Hail combiner
    treats any record with a defined END as a reference block and errors if such a record has a
    non-hom-ref genotype. A record is a reference block iff its ALT is exactly <NON_REF>; every
    other record is a variant and must have END absent so the combiner reads it correctly.
    """
    plain = str(tmp / 'stripped.vcf')
    with gzip.open(raw_path, 'rt') as fin, Path(plain).open('w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            columns = line.rstrip('\n').split('\t')
            alt, info = columns[4], columns[7]
            if alt == '<NON_REF>' or 'END=' not in info:
                fout.write(line)
                continue
            kept = [field for field in info.split(';') if not field.startswith('END=')]
            columns[7] = ';'.join(kept) if kept else '.'
            fout.write('\t'.join(columns) + '\n')
    pysam.tabix_compress(plain, out_path, force=True)


def _next_contig(mother_iter: PeekIter, father_iter: PeekIter, contig_rank: dict[str, int]) -> str:
    """Pick the next contig to process - the lower-ranked of the two streams' current heads."""
    heads = [it.peek() for it in (mother_iter, father_iter) if it.peek() is not None]
    return min((h.contig for h in heads), key=lambda c: contig_rank.get(c, len(contig_rank)))


def _on_contig(iterator: PeekIter, contig: str) -> bool:
    """Is the iterator's current record on the given contig?"""
    head = iterator.peek()
    return head is not None and head.contig == contig


def _skip_contig(iterator: PeekIter, contig: str) -> None:
    """Discard all leading records on the given contig from the iterator."""
    while _on_contig(iterator, contig):
        iterator.pop()


def _copy_contig(
    iterator: PeekIter,
    contig: str,
    out_vcf: pysam.VariantFile,
    header: pysam.VariantHeader,
    sample_name: str,
) -> None:
    """Clone every leading record on the given contig into the output, renamed to the dummy sample."""
    while _on_contig(iterator, contig):
        out_vcf.write(clone_record(header, sample_name, iterator.pop()))


def _localise(path: str, tmp: Path) -> str:
    """Download a remote (gs://) gVCF to a local temp path; pass local paths through unchanged."""
    if '://' not in path:
        return path
    local = tmp / Path(path).name
    logger.info(f'Localising {path} -> {local}')
    to_path(path).copy(local)
    return str(local)


def _delocalise(local_path: str, out_path: str) -> None:
    """Upload a local file to a remote (gs://) destination, or copy locally."""
    if '://' in out_path:
        to_path(out_path).upload_from(local_path)
    elif local_path != out_path:
        shutil.copy(local_path, out_path)


def cli_main():
    """CLI entrypoint."""
    parser = ArgumentParser(description='Synthesize a dummy male proband gVCF from two parental gVCFs')
    parser.add_argument('--mother_gvcf', required=True, help='Path to the maternal gVCF (.g.vcf.gz)')
    parser.add_argument('--father_gvcf', required=True, help='Path to the paternal gVCF (.g.vcf.gz)')
    parser.add_argument('--output', required=True, help='Output dummy gVCF path (.g.vcf.gz)')
    parser.add_argument('--sample_name', required=True, help='Sample (SM) name for the dummy proband')
    args = parser.parse_args()
    synthesize_dummy_gvcf(
        mother_gvcf=args.mother_gvcf,
        father_gvcf=args.father_gvcf,
        sample_name=args.sample_name,
        out_path=args.output,
    )


if __name__ == '__main__':
    cli_main()
