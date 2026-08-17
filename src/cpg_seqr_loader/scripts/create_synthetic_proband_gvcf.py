"""
Synthesize a synthetic male proband gVCF from a pair of parental gVCFs (RD-1168).

For cohorts containing only unaffected parental duos (no probands), we can enable
trio analysis in seqr by creating a synthetic proband that inherits, at worst case, every variant
its parents carry:

    autosomes / chrX PAR : maternal allele = ALT iff mother carries any ALT
                           paternal allele = ALT iff father carries any ALT
                           e.g. 0/1 + 0/1 -> 1/1 ; 0/0 + 0/1 -> 0/1 ; 0/0 + 0/0 -> 0/0
    non-PAR chrX (male)  : hemizygous, maternal-only -> haploid 1 iff mother carries ALT
    chrY                 : copy the father's records verbatim (already hemizygous)
    chrM                 : copy the mother's records verbatim (maternal inheritance)

The output is a spec-valid GATK-style gVCF *with <NON_REF> reference blocks*, so that once it is
combined with the real samples it reads as confident hom-ref (0/0) - not missing - at sites that
are variant only in other families. The synthetic proband's reference blocks are the intersection
of the two parents' hom-ref bands, split at the union of their variant sites.

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

# GRCh38 pseudoautosomal regions on chrX (1-based, inclusive). Inside these the synthetic proband
# is diploid; outside them (non-PAR chrX) the male synthetic proband is hemizygous and inherits
# only from the mother.
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


def carried_alts(rec: pysam.VariantRecord | None) -> list[str]:
    """Return every real ALT allele the single sample carries, in GT order.

    A reference block (GT 0/0) returns [], which is exactly the 'contributes a ref allele'
    behaviour the worst-case rule needs. A 1/2 heterozygous call returns both non-ref alleles;
    a 1/1 returns the ALT twice (kept as-is; callers can dedupe if they need to).
    """
    if rec is None:
        return []
    gt = rec.samples[0].get('GT')
    if gt is None:
        return []
    alts: list[str] = []
    for allele_index in gt:
        if allele_index in (0, None):
            continue
        allele = rec.alleles[allele_index]
        if allele in NON_REF_ALLELES:
            continue
        alts.append(allele)
    return alts


def carried_alt(rec: pysam.VariantRecord | None) -> str | None:
    """Convenience wrapper: the first carried real ALT allele, or None. Used on non-PAR chrX where
    inheritance is single-parent (maternal), so picking one alt is unambiguous.
    """
    alts = carried_alts(rec)
    return alts[0] if alts else None


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


def gq_from_pl(rec: pysam.VariantRecord) -> int:
    """Return the sample's PL-derived GQ (second-smallest PL minus smallest PL).

    The Hail gVCF combiner (and therefore seqr's loader) recomputes GQ from PL during
    sparse_split_multi, discarding whatever value is in the stored GQ field. GATK's
    CalculateGenotypePosteriors can inflate the stored GQ above the raw-PL value by folding in
    posterior probabilities (GP), so using stored GQ here would make the synthetic proband appear
    more confident than the contributor parent will be shown as in seqr. Deriving from PL directly
    keeps the synthetic proband's displayed GQ aligned with the parent's.

    Falls back to the stored GQ field if PL is missing or has fewer than 2 entries.
    """
    pl = rec.samples[0].get('PL')
    if pl is None:
        return gq_of(rec)
    values = sorted(int(x) for x in pl if x is not None)
    if len(values) < 2:
        return gq_of(rec)
    return max(0, values[1] - values[0])


def dp_of(rec: pysam.VariantRecord) -> int:
    sample = rec.samples[0]
    return _int(sample.get('DP') or sample.get('MIN_DP'))


def min_dp_of(rec: pysam.VariantRecord) -> int:
    sample = rec.samples[0]
    return _int(sample.get('MIN_DP') or sample.get('DP'))


class PeekIter:
    """Buffered iterator with a small look-ahead window.

    GATK HaplotypeCaller can emit overlapping variant records within a single sample's gVCF at
    complex tandem-repeat loci - typically an outer spanning deletion plus phased interior
    insertions/SNVs. A single-record peek would only see the outer record while it's active, so
    the interior contribution to the synthetic proband is silently dropped (the outer pops after the sweep
    has already advanced past the interior's position). The buffered window lets `peek_starter`
    look ahead a few records to find an interior variant that starts at the current segment.

    Window of 16 comfortably covers typical HaplotypeCaller TR clusters; scan cost is bounded by
    the early-exit on sorted POS, so a larger window has no runtime penalty in practice.
    """

    def __init__(self, iterator, window: int = 16) -> None:
        self._iter = iter(iterator)
        self._buffer: list[pysam.VariantRecord] = []
        self._window = window
        self._fill()

    def _fill(self) -> None:
        while len(self._buffer) < self._window:
            nxt = next(self._iter, None)
            if nxt is None:
                break
            self._buffer.append(nxt)

    def peek(self) -> pysam.VariantRecord | None:
        return self._buffer[0] if self._buffer else None

    def peek_starter(self, contig: str, pos: int) -> tuple[int, pysam.VariantRecord] | None:
        """Find a buffered variant record whose POS matches `pos` on `contig`.

        Returns `(buffer_index, record)` for the first hit, else `None`. Ref blocks are skipped
        (they can't contribute an ALT). Uses the sort-by-POS invariant to early-exit once we've
        passed `pos`.
        """
        for idx, rec in enumerate(self._buffer):
            if rec.contig != contig:
                break
            if rec.pos > pos:
                break
            if rec.pos == pos and not is_ref_block(rec):
                return idx, rec
        return None

    def pop(self) -> pysam.VariantRecord | None:
        if not self._buffer:
            return None
        head = self._buffer.pop(0)
        self._fill()
        return head

    def pop_at(self, index: int) -> None:
        """Remove a specific buffered record (used to consume out-of-order starters)."""
        if 0 <= index < len(self._buffer):
            self._buffer.pop(index)
            self._fill()


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
    the mother-derived output header). END is preserved via stop=, matching the rest of the
    synthetic gVCF, which also carries no INFO. Only FORMAT fields the output header declares are copied.
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
    """Write one synthetic variant record applying the worst-case inheritance rule at this site."""
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
    """Write one synthetic <NON_REF> reference block spanning [start, end] (1-based inclusive)."""
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


def merge_contig(  # noqa: PLR0915 - complex two-pointer sweep with several interior branches; splitting hurts readability more than it helps
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
        # _on_contig guarantees peek() is not None; raise explicitly so mypy can narrow the
        # types (and so this survives `python -O` unlike an assert).
        if mother is None or father is None:
            raise RuntimeError('_on_contig invariant violated: peek returned None')
        m_start, m_end = mother.pos, mother.stop
        f_start, f_end = father.pos, father.stop

        seg_start = max(m_start, f_start)
        seg_end = min(m_end, f_end)

        # no overlap (a single-parent coverage gap) - drop the earlier interval, synthetic proband is missing there
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
                    gq_from_pl(mother),
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

            # Look-ahead: if a parent's peek is a spanning outer variant that doesn't start here,
            # check their buffer for an interior variant that does start here. GATK HaplotypeCaller
            # emits overlapping records in complex TR loci (outer deletion + phased interior
            # insertion/SNV); without look-ahead the interior would be silently dropped because
            # the outer pops after the sweep has already advanced past the interior's position.
            mother_starter_result = None if mother_starts_here else mother_iter.peek_starter(contig, seg_start)
            father_starter_result = None if father_starts_here else father_iter.peek_starter(contig, seg_start)
            mother_contributor = (
                mother if mother_starts_here else (mother_starter_result[1] if mother_starter_result else None)
            )
            father_contributor = (
                father if father_starts_here else (father_starter_result[1] if father_starter_result else None)
            )

            if mother_contributor is not None or father_contributor is not None:
                # A parent contributes its ALT only if a variant record of theirs starts here (via
                # the current peek or a buffered interior). A parent whose indel merely spans this
                # position was already emitted at its own start; including it here would create a
                # spurious (often REF==ALT) allele.
                # The common REF is the longest contributing REF; each ALT is re-anchored to it
                # (see _reanchor_alt) so a multiallelic merge of differently-anchored indels stays
                # normalised and left-aligned for sparse_split_multi.
                contributing_refs = [rec.ref for rec in (mother_contributor, father_contributor) if rec is not None]
                common_ref = max(contributing_refs, key=len)

                # Collect each contributing parent's real ALTs (re-anchored to common_ref), then
                # prefer a shared ALT when both parents contribute one. This handles the case where
                # a parent is 1/2 (carries two different ALTs) and the "correct" worst-case pairing
                # is with the parent's ALT that matches the other parent - not just their first ALT.
                # Falls back to the first-listed ALT when there is no overlap.
                m_alts = (
                    [_reanchor_alt(a, mother_contributor.ref, common_ref) for a in carried_alts(mother_contributor)]
                    if mother_contributor is not None
                    else []
                )
                f_alts = (
                    [_reanchor_alt(a, father_contributor.ref, common_ref) for a in carried_alts(father_contributor)]
                    if father_contributor is not None
                    else []
                )
                shared = next((a for a in m_alts if a in f_alts), None)
                maternal_alt: str | None
                paternal_alt: str | None
                if shared is not None:
                    maternal_alt = paternal_alt = shared
                else:
                    maternal_alt = m_alts[0] if m_alts else None
                    paternal_alt = f_alts[0] if f_alts else None

                # GQ/DP: only the contributing parents' variant-record values. min()-ing across
                # a non-contributing parent's ref block drags the synthetic proband's quality down to
                # the block's band-wide minimum (systematically lower than the contributor's site-
                # specific GQ/DP), which was producing spuriously low GQ/DP on 0/0+0/1 sites.
                # GQ is taken from PL (see gq_from_pl) so the synthetic proband's displayed GQ in seqr
                # matches what seqr shows for the contributing parent.
                contrib_gqs: list[int] = []
                contrib_dps: list[int] = []
                if mother_contributor is not None:
                    contrib_gqs.append(gq_from_pl(mother_contributor))
                    contrib_dps.append(dp_of(mother_contributor))
                if father_contributor is not None:
                    contrib_gqs.append(gq_from_pl(father_contributor))
                    contrib_dps.append(dp_of(father_contributor))

                emit_variant(
                    out_vcf,
                    header,
                    sample_name,
                    contig,
                    seg_start,
                    common_ref,
                    maternal_alt,
                    paternal_alt,
                    min(contrib_gqs),
                    min(contrib_dps),
                    haploid=False,
                )

                # Consume any buffered starter records so they aren't re-processed later.
                # (The outer peek stays put and pops via the usual `m_end == seg_end` check below.)
                if mother_starter_result:
                    mother_iter.pop_at(mother_starter_result[0])
                if father_starter_result:
                    father_iter.pop_at(father_starter_result[0])
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

    # discard any single-parent tail left on this contig (synthetic proband is missing there)
    _skip_contig(mother_iter, contig)
    _skip_contig(father_iter, contig)


def create_synthetic_gvcf(mother_gvcf: str, father_gvcf: str, sample_name: str, out_path: str) -> None:
    """Generate the synthetic proband gVCF at out_path from the two parental gVCFs."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        local_mother = _localise(mother_gvcf, tmp)
        local_father = _localise(father_gvcf, tmp)
        # pysam forces INFO/END onto every record it writes (END is declared in the header), so we
        # write to a raw file first then strip END off variant lines in a text pass (see below).
        local_raw = str(tmp / 'raw.g.vcf.gz')
        local_out = str(tmp / 'synthetic.g.vcf.gz')

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
        logger.info(f'Wrote synthetic proband gVCF for {sample_name} to {out_path}')


def _strip_end_from_variants(raw_path: str, out_path: str, tmp: Path) -> None:
    """Remove INFO/END from variant records, keeping it only on <NON_REF>-only reference blocks.

    pysam writes INFO/END on every record (it is declared in the header), but the Hail combiner
    treats any record with a defined END as a reference block and errors if such a record has a
    non-hom-ref genotype. A record is a reference block iff its ALT is exactly <NON_REF>; every
    other record is a variant and must have END absent so the combiner reads it correctly.

    Streams bgzipped raw -> bgzipped output line by line; no plain-text intermediate. The previous
    design wrote an uncompressed VCF as an intermediate (~15-20x bgzip expansion), which on a WGS
    run exhausted the job's disk budget.
    """
    del tmp  # no scratch file needed with the streaming design
    with gzip.open(raw_path, 'rt') as fin, pysam.BGZFile(out_path, 'wb') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line.encode())
                continue
            columns = line.rstrip('\n').split('\t')
            alt, info = columns[4], columns[7]
            if alt == '<NON_REF>' or 'END=' not in info:
                fout.write(line.encode())
                continue
            kept = [field for field in info.split(';') if not field.startswith('END=')]
            columns[7] = ';'.join(kept) if kept else '.'
            fout.write(('\t'.join(columns) + '\n').encode())


def _next_contig(mother_iter: PeekIter, father_iter: PeekIter, contig_rank: dict[str, int]) -> str:
    """Pick the next contig to process - the lower-ranked of the two streams' current heads."""
    # Capture each peek exactly once so the None-filter narrows the value we actually use.
    heads = [head for head in (mother_iter.peek(), father_iter.peek()) if head is not None]
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
    """Clone every leading record on the given contig into the output, renamed to the synthetic proband sample."""
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
    parser = ArgumentParser(description='Synthesize a synthetic male proband gVCF from two parental gVCFs')
    parser.add_argument('--mother_gvcf', required=True, help='Path to the maternal gVCF (.g.vcf.gz)')
    parser.add_argument('--father_gvcf', required=True, help='Path to the paternal gVCF (.g.vcf.gz)')
    parser.add_argument('--output', required=True, help='Output synthetic gVCF path (.g.vcf.gz)')
    parser.add_argument('--sample_name', required=True, help='Sample (SM) name for the synthetic proband')
    args = parser.parse_args()
    create_synthetic_gvcf(
        mother_gvcf=args.mother_gvcf,
        father_gvcf=args.father_gvcf,
        sample_name=args.sample_name,
        out_path=args.output,
    )


if __name__ == '__main__':
    cli_main()
