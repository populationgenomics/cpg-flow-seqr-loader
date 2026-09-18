"""Annotate a densified cohort MT with row-level annotations from a global annotate_cohort.mt.

Used by workflows that skip the standard annotation stack (VEP / VQSR / gnomAD / clinvar joins)
because their cohort includes synthetic samples whose AC/AN/AF would pollute the global stats.
Instead of re-annotating from scratch, we join row-wise against the standard seqr-loader's
annotate_cohort.mt (produced from real samples only) and copy every row annotation across.

Some input_mt rows fail exact (locus, alleles) match against the global because Hail's
sparse_split_multi (in densify_VDS_to_MT) leaves non-minimal padding on split alleles, and the
input_mt and global went through independent combiner runs. Rather than shuffle the (very large)
global MT to re-normalise it, we recover those mismatches by driver-side locus-window search plus
allele trim - cost is bounded by the (small) mismatch count, not by global MT size.

To avoid a full-MT distributed sort (which hangs on data skew for genome-scale MTs), we split the
input into:
  - "already_ok": rows whose exact keys already match global - joined directly, no re-key.
  - "needs_rewrite": tiny subset (~100s of rows) that need re-keying to the global's non-minimal
    form. Coalesced to one partition first, then re-keyed - trivial shuffle at that scale.
Then union_rows the two halves. Both are sorted, so union_rows does a linear merge.

Invariant: every input row must exist in the global (directly or via trim-recovery). If any row
can't be recovered, the recovery step raises loudly - either create_synthetic_proband_gvcf.py has
produced a variant no real sample carries, or temporal drift has moved the parental gVCFs beyond
what the global captured.
"""

import argparse

import loguru
from cpg_utils import hail_batch

import hail as hl


def annotate_from_global_callset(
    input_mt_path: str,
    global_mt_path: str,
    output_mt_path: str,
) -> None:
    """Row-join the input densified MT against a global annotate_cohort.mt and write the result."""
    loguru.logger.info(f'Reading input MT: {input_mt_path}')
    input_mt = hl.read_matrix_table(input_mt_path)

    loguru.logger.info(f'Reading global annotate_cohort MT: {global_mt_path}')
    global_mt = hl.read_matrix_table(global_mt_path)
    global_rows = global_mt.rows()

    # Drop row annotations produced by densify - their AC/AN/AF etc. are computed over the
    # synthetic cohort and would be misleading. Global values replace them.
    input_mt = input_mt.drop('info', 'site_dp', 'ANS')

    # Find rewrites for rows whose keys don't exact-match global. Raises loudly if any row is
    # unrecoverable (invented variant or temporal drift beyond the global's build).
    rewrite_ht = _build_rewrite_table(input_mt, global_rows)

    loguru.logger.info('Joining global row annotations onto input MT')
    annotated_mt = _annotate_via_split_union(input_mt, global_rows, rewrite_ht)

    global_metadata = global_mt.index_globals()
    annotated_mt = annotated_mt.annotate_globals(
        sourceFilePath=input_mt_path,
        genomeVersion=global_metadata.genomeVersion,
        sampleType=global_metadata.sampleType,
        hail_version=hl.__version__,
    )

    loguru.logger.info(f'Writing annotated MT: {output_mt_path}')
    annotated_mt.write(output_mt_path, overwrite=True)


def _trim_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Trim common right bases then common left (shifting position). Mirrors hl.min_rep semantics."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    shift = 0
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        shift += 1
    return pos + shift, ref, alt


def _build_rewrite_table(
    input_mt: hl.MatrixTable,
    global_rows: hl.Table,
    window_bp: int = 50,
) -> hl.Table | None:
    """Find rewrites for input_mt rows that fail exact-match against global.

    For each input row not in global, search nearby global rows (locus +/- window_bp), trim their
    alleles common-bases-style, and find one whose trimmed form matches the input row's key.
    Returns a keyed Table mapping (orig_contig, orig_pos, orig_ref, orig_alt) -> new (contig, pos,
    alleles). Returns None if no rows need rewriting.

    Raises ValueError if any row is unrecoverable - the workflow must not silently drop rows the
    user expected to appear in seqr.

    The recovery walk is driver-side but bounded by the (small) mismatch count. Each per-mismatch
    lookup uses the global's row-key partition index, so cost is O(mismatches), not O(global size).
    """
    missing_rows = input_mt.rows().anti_join(global_rows).select().collect()
    if not missing_rows:
        loguru.logger.info('All input MT rows exact-match global - no key recovery needed')
        return None

    loguru.logger.info(
        f'Attempting key recovery via locus-window search for {len(missing_rows)} unmatched rows',
    )

    rewrites: list[dict] = []
    unrecovered: list[tuple[str, int, str, str]] = []
    for row in missing_rows:
        contig = row.locus.contig
        pos = row.locus.position
        ref, alt = row.alleles[0], row.alleles[1]
        window = hl.locus_interval(
            contig,
            max(1, pos - window_bp),
            pos + window_bp,
            reference_genome='GRCh38',
        )
        candidates = global_rows.filter(window.contains(global_rows.locus)).select().collect()
        matched = False
        for cand in candidates:
            if cand.locus.contig != contig:
                continue
            cand_pos, cand_ref, cand_alt = _trim_alleles(
                cand.locus.position,
                cand.alleles[0],
                cand.alleles[1],
            )
            if cand_pos == pos and cand_ref == ref and cand_alt == alt:
                rewrites.append(
                    {
                        'orig_contig': contig,
                        'orig_pos': pos,
                        'orig_ref': ref,
                        'orig_alt': alt,
                        'new_contig': cand.locus.contig,
                        'new_pos': cand.locus.position,
                        'new_alleles': list(cand.alleles),
                    },
                )
                matched = True
                break
        if not matched:
            unrecovered.append((contig, pos, ref, alt))

    loguru.logger.info(
        f'Recovered keys for {len(rewrites)}/{len(missing_rows)} unmatched rows via trim-normalisation',
    )

    if unrecovered:
        example_strs = [f'{c}:{p} {r}>{a}' for c, p, r, a in unrecovered[:10]]
        raise ValueError(
            f'AnnotateFromGlobalCallset invariant violation: {len(unrecovered)} variant rows in the '
            f'input MT do not exist in the global annotate_cohort.mt, even after key-recovery via '
            f'trim-normalisation. Either create_synthetic_proband_gvcf.py has produced variants no '
            f'real sample carries, or temporal drift has moved the parental gVCFs beyond what the '
            f'global captured. Investigate one of the unrecovered variants in the parental gVCFs.\n'
            f'First {len(example_strs)} unrecovered variants:\n  ' + '\n  '.join(example_strs),
        )

    if not rewrites:
        return None

    return hl.Table.parallelize(
        rewrites,
        schema=hl.tstruct(
            orig_contig=hl.tstr,
            orig_pos=hl.tint32,
            orig_ref=hl.tstr,
            orig_alt=hl.tstr,
            new_contig=hl.tstr,
            new_pos=hl.tint32,
            new_alleles=hl.tarray(hl.tstr),
        ),
        key=['orig_contig', 'orig_pos', 'orig_ref', 'orig_alt'],
    )


def _annotate_via_split_union(
    input_mt: hl.MatrixTable,
    global_rows: hl.Table,
    rewrite_ht: hl.Table | None,
) -> hl.MatrixTable:
    """Join global row annotations onto input_mt without re-keying the whole MT.

    If no rewrites are needed, a direct annotate_rows suffices. Otherwise we split into two
    halves - one that matches global directly (kept at its original keys, no re-key) and one
    that needs re-keying (~100s of rows, tiny sort) - annotate each, then union_rows.

    Avoids the full-MT distributed sort that caused workers to hang on data skew.
    """
    if rewrite_ht is None:
        return input_mt.annotate_rows(**global_rows[input_mt.row_key])

    input_mt = input_mt.annotate_rows(
        rewrite_tmp=rewrite_ht[
            input_mt.locus.contig,
            input_mt.locus.position,
            input_mt.alleles[0],
            input_mt.alleles[1],
        ],
    )

    # Split into rows that match global directly vs rows that need re-keying.
    already_ok = input_mt.filter_rows(~hl.is_defined(input_mt.rewrite_tmp)).drop('rewrite_tmp')
    needs_rewrite = input_mt.filter_rows(hl.is_defined(input_mt.rewrite_tmp))

    # Direct join for the bulk of rows - no re-key, keys already match global.
    already_annotated = already_ok.annotate_rows(**global_rows[already_ok.row_key])

    # For the tiny subset that needs re-keying: coalesce to a single partition so the
    # subsequent sort works on a small compact input, then re-key to the global's form.
    # Bind naive_coalesce to a variable before referencing rewrite_tmp - chaining these
    # calls binds the field expressions to the pre-coalesce MT identity, and Hail refuses
    # to mix expressions from different-identity sources even when the schema matches.
    needs_rewrite = needs_rewrite.naive_coalesce(1)
    needs_rewrite = needs_rewrite.key_rows_by(
        locus=hl.locus(
            needs_rewrite.rewrite_tmp.new_contig,
            needs_rewrite.rewrite_tmp.new_pos,
            reference_genome='GRCh38',
        ),
        alleles=needs_rewrite.rewrite_tmp.new_alleles,
    ).drop('rewrite_tmp')
    newly_annotated = needs_rewrite.annotate_rows(**global_rows[needs_rewrite.row_key])

    # Both halves are sorted by (locus, alleles). union_rows does a linear merge over sorted
    # inputs - no distributed sort required.
    return already_annotated.union_rows(newly_annotated)


def cli_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to the densified input MT')
    parser.add_argument('--global_mt', required=True, help='Path to the global annotate_cohort.mt to join against')
    parser.add_argument('--output', required=True, help='Path to write the annotated MT')
    args = parser.parse_args()

    hail_batch.init_batch()
    annotate_from_global_callset(
        input_mt_path=args.input,
        global_mt_path=args.global_mt,
        output_mt_path=args.output,
    )


if __name__ == '__main__':
    cli_main()
