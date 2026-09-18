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

The recovered mismatches are handled via an "effective lookup key" computed per row: for the
~124 recovered rows the global's non-minimal (locus, alleles) is used to fetch annotations; for
every other row the original key is used. Global annotations still land on every input row and
the input MT is never re-keyed.

Consequence: the output MT keeps min-rep keys throughout, including for the 124 recovered rows.
Those variants appear in seqr under min-rep form rather than the global's non-minimal padding.
This is arguably the more canonical representation (matches ClinVar / gnomAD) and no variants
are lost.

Before the join we checkpoint the input MT (with the added lookup-key fields) so QoB doesn't
have to fuse the entire pipeline into one big shuffle and so partial-failure recovery is cheap.
QoB resourcing (driver_cores / worker_cores) is exposed via config keys under
`annotate_from_global_callset` for the retry ladder (start 2/1, then 2/2, then 4/2).

Invariant: every input row must exist in the global (directly or via trim-recovery). If any row
can't be recovered, the recovery step raises loudly - either create_synthetic_proband_gvcf.py has
produced a variant no real sample carries, or temporal drift has moved the parental gVCFs beyond
what the global captured.
"""

import argparse

import loguru
from cpg_utils import config, hail_batch

import hail as hl


def annotate_from_global_callset(
    input_mt_path: str,
    global_mt_path: str,
    output_mt_path: str,
    checkpoint_path: str,
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
    annotated_mt = _annotate_via_effective_key(input_mt, global_rows, rewrite_ht, checkpoint_path)

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


def _annotate_via_effective_key(
    input_mt: hl.MatrixTable,
    global_rows: hl.Table,
    rewrite_ht: hl.Table | None,
    checkpoint_path: str,
) -> hl.MatrixTable:
    """Join global row annotations onto input_mt without re-keying it.

    Computes an "effective lookup key" per row (rewrite -> global's non-minimal (locus,
    alleles) for the ~124 rewritten rows; original key for everything else) and indexes
    global_rows by that computed key. Output MT keeps the input's min-rep keys throughout.

    Before the final annotation join, we checkpoint the input MT with its lookup-key fields
    added. Two reasons (per Ed): (1) makes recovery cheap if the join step fails, (2) splits
    Hail's query graph so QoB doesn't have to fuse the entire pipeline into one shuffle.
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
    input_mt = input_mt.annotate_rows(
        lookup_locus_tmp=hl.if_else(
            hl.is_defined(input_mt.rewrite_tmp),
            hl.locus(
                input_mt.rewrite_tmp.new_contig,
                input_mt.rewrite_tmp.new_pos,
                reference_genome='GRCh38',
            ),
            input_mt.locus,
        ),
        lookup_alleles_tmp=hl.if_else(
            hl.is_defined(input_mt.rewrite_tmp),
            input_mt.rewrite_tmp.new_alleles,
            input_mt.alleles,
        ),
    )

    loguru.logger.info(f'Checkpointing input MT with lookup keys to {checkpoint_path}')
    input_mt = input_mt.checkpoint(checkpoint_path, overwrite=True)

    # Index global_rows by the computed lookup key. Hail evaluates global_rows[<expr>, <expr>]
    # as a keyed lookup and delivers the matching row struct.
    annotated_mt = input_mt.annotate_rows(
        **global_rows[input_mt.lookup_locus_tmp, input_mt.lookup_alleles_tmp],
    )
    return annotated_mt.drop('rewrite_tmp', 'lookup_locus_tmp', 'lookup_alleles_tmp')


def cli_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to the densified input MT')
    parser.add_argument('--global_mt', required=True, help='Path to the global annotate_cohort.mt to join against')
    parser.add_argument('--output', required=True, help='Path to write the annotated MT')
    parser.add_argument('--checkpoint', required=True, help='Path to checkpoint the input MT before the join')
    args = parser.parse_args()

    # driver_cores + worker_cores tune QoB resourcing for the shuffle-heavy join step.
    # Defaults per Ed's guidance for large-MT joins; override via config for retry ladder.
    hail_batch.init_batch(
        driver_cores=config.config_retrieve(['annotate_from_global_callset', 'driver_cores'], 2),
        worker_cores=config.config_retrieve(['annotate_from_global_callset', 'worker_cores'], 1),
    )
    annotate_from_global_callset(
        input_mt_path=args.input,
        global_mt_path=args.global_mt,
        output_mt_path=args.output,
        checkpoint_path=args.checkpoint,
    )


if __name__ == '__main__':
    cli_main()
