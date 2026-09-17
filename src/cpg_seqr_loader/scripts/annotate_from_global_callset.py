"""Annotate a densified cohort MT with row-level annotations from a global annotate_cohort.mt.

Used by workflows that skip the standard annotation stack (VEP / VQSR / gnomAD / clinvar joins)
because their cohort includes synthetic samples whose AC/AN/AF would pollute the global stats.
Instead of re-annotating from scratch, we join row-wise against the standard seqr-loader's
annotate_cohort.mt (produced from real samples only) and copy every row annotation across.

Some input_mt rows fail exact (locus, alleles) match against the global because Hail's
sparse_split_multi (in densify_VDS_to_MT) leaves non-minimal padding on split alleles, and the
input_mt and global went through independent combiner runs. Rather than shuffle the (very large)
global MT to re-normalise it, we recover those mismatches by driver-side locus-window search plus
allele trim - cost is bounded by the (small) mismatch count, not by global MT size. The recovered
rows are re-keyed to the global's non-minimal form so the annotation join lands.

Invariant: after key recovery, every input row must exist in the global. If not, either
create_synthetic_proband_gvcf.py has produced a variant no real sample carries, or temporal drift
has moved the parental gVCFs beyond what the global captured.
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

    # Rewrite input keys to the global's form for rows whose only mismatch is
    # trim-normalisation padding on the global side.
    input_mt = _recover_mismatched_keys_from_global(input_mt, global_rows)

    # After recovery, any remaining unmatched rows are real invariant violations.
    _assert_every_row_present_in_global(input_mt, global_rows)

    # Drop row annotations produced by densify - their AC/AN/AF etc. are computed over the
    # synthetic cohort and would be misleading. Global values replace them.
    input_mt = input_mt.drop('info', 'site_dp', 'ANS')

    loguru.logger.info('Joining global row annotations onto input MT')
    annotated_mt = input_mt.annotate_rows(**global_rows[input_mt.row_key])

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


def _recover_mismatched_keys_from_global(
    input_mt: hl.MatrixTable,
    global_rows: hl.Table,
    window_bp: int = 50,
) -> hl.MatrixTable:
    """Rewrite input_mt row keys that fail exact-match against global to the global's own key.

    For each input row not in global, search nearby global rows (locus +/- window_bp), trim their
    alleles common-bases-style, and check if any trimmed form matches the input row's key. If so,
    the input row is re-keyed to the global's (non-minimal) key so annotate_rows will land.

    The recovery walk is driver-side but bounded by the (small) mismatch count. Each per-mismatch
    lookup uses the global's row-key partition index, so cost is O(mismatches), not O(global size).
    """
    missing_rows = input_mt.rows().anti_join(global_rows).select().collect()
    if not missing_rows:
        loguru.logger.info('All input MT rows exact-match global - no key recovery needed')
        return input_mt

    loguru.logger.info(
        f'Attempting key recovery via locus-window search for {len(missing_rows)} unmatched rows',
    )

    rewrites: list[dict] = []
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
                break

    loguru.logger.info(
        f'Recovered keys for {len(rewrites)}/{len(missing_rows)} unmatched rows via trim-normalisation',
    )

    if not rewrites:
        return input_mt

    # Broadcast the rewrite table and use it to re-key the affected rows.
    rewrite_ht = hl.Table.parallelize(
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

    input_mt = input_mt.annotate_rows(
        key_rewrite_tmp=rewrite_ht[
            input_mt.locus.contig,
            input_mt.locus.position,
            input_mt.alleles[0],
            input_mt.alleles[1],
        ],
    )
    input_mt = input_mt.annotate_rows(
        new_locus_tmp=hl.if_else(
            hl.is_defined(input_mt.key_rewrite_tmp),
            hl.locus(
                input_mt.key_rewrite_tmp.new_contig,
                input_mt.key_rewrite_tmp.new_pos,
                reference_genome='GRCh38',
            ),
            input_mt.locus,
        ),
        new_alleles_tmp=hl.if_else(
            hl.is_defined(input_mt.key_rewrite_tmp),
            input_mt.key_rewrite_tmp.new_alleles,
            input_mt.alleles,
        ),
    )
    input_mt = input_mt.key_rows_by(locus=input_mt.new_locus_tmp, alleles=input_mt.new_alleles_tmp)
    return input_mt.drop('key_rewrite_tmp', 'new_locus_tmp', 'new_alleles_tmp')


def _assert_every_row_present_in_global(input_mt: hl.MatrixTable, global_rows: hl.Table) -> None:
    """Fail loud if the input MT (post key-recovery) has variants missing from the global callset.

    Uses head(10) as a cheap short-circuit - if there are no misses, we skip the full count().
    """
    missing = input_mt.rows().anti_join(global_rows)
    missing_examples = missing.head(10).collect()

    if not missing_examples:
        return

    missing_count = missing.count()
    example_strs = [f'{row.locus}:{row.alleles[0]}>{",".join(row.alleles[1:])}' for row in missing_examples]
    raise ValueError(
        f'AnnotateFromGlobalCallset invariant violation: {missing_count} variant rows in the '
        f'input MT do not exist in the global annotate_cohort.mt, even after key-recovery via '
        f'trim-normalisation. Either create_synthetic_proband_gvcf.py has produced variants no '
        f'real sample carries, or temporal drift has moved the parental gVCFs beyond what the '
        f'global captured. Investigate one of the missing variants in the parental gVCFs.\n'
        f'First {len(example_strs)} missing variants:\n  ' + '\n  '.join(example_strs),
    )


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
