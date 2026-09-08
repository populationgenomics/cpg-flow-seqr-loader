"""Annotate a densified cohort MT with row-level annotations from a global annotate_cohort.mt.

Used by workflows that skip the standard annotation stack (VEP / VQSR / gnomAD / clinvar joins)
because their cohort includes synthetic samples whose AC/AN/AF would pollute the global stats.
Instead of re-annotating from scratch, we join row-wise against the standard seqr-loader's
annotate_cohort.mt (produced from real samples only) and copy every row annotation across.

Invariant: every row in the input MT must exist in the global MT. If not, the synthetic-proband
gVCF logic (create_synthetic_proband_gvcf.py) has produced a variant that no real sample carries,
which shouldn't be possible under its worst-case-inheritance rules.
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


def _assert_every_row_present_in_global(input_mt: hl.MatrixTable, global_rows: hl.Table) -> None:
    """Fail loud if the input MT has variants missing from the global callset.

    Uses head(10) as a cheap short-circuit — if there are no misses, we skip the full count().
    """
    missing = input_mt.rows().anti_join(global_rows)
    missing_examples = missing.head(10).collect()

    if not missing_examples:
        return

    missing_count = missing.count()
    example_strs = [f'{row.locus}:{row.alleles[0]}>{",".join(row.alleles[1:])}' for row in missing_examples]
    raise ValueError(
        f'AnnotateFromGlobalCallset invariant violation: {missing_count} variant rows in the '
        f'input MT do not exist in the global annotate_cohort.mt.\n'
        f'Every input variant must also appear in the global callset (via the real parental '
        f'samples). If this fires, create_synthetic_proband_gvcf.py has produced variants not '
        f'present in any real sample - investigate that script.\n'
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
