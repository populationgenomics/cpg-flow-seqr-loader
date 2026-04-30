"""Normalize variants in a Hail MatrixTable to their minimal representation.

Strategy (all data stays in Hail — no collect to driver):
  1. Annotate every row with its normalized locus/alleles via min_rep
  2. Filter to only the rows whose representation changed
  3. Checkpoint that small subset to GCP so it can be reused cheaply
  4. anti_join_rows to remove the old unnormalized entries from the original MT
  5. Re-key the subset with normalized coordinates, dedup
  6. Cross-collision check against existing correct rows via semi_join
  7. union_rows + write
"""

import argparse

import loguru
from cpg_utils import config, hail_batch

import hail as hl


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Path to input MatrixTable')
    parser.add_argument('output', help='Path to write normalized MatrixTable')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output if it exists')

    args = parser.parse_args()

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        worker_cores=config.config_retrieve(['combiner', 'worker_cores']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),
    )

    mt = hl.read_matrix_table(args.input)

    mr = hl.min_rep(mt.locus, mt.alleles)
    mt_annotated = mt.annotate_rows(
        normal_locus=mr.locus,
        normal_alleles=mr.alleles,
    )

    # Filter to only rows whose representation changed
    mt_changed = mt_annotated.filter_rows(
        (mt_annotated.normal_locus != mt_annotated.locus) | (mt_annotated.normal_alleles != mt_annotated.alleles),
    )

    # Checkpoint
    changed_path = args.output.rstrip('/') + '_norm_changed.mt'
    mt_changed = mt_changed.checkpoint(changed_path, overwrite=True)

    n_changed = mt_changed.count_rows()
    loguru.logger.info(f'Found {n_changed} rows needing normalization')
    loguru.logger.info(f'Saved the small subset of changes to be made, {changed_path}')

    if n_changed == 0:
        loguru.logger.info('Done. No normalization needed. Happy Days')
        return

    # Remove unnormalized rows from original MT
    old_keys = mt_changed.rows().select().key_by('locus', 'alleles')
    mt_clean = mt.anti_join_rows(old_keys)

    # Re-key the changed subset with normalized coordinates
    mt_subset = mt_changed.key_rows_by(
        locus=mt_changed.normal_locus,
        alleles=mt_changed.normal_alleles,
    )
    mt_subset = mt_subset.drop('normal_locus', 'normal_alleles')
    mt_subset = mt_subset.distinct_by_row()

    # Check for collisions: do any new keys already exist in the clean data?
    new_keys = (
        mt_changed.rows()
        .key_by(
            locus=mt_changed.rows().normal_locus,
            alleles=mt_changed.rows().normal_alleles,
        )
        .select()
    )
    n_collisions = mt_clean.semi_join_rows(new_keys).count_rows()

    if n_collisions > 0:
        loguru.logger.info(f'{n_collisions} normalized keys collide with existing rows; dropping duplicates')
        collision_ht = mt_clean.semi_join_rows(new_keys).rows().select().key_by('locus', 'alleles')
        mt_subset = mt_subset.anti_join_rows(collision_ht)

    mt_clean.union_rows(mt_subset).write(args.output, overwrite=args.overwrite)
    loguru.logger.info(f'Done. Wrote normalized MatrixTable to {args.output}')


if __name__ == '__main__':
    main()
