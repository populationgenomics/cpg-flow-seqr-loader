#!/usr/bin/env python3
"""Normalize variants in a Hail MatrixTable to their minimal representation.

Strategy (avoids full-table shuffle and doubled DAG):
  1. Single scan to count unnormalized rows, then collect their old+new coordinates
  2. Build tiny Hail Tables from the collected Python data
  3. anti_join_rows with tiny old-keys table to remove unnormalized entries (broadcast join)
  4. filter_intervals + semi_join to extract the subset (partition-pruned), re-key, dedup
  5. Cross-collision check against existing correct rows (partition-pruned)
  6. union_rows + write
"""

import argparse
import hail as hl
from cpg_utils import config, hail_batch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Path to input MatrixTable')
    parser.add_argument('output', help='Path to write normalized MatrixTable')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output if it exists')

    args = parser.parse_args()

    hl.init_batch(worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        worker_cores=config.config_retrieve(['combiner', 'worker_cores']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),)

    mt = hl.read_matrix_table(args.input)


    mr = hl.min_rep(mt.locus, mt.alleles)
    needs_norm = (mr.locus != mt.locus) | (mr.alleles != mt.alleles)

    collected = mt.aggregate_rows(
        hl.agg.filter(
            needs_norm,
            hl.agg.collect(
                hl.struct(
                    old_locus=mt.locus,
                    old_alleles=mt.alleles,
                    new_locus=mr.locus,
                    new_alleles=mr.alleles,
                )
            ),
        )
    )
    n_changed = len(collected)
    print(f'Found {n_changed} rows needing normalization')

    if n_changed == 0:
        mt.write(args.output, overwrite=args.overwrite)
        print(f'Done. No normalization needed. Wrote to {args.output}')
        return


    locus_schema = hl.tstruct(locus=hl.tlocus('GRCh38'), alleles=hl.tarray(hl.tstr))

    old_keys_ht = hl.Table.parallelize(
        [hl.Struct(locus=c.old_locus, alleles=c.old_alleles) for c in collected],
        schema=locus_schema,
    ).key_by('locus', 'alleles')

    new_keys_ht = hl.Table.parallelize(
        [hl.Struct(locus=c.new_locus, alleles=c.new_alleles) for c in collected],
        schema=locus_schema,
    ).key_by('locus', 'alleles')


    mt_clean = mt.anti_join_rows(old_keys_ht)


    old_intervals = [
        hl.parse_locus_interval(f'[{loc}-{loc}]', reference_genome='GRCh38')
        for loc in {str(c.old_locus) for c in collected}
    ]
    mt_subset = hl.filter_intervals(mt, old_intervals)
    mt_subset = mt_subset.semi_join_rows(old_keys_ht)

    mr_sub = hl.min_rep(mt_subset.locus, mt_subset.alleles)
    mt_subset = mt_subset.key_rows_by(locus=mr_sub.locus, alleles=mr_sub.alleles)
    mt_subset = mt_subset.distinct_by_row()


    new_intervals = [
        hl.parse_locus_interval(f'[{loc}-{loc}]', reference_genome='GRCh38')
        for loc in {str(c.new_locus) for c in collected}
    ]
    mt_collisions = hl.filter_intervals(mt_clean, new_intervals)
    mt_collisions = mt_collisions.semi_join_rows(new_keys_ht)
    n_collisions = mt_collisions.count_rows()

    if n_collisions > 0:
        print(f'{n_collisions} normalized keys collide with existing rows; dropping duplicates')
        collision_ht = mt_collisions.rows().select().key_by('locus', 'alleles')
        mt_subset = mt_subset.anti_join_rows(collision_ht)


    mt_clean.union_rows(mt_subset).write(args.output, overwrite=args.overwrite)
    print(f'Done. Wrote normalized MatrixTable to {args.output}')


if __name__ == '__main__':
    main()
