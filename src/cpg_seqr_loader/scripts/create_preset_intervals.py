#!/usr/bin/env python3

"""
Script to read in the gnomAD 4.1 data table, use it to obtain evenly distributed (variant-balanced) intervals

n.b. requires an additional Mito interval... (does it? we should produce with/without, we have a mito calling workflow)

requires some additional re-splitting of the resulting intervals e.g. to have clean breaks at contig ends

target interval counts - 50, 100, 500, 1000, 2000

some inspo here
- https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/scripts/generate_new_intervals.py
- https://github.com/populationgenomics/production-pipelines/pull/883
- meres: gs://cpg-common-main/references/hg38/v0/hg38.telomeresAndMergedCentromeres.interval_list
"""

from argparse import ArgumentParser

from cpg_utils import Path, to_path, hail_batch


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--gnomad', required=True, type=str, help='Path to the gnomAD HT.')
    parser.add_argument('--count', nargs='+', type=int, required=True, help='Interval counts to generate.', )
    parser.add_argument('--output', required=True, type=str, help='Dir to write output files.', )
    parser.add_argument('--meres', default='gs://cpg-common-main/references/hg38/v0/hg38.telomeresAndMergedCentromeres.interval_list', type=str, help='Path to centro/telomere regions file.',
                        )
    args = parser.parse_args()
    print(args)
