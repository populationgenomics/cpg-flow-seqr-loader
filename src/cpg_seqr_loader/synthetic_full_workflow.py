#!/usr/bin/env python3

"""End-to-end entry point for the synthetic-proband seqr-loader pipeline.

Runs the full DAG from synthetic gVCF invention through ES index load:
  GenerateSyntheticProbandGvcfs
    -> GenerateSyntheticProbandCombinerInputs
    -> CombineGvcfsIntoVdsFromManifest
    -> CreateDenseMtFromVdsWithHail (with VCF-fragment outputs disabled via config)
    -> AnnotateFromGlobalCallset
    -> SubsetMtToDatasetFromGlobalCallset (only when multi-dataset or only_families set)
    -> AnnotateDatasetFromGlobalCallset
    -> ExportMtAsEsIndexFromGlobalCallset

Only the terminal stage is listed - cpg-flow walks required_stages transitively to
build the full DAG.
"""

import argparse

from cpg_flow import workflow

from cpg_seqr_loader.synthetic_proband_stages import ExportMtAsEsIndexFromGlobalCallset


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(
        name='seqr_loader',
        stages=[ExportMtAsEsIndexFromGlobalCallset],
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
