#!/usr/bin/env python3

"""
This is the main entry point for the secondary workflow, running on VCF fragments, and resulting in per-dataset MTs.
"""

import argparse

from cpg_flow import workflow

from cpg_seqr_loader.stages import AnnotatedDatasetMtToVcf, ExportMtAsEsIndex


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(stages=[ExportMtAsEsIndex, AnnotatedDatasetMtToVcf], dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
