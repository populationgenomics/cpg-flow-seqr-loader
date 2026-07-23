#!/usr/bin/env python3

"""
The entry point for the pre-first workflow, creating dummy proband gvcfs to go in to the combiner (first_workflow)
"""

import argparse

from cpg_flow import workflow

from cpg_seqr_loader.stages import GenerateDummyProbands, GenerateDummyProbandsCombinerInputs


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(
        name='seqr_loader',
        stages=[GenerateDummyProbands, GenerateDummyProbandsCombinerInputs],
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
