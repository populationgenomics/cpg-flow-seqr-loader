#!/usr/bin/env python3

"""
This is the entry point for the first workflow, running from gVCFs, through VDS, to a dense MatrixTable and VCFs
"""

import argparse

from cpg_flow import workflow

from cpg_seqr_loader.stages import CreateDenseMtFromVdsWithHail, DeleteCombinerTemp


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(name='seqr_loader', stages=[DeleteCombinerTemp, CreateDenseMtFromVdsWithHail], dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
