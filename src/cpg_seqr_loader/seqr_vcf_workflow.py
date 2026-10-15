#!/usr/bin/env python3

"""
This is the entry point for the seqr VCF workflow, the alternative second phase following first_workflow.
Runs VQSR on the VCF fragments, and writes per-dataset seqr VCFs from the dense MT, skipping VEP and all annotation.
"""

import argparse

from cpg_flow import workflow

from cpg_seqr_loader.stages import SeqrVcfFromDenseMt


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(
        name='seqr_loader',
        stages=[SeqrVcfFromDenseMt],
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
