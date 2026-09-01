"""
Register a synthetic proband gVCF as an Analysis record in metamist.

Invoked as a Batch job by the GenerateSyntheticProbandGvcfs stage, after the corresponding
synthesis job has produced the gVCF (Batch dependency handles the ordering).

Uses `cpg_utils.metamist_registration.create_new`, which is the canonical entry point for
metamist analysis registration in this codebase (see also cpg-flow-align-genotype's
sg_qc_report.py). It writes analysis records with the modern `outputs` block, enforces that the
primary output is a `gs://` path that actually exists, and lets cpg-utils own the GraphQL call.

Idempotency is handled entirely by cpg-flow's REUSE mechanism, not this script: the Stage 1
`_registered.txt` sentinel is what stops cpg-flow from queuing this script on subsequent runs
once a family has been registered. This script itself is deliberately not idempotent - matching
the rest of the codebase's `analysis_type=...`-based stages, which also call the metamist
`create` primitive unconditionally when queued.
"""

from argparse import ArgumentParser

from cpg_utils import to_path
from cpg_utils.config import dataset_for_access_level
from cpg_utils.metamist_registration import create_new
from loguru import logger

from cpg_seqr_loader.utils import SYNTHETIC_GVCF_ANALYSIS_TYPE


def cli_main():
    parser = ArgumentParser(description='Register a synthetic proband gVCF in metamist')
    parser.add_argument('--project', required=True, help='Metamist project (base name, no -test suffix)')
    parser.add_argument('--gvcf_path', required=True, help='gs:// path to the synthetic gVCF')
    parser.add_argument('--sample_name', required=True, help='Sample name embedded in the gVCF header')
    parser.add_argument('--family_id', required=True, help='External family ID')
    parser.add_argument('--mother_sg_id', required=True, help='Mother sequencing_group id')
    parser.add_argument('--father_sg_id', required=True, help='Father sequencing_group id')
    parser.add_argument('--mother_source_gvcf', required=True, help='Path to the mother gVCF used as input')
    parser.add_argument('--father_source_gvcf', required=True, help='Path to the father gVCF used as input')
    parser.add_argument(
        '--marker_path',
        required=True,
        help=(
            'Path (typically gs://) to a per-family sentinel file that is written when this '
            'invocation completes successfully. Its presence is what cpg_flow uses to decide '
            'whether Stage 1 can be reused on subsequent runs - without it, the stage-wide REUSE '
            'check would skip queue_jobs (and thus registration) whenever the gVCF already exists.'
        ),
    )
    args = parser.parse_args()

    create_new(
        project=dataset_for_access_level(args.project),
        output=args.gvcf_path,
        analysis_type=SYNTHETIC_GVCF_ANALYSIS_TYPE,
        sgs=[args.mother_sg_id, args.father_sg_id],
        meta={
            'family_id': args.family_id,
            'sample_name_in_gvcf': args.sample_name,
            'source_mother_gvcf': args.mother_source_gvcf,
            'source_father_gvcf': args.father_source_gvcf,
        },
    )

    # Sentinel: write only after registration succeeds (any exception above skips this).
    # The empty file's existence is what unblocks cpg_flow's stage-level REUSE on future runs.
    with to_path(args.marker_path).open('w') as marker:
        marker.write('')
    logger.info(f'Wrote registration marker to {args.marker_path}')


if __name__ == '__main__':
    cli_main()
