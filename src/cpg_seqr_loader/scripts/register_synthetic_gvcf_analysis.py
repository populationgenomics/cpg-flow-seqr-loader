"""
Register a synthetic proband gVCF as an Analysis record in metamist, or refresh the existing
registration if it's stale.

Invoked as a Batch job by the GenerateSyntheticProbandGvcfs stage, after the corresponding
synthesis job has produced the gVCF (Batch dependency handles the ordering).

Idempotency + staleness policy
------------------------------
  - No active `synthetic_gvcf` Analysis exists for the two parental SG IDs  -> create one.
  - One exists whose meta source-gVCF paths match the current parental gVCFs  -> skip.
  - One exists but the recorded source paths differ from the current parental gVCFs  -> the
    existing Analysis is stale (a parent's gVCF was regenerated upstream), so deactivate the
    old record and create a fresh one.

The source-gVCF paths in `meta` are the mechanism for staleness detection - matching on SG IDs
alone would miss upstream gVCF regeneration since SG IDs stay stable when their gVCF is
re-produced.
"""

from argparse import ArgumentParser

from cpg_flow.metamist import get_metamist
from loguru import logger
from metamist import models
from metamist.apis import AnalysisApi
from metamist.graphql import gql, query

from cpg_seqr_loader.utils import SYNTHETIC_GVCF_ANALYSIS_TYPE


# GraphQL: list every active synthetic_gvcf Analysis in the project. We filter to the specific
# parental SG pair in Python (a `sequencingGroupIds` filter that requires exact-set match is not
# directly supported by the API - see AnalysisQueryModel.sequencing_group_ids, which is an
# "any-of" contains-filter, not "equals-set").
EXISTING_SYNTHETIC_ANALYSES_QUERY = gql(
    """
    query ExistingSyntheticAnalyses($project: String!, $analysis_type: String!) {
        project(name: $project) {
            analyses(active: {eq: true}, type: {eq: $analysis_type}) {
                id
                output
                meta
                sequencingGroups {
                    id
                }
            }
        }
    }
    """,
)


def find_existing_registration(project: str, mother_sg_id: str, father_sg_id: str) -> dict | None:
    """Return the active synthetic_gvcf Analysis linked to exactly these two SGs, or None."""
    result = query(
        EXISTING_SYNTHETIC_ANALYSES_QUERY,
        variables={'project': project, 'analysis_type': SYNTHETIC_GVCF_ANALYSIS_TYPE},
    )
    target = {mother_sg_id, father_sg_id}

    matches = [
        analysis
        for analysis in result['project']['analyses']
        if {sg['id'] for sg in analysis['sequencingGroups']} == target
    ]

    if not matches:
        return None
    if len(matches) > 1:
        # Shouldn't happen if the workflow always deactivates old records, but log defensively.
        logger.warning(
            f'Found {len(matches)} active synthetic_gvcf analyses for SGs {sorted(target)}; '
            f'using the last one (id={matches[-1]["id"]}) and leaving the others alone',
        )
    return matches[-1]


def deactivate_analysis(analysis_id: int) -> None:
    """Mark an existing Analysis as inactive (soft-delete)."""
    AnalysisApi().update_analysis(
        analysis_id=analysis_id,
        analysis_update_model=models.AnalysisUpdateModel(active=False),
    )
    logger.info(f'Marked Analysis id={analysis_id} inactive (stale source gVCFs)')


def create_registration(
    project: str,
    gvcf_path: str,
    sample_name: str,
    family_id: str,
    mother_sg_id: str,
    father_sg_id: str,
    mother_source_gvcf: str,
    father_source_gvcf: str,
) -> int | None:
    """Create a fresh synthetic_gvcf Analysis. Returns the new analysis id.

    Uses cpg_flow.metamist.get_metamist().create_analysis so that the metamist project name is
    adjusted for the current access level (bare `ravenscroft-rpl` becomes `ravenscroft-rpl-test`
    at test access) and retries are handled by the framework.
    """
    return get_metamist().create_analysis(
        output=gvcf_path,
        type_=SYNTHETIC_GVCF_ANALYSIS_TYPE,
        status='completed',
        sequencing_group_ids=[mother_sg_id, father_sg_id],
        dataset=project,
        meta={
            'family_id': family_id,
            'sample_name_in_gvcf': sample_name,
            'source_mother_gvcf': mother_source_gvcf,
            'source_father_gvcf': father_source_gvcf,
        },
    )


def cli_main():
    parser = ArgumentParser(description='Register or refresh a synthetic proband gVCF in metamist')
    parser.add_argument('--project', required=True, help='Metamist project (e.g. ravenscroft-rpl)')
    parser.add_argument('--gvcf_path', required=True, help='gs:// path to the synthetic gVCF')
    parser.add_argument('--sample_name', required=True, help='Sample name embedded in the gVCF header')
    parser.add_argument('--family_id', required=True, help='External family ID')
    parser.add_argument('--mother_sg_id', required=True, help='Mother sequencing_group id')
    parser.add_argument('--father_sg_id', required=True, help='Father sequencing_group id')
    parser.add_argument('--mother_source_gvcf', required=True, help='Path to the mother gVCF used as input')
    parser.add_argument('--father_source_gvcf', required=True, help='Path to the father gVCF used as input')
    args = parser.parse_args()

    existing = find_existing_registration(
        project=args.project,
        mother_sg_id=args.mother_sg_id,
        father_sg_id=args.father_sg_id,
    )

    if existing is None:
        create_registration(
            project=args.project,
            gvcf_path=args.gvcf_path,
            sample_name=args.sample_name,
            family_id=args.family_id,
            mother_sg_id=args.mother_sg_id,
            father_sg_id=args.father_sg_id,
            mother_source_gvcf=args.mother_source_gvcf,
            father_source_gvcf=args.father_source_gvcf,
        )
        return

    existing_meta = existing.get('meta') or {}
    if (
        existing_meta.get('source_mother_gvcf') == args.mother_source_gvcf
        and existing_meta.get('source_father_gvcf') == args.father_source_gvcf
    ):
        logger.info(
            f'Analysis id={existing["id"]} for family {args.family_id} is already registered with '
            f'matching source gVCFs; nothing to do.',
        )
        return

    logger.info(
        f'Analysis id={existing["id"]} for family {args.family_id} has stale source gVCF paths; '
        f'deactivating and re-registering.',
    )
    deactivate_analysis(existing['id'])
    create_registration(
        project=args.project,
        gvcf_path=args.gvcf_path,
        sample_name=args.sample_name,
        family_id=args.family_id,
        mother_sg_id=args.mother_sg_id,
        father_sg_id=args.father_sg_id,
        mother_source_gvcf=args.mother_source_gvcf,
        father_source_gvcf=args.father_source_gvcf,
    )


if __name__ == '__main__':
    cli_main()
