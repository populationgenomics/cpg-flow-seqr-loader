"""
suggested location for any utility methods or constants used across multiple stages
"""

from typing import TYPE_CHECKING
import datetime
import functools
import hashlib

import loguru
import hail as hl

from cpg_utils import config, hail_batch, Path
from metamist.graphql import gql, query
from cpg_flow import targets


if TYPE_CHECKING:
    from hailtop.batch.resource import ResourceGroup

DATE_STRING: str = datetime.datetime.now().strftime('%y-%m')  # noqa: DTZ005


TRAINING_PER_JOB: int = config.config_retrieve(['vqsr', 'vqsr_training_fragments_per_job'])
RECALIBRATION_PER_JOB: int = config.config_retrieve(['vqsr', 'vqsr_apply_fragments_per_job'])
INDEL_RECAL_DISC_SIZE: int = config.config_retrieve(['vqsr', 'indel_recal_disc_size'])
SNPS_RECAL_DISC_SIZE: int = config.config_retrieve(['vqsr', 'snps_recal_disc_size'])
SNPS_GATHER_DISC_SIZE: int = config.config_retrieve(['vqsr', 'snps_gather_disc_size'])

# some file extension constants
VCF_BGZ = 'vcf.bgz'
VCF_BGZ_TBI = 'vcf.bgz.tbi'
VCF_GZ = 'vcf.gz'
VCF_GZ_TBI = 'vcf.gz.tbi'

STANDARD_FEATURES = [
    'ReadPosRankSum',
    'MQRankSum',
    'QD',
    'FS',
    'SOR',
]
SNP_STANDARD_FEATURES = [*STANDARD_FEATURES, 'MQ']
INDEL_STANDARD_FEATURES = STANDARD_FEATURES

ALLELE_SPECIFIC_FEATURES = [
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
    'AS_FS',
    'AS_SOR',
    # Not using depth for the following reasons:
    # 1. The Broad pipelines don't use it;
    # 2. -G AS_StandardAnnotation flag to GenotypeGVCFs doesn't include it;
    # 3. For exomes, depth is an irrelevant feature and should be skipped:
    # 'AS_VarDP'
    # Note that for consistency, we also skip it for WGS.
]
SNP_ALLELE_SPECIFIC_FEATURES = [*ALLELE_SPECIFIC_FEATURES, 'AS_MQ']
INDEL_ALLELE_SPECIFIC_FEATURES = ALLELE_SPECIFIC_FEATURES

SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]


LATEST_ANALYSIS_QUERY = gql(
    """
    query LatestAnalysisEntry($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                meta
                output
                sequencingGroups {
                    id
                }
                timestampCompleted
            }
        }
    }
""",
)

SPECIFIC_VDS_QUERY = gql(
    """
    query getVDSByAnalysisId($vds_id: Int!) {
        analyses(id: {eq: $vds_id}) {
            output
            sequencingGroups {
                id
            }
        }
    }
""",
)


def query_for_specific_vds(vds_id: int) -> tuple[str, set[str]] | None:
    """
    query for a specific analysis of type entry_type for a dataset
    if found, return the set of SG IDs in the VDS (using the metadata)

    - stolen from the cpg_workflows.large_cohort.combiner Stage, but duplicated here so we can split pipelines without
      further code changes

    Args:
        vds_id (int): analysis id to query for

    Returns:
        either None if the analysis wasn't found, or a set of SG IDs in the VDS
    """

    # query for the exact, single analysis entry
    query_results: dict[str, dict] = query(SPECIFIC_VDS_QUERY, variables={'vds_id': vds_id})

    if not query_results['analyses']:
        return None
    vds_path: str = query_results['analyses'][0]['output']
    sg_ids = {sg['id'] for sg in query_results['analyses'][0]['sequencingGroups']}
    return vds_path, sg_ids


def query_for_latest_vds(dataset: str, entry_type: str = 'combiner') -> dict | None:
    """
    query for the latest analysis of type entry_type for a dataset
    Args:
        dataset (str): project to query for
        entry_type (str): type of analysis entry to query for
    Returns:
        str, the path to the latest analysis
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(LATEST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': entry_type})

    analyses_by_date = {}

    for analysis in result['project']['analyses']:
        if analysis['output'] and (
            analysis['meta']['sequencing_type'] == config.config_retrieve(['workflow', 'sequencing_type'])
        ):
            analyses_by_date[analysis['timestampCompleted']] = analysis

    if not analyses_by_date:
        loguru.logger.warning(f'No analysis of type {entry_type} found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort as strings
    return analyses_by_date[sorted(analyses_by_date)[-1]]


@functools.lru_cache(1)
def get_localised_resources_for_vqsr() -> dict[str, 'ResourceGroup']:
    """
    get the resources required for VQSR, once per run
    Returns:
        the dictionary of resources and their names
    """

    return {
        key: hail_batch.get_batch().read_input_group(
            base=config.reference_path(f'broad/{key}_vcf'),
            index=config.reference_path(f'broad/{key}_vcf_index'),
        )
        for key in [
            'axiom_poly',
            'dbsnp',
            'hapmap',
            'mills',
            'omni',
            'one_thousand_genomes',
        ]
    }


@functools.lru_cache(2)
def get_all_fragments_from_manifest(manifest_file: Path) -> list['ResourceGroup']:
    """
    read the manifest file, and return all the fragment resources as an ordered list
    this is a cached method as we don't want to localise every fragment once per task

    Args:
        manifest_file ():

    Returns:
        an ordered list of all the fragment VCFs and corresponding indices
    """

    resource_objects = []
    manifest_folder: Path = manifest_file.parent
    with manifest_file.open() as f:
        for line in f:
            vcf_path = manifest_folder / line.strip()
            resource_objects.append(
                hail_batch.get_batch().read_input_group(
                    **{
                        VCF_GZ: vcf_path,
                        VCF_GZ_TBI: f'{vcf_path}.tbi',
                    }
                ),
            )
    return resource_objects


@functools.cache
def get_family_sequencing_groups(dataset: targets.Dataset) -> dict | None:
    """
    Get the subset of sequencing groups that are in the specified families for a dataset
    Returns a dict containing the sequencing groups and a name suffix for the outputs
    """
    if not config.config_retrieve(['workflow', dataset.name, 'only_families'], []):
        return None
    only_family_ids = set(config.config_retrieve(['workflow', dataset.name, 'only_families'], []))
    # keep only the SG IDs for the families in the only_families list
    loguru.logger.info(f'Finding sequencing groups for families {only_family_ids} in dataset {dataset.name}')
    family_sg_ids = [sg.id for sg in dataset.get_sequencing_groups() if sg.pedigree.fam_id in only_family_ids]
    if not family_sg_ids:
        raise ValueError(f'No sequencing groups found for families {only_family_ids} in dataset {dataset.name}.')
    loguru.logger.info(f'Keeping only {len(family_sg_ids)} SGs from families {len(only_family_ids)} in {dataset}:')
    loguru.logger.info(only_family_ids)
    loguru.logger.info(family_sg_ids)

    h = hashlib.sha256(''.join(sorted(family_sg_ids)).encode()).hexdigest()[:4]
    name_suffix = f'{len(family_sg_ids)}_sgs-{len(only_family_ids)}_families-{h}'

    return {'family_sg_ids': family_sg_ids, 'name_suffix': name_suffix}


def manually_find_ids_from_vds(vds_path: str) -> set[str]:
    """
    during development and the transition to input_cohorts over input_datasets, there are some instances
    where we have VDS entries in Metamist, but the analysis entry contains SG IDs which weren't combined into the VDS

    This check bypasses the quick "get all SG IDs in the VDS analysis entry" check,
    and instead checks the exact contents of the VDS

    Args:
        vds_path (str): path to the VDS. Assuming it exists, this will be checked before calling this method

    Returns:
        set[str]: the set of sample IDs in the VDS
    """
    hail_batch.init_batch()
    vds = hl.vds.read_vds(vds_path)

    # find the samples in the Variant Data MT
    return set(vds.variant_data.s.collect())
