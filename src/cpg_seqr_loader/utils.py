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

from cpg_flow import targets


if TYPE_CHECKING:
    from hailtop.batch.resource import ResourceGroup

DATE_STRING: str = datetime.datetime.now().strftime('%y-%m')  # noqa: DTZ005


TRAINING_PER_JOB: int = config.config_retrieve(['rd_combiner', 'vqsr_training_fragments_per_job'], 100)
RECALIBRATION_PER_JOB: int = config.config_retrieve(['rd_combiner', 'vqsr_apply_fragments_per_job'], 60)
INDEL_RECAL_DISC_SIZE: int = config.config_retrieve(['rd_combiner', 'indel_recal_disc_size'], 20)
SNPS_RECAL_DISC_SIZE: int = config.config_retrieve(['rd_combiner', 'snps_recal_disc_size'], 20)
SNPS_GATHER_DISC_SIZE: int = config.config_retrieve(['rd_combiner', 'snps_gather_disc_size'], 10)

# some file extension constants
VCF_BGZ = 'vcf.bgz'
VCF_BGZ_TBI = 'vcf.bgz.tbi'
VCF_GZ = 'vcf.gz'
VCF_GZ_TBI = 'vcf.gz.tbi'


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
