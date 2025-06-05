"""
suggested location for any utility methods or constants used across multiple stages
"""

from typing import TYPE_CHECKING
import functools
import datetime

from cpg_utils import config, hail_batch, Path

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
