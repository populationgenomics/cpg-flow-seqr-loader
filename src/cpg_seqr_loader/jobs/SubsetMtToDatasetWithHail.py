from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

from cpg_seqr_loader import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_subset_mt_job(
    dataset: targets.Dataset,
    input_mt: str,
    id_file: Path,
    output_mt: Path,
    job_attrs: dict[str, str],
) -> 'BashJob':
    family_sg_ids = family_sgs['family_sg_ids'] if (family_sgs := utils.get_family_sequencing_groups(dataset)) else None

    # write a list of all the SG IDs to retain
    # don't try and extract samples which didn't have a gVCF
    if not config.config_retrieve(['workflow', 'dry_run'], False):
        with id_file.open('w') as f:
            for sg in dataset.get_sequencing_groups():
                if ((family_sg_ids is None) or sg.id in family_sg_ids) and (sg.gvcf is not None):
                    f.write(f'{sg.id}\n')

    job = hail_batch.get_batch().new_bash_job(
        name=f'Subset MT for {dataset.name}',
        attributes=job_attrs,
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('highmem').storage('10Gi')
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.subset_mt_to_dataset \\
            --input {input_mt} \\
            --output {output_mt!s} \\
            --sg_id_file {id_file!s}
        """
    )
    return job
