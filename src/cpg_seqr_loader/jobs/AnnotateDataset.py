from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_annotate_dataset_job(
    dataset: targets.Dataset,
    input_mt: str,
    output_mt: Path,
    job_attrs: dict[str, str],
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job(
        f'AnnotateDataset for {dataset.name}',
        attributes=job_attrs | {'tool': 'hail'},
    )

    job.cpu(2).memory('highmem').storage('10Gi')

    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.annotate_dataset \\
            --input {input_mt} \\
            --output {output_mt!s}
        """
    )
    return job
