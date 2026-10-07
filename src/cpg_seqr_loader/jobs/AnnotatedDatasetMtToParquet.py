from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_mt_to_parquet_job(
    input_mt: str,
    json_dir: str,
    output_parquet: Path,
    job_attrs: dict,
) -> 'BashJob':
    """Take the single-dataset MT, and write to a long-format Parquet file."""

    job = hail_batch.get_batch().new_bash_job('Parquet from dataset MT', job_attrs | {'tool': 'hail query'})
    # Hail exports the MT as JSON shards, which are then streamed into Parquet in bounded blocks within this job
    job.cpu(4).memory('highmem').storage('10Gi')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.mt_to_parquet \\
            --input {input_mt} \\
            --json_dir {json_dir} \\
            --output {output_parquet!s}
        """
    )
    return job
