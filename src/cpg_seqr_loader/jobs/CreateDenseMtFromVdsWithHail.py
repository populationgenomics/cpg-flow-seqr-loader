from typing import TYPE_CHECKING

from cpg_utils import hail_batch, config, Path
from cpg_flow import targets

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_densify_jobs(
        input_vds: str,
        output_mt: Path,
        output_sites_only: str,
        output_separate_header: str,
        job_attrs: dict[str, str],
) -> 'BashJob':

    job = hail_batch.get_batch().new_bash_job('Densify VDS and export MT', attributes=job_attrs | {'tool': 'hail'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.densify_VDS_to_MT \\
        --input {input_vds!s} \\
        --output {output_mt!s} \\
        --sites_only {output_sites_only!s} \\
        --separate_header {output_separate_header!s}
        """
    )

    return job

