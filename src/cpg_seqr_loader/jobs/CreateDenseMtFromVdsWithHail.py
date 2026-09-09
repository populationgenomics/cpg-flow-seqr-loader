from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_densify_jobs(
    input_vds: str,
    output_mt: Path,
    output_sites_only: str | None,
    output_separate_header: str | None,
    checkpoint: str,
    job_attrs: dict[str, str],
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job('Densify VDS and export MT', attributes=job_attrs | {'tool': 'hail'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    job.spot(False)

    # Keep conditionally-empty args on the same line as a real arg, otherwise a
    # bare `--checkpoint <path> \` followed by a whitespace-only continuation
    # line breaks bash parsing when both optional args are absent.
    sites_only_arg = f'--sites_only {output_sites_only!s}' if output_sites_only else ''
    separate_header_arg = f'--separate_header {output_separate_header!s}' if output_separate_header else ''

    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.densify_VDS_to_MT \\
            --input {input_vds!s} \\
            --output {output_mt!s} \\
            --checkpoint {checkpoint!s} {sites_only_arg} {separate_header_arg}
        """
    )

    return job
