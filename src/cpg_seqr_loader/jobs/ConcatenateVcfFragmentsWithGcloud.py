from typing import TYPE_CHECKING

from cpg_utils import hail_batch, config, Path
from cpg_flow import targets

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_and_run_compose_script(
    multicohort: targets.MultiCohort,
    manifest_file: str,
    manifest_dir: str,
    output: Path,
    tmp_dir: Path,
    job_attrs: dict,
) -> list['BashJob']:
    local_manifest = hail_batch.get_batch().read_input(manifest_file)

    # generate a bash script to do the composition
    job_1 = hail_batch.get_batch().new_bash_job(f'Create Compose Script: {multicohort.name}', attributes=job_attrs)
    job_1.image(config.config_retrieve(['workflow', 'driver_image']))
    job_1.command(
        f"""
        python -m cpg_seqr_loader.scripts.write_gcloud_compose_script \\
        --input {local_manifest} \\
        --vcf_dir {manifest_dir} \\
        --output {output!s} \\
        --script {job_1.output} \\
        --tmp {tmp_dir / 'compose_intermediates' / multicohort.name!s}
        """,
    )

    job_2 = hail_batch.get_batch().new_bash_job(f'Run GCloud Compose: {multicohort.name}', attributes=job_attrs)
    job_2.image(config.config_retrieve(['images', 'cpg_hail_gcloud']))
    job_2.command(f'bash {job_1.output}')

    return [job_1, job_2]
