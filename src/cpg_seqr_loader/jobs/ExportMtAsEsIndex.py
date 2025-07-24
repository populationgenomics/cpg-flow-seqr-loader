from typing import TYPE_CHECKING

from cpg_flow import utils as cpg_flow_utils
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_es_export_job(
    index_name: str,
    done_flag: Path,
    mt_path: str,
    sample_size: int,
    job_attrs: dict[str, str],
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job(
        f'Generate {index_name} from {mt_path}',
        attributes=job_attrs | {'tool': 'elasticsearch'},
    )

    # prevent job failures when constructing HUGE file/directory names using Hail default behaviour
    job._dirname = f'{index_name}-{job._token}'  # noqa: SLF001

    # Use a non-preemptible instance if spot_instance is False in the config
    job = job.spot(is_spot=config.config_retrieve(['elasticsearch', 'spot_instance']))

    req_storage = cpg_flow_utils.tshirt_mt_sizing(
        sequencing_type=config.config_retrieve(['workflow', 'sequencing_type']),
        cohort_size=sample_size,
    )

    job.cpu(4).storage(f'{req_storage}Gi').memory('lowmem').image(config.config_retrieve(['workflow', 'driver_image']))

    # localise the MT
    job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

    # and just the name, used after localisation
    mt_name = mt_path.split('/')[-1]

    # run the export from the localised MT - this job writes no new data, just transforms and exports over network
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.mt_to_es \\
            --mt_path "${{BATCH_TMPDIR}}/{mt_name}" \\
            --index {index_name} \\
            --flag {done_flag!s}
        """
    )

    return job
