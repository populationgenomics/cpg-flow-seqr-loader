from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_annotate_from_global_callset_job(
    input_mt: str,
    global_mt: str,
    output_mt: Path,
    checkpoint_path: str,
    job_attrs: dict[str, str],
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job(
        'AnnotateFromGlobalCallset; join cohort MT against global annotate_cohort.mt',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('highmem').storage('10Gi')
    job.spot(False)
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.annotate_from_global_callset \\
            --input {input_mt!s} \\
            --global_mt {global_mt!s} \\
            --output {output_mt!s} \\
            --checkpoint {checkpoint_path!s}
        """
    )
    return job
