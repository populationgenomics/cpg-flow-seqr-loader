from typing import TYPE_CHECKING

from cpg_utils import hail_batch, config, Path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_annotate_cohort_job(
    output_mt: Path,
    vep_ht: str,
    vqsr_vcf: str,
    variant_mt: str,
    checkpoint_prefix: Path,
    job_attrs: dict[str, str],
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job(
        'AnnotateCohort; join Vars, VEP, and VQSR',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(2).memory('highmem').storage('10Gi')
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.annotate_cohort \\
            --input {variant_mt} \\
            --output {output_mt!s} \\
            --vep {vep_ht!s} \\
            --checkpoint {checkpoint_prefix!s} \\
            --vqsr {vqsr_vcf}
        """
    )
    return job
