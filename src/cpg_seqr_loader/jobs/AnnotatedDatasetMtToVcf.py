from cpg_utils import Path, config, hail_batch


def cohort_to_vcf_job(
    input_mt: str,
    output_vcf: Path,
    job_attrs: dict,
):
    """Take the single-dataset MT, and write to a VCF."""

    vcf_j = hail_batch.get_batch().new_bash_job('VCF from dataset MT', job_attrs | {'tool': 'hail query'})
    vcf_j.image(config.config_retrieve(['workflow', 'driver_image']))
    vcf_j.command(f'python -m cpg_seqr_loader.scripts.vcf_from_mt_subset --input {input_mt} --output {output_vcf!s}')
    return vcf_j
