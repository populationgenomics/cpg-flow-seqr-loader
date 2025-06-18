from typing import TYPE_CHECKING

from cpg_utils import hail_batch, config
from cpg_seqr_loader import utils
from cpg_flow import resources

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def apply_recalibration_indels(
    snp_annotated_vcf: str,
    indel_recalibration: str,
    indel_tranches: str,
    output_path: str,
    job_attrs: dict,
) -> 'BashJob':
    """
    Apply indel recalibration to the annotated SNP VCF.
    """

    snp_vcf_in_batch = hail_batch.get_batch().read_input_group(
        vcf=snp_annotated_vcf,
        vcf_index=f'{snp_annotated_vcf}.tbi',
    )
    indel_tranches_in_batch = hail_batch.get_batch().read_input(indel_tranches)
    indel_recalibration_in_batch = hail_batch.get_batch().read_input_group(
        recal=indel_recalibration,
        recal_idx=f'{indel_recalibration}.idx',
    )

    job = hail_batch.get_batch().new_bash_job(f'RunTrainedIndelVqsrOnCombinedVcf on {snp_annotated_vcf}', job_attrs)
    job.image(config.config_retrieve(['images', 'gatk']))
    res = resources.STANDARD.set_resources(
        job,
        ncpu=2,
        storage_gb=utils.INDEL_RECAL_DISC_SIZE,
    )

    job.declare_resource_group(
        output={
            utils.VCF_GZ: '{root}.vcf.gz',
            utils.VCF_GZ_TBI: '{root}.vcf.gz.tbi',
        }
    )

    filter_level = config.config_retrieve(['vqsr', 'indel_filter_level'])

    job.command(
        f"""
    gatk --java-options "{res.java_mem_options()}" \\
        ApplyVQSR \\
        --tmp-dir $BATCH_TMPDIR \\
        -O {job.output[utils.VCF_GZ]} \\
        -V {snp_vcf_in_batch.vcf} \\
        --recal-file {indel_recalibration_in_batch.recal} \\
        --tranches-file {indel_tranches_in_batch} \\
        --truth-sensitivity-filter-level {filter_level} \\
        --use-allele-specific-annotations \\
        -mode INDEL
    tabix -p vcf -f {job.output[utils.VCF_GZ]}
    """,
    )
    hail_batch.get_batch().write_output(job.output, output_path.removesuffix('.vcf.gz'))
    return job
