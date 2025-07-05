from typing import TYPE_CHECKING

from cpg_flow import resources
from cpg_utils import config, hail_batch

from cpg_seqr_loader import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def train_vqsr_snp_model(
    sites_only_vcf: str,
    snp_model: str,
    job_attrs: dict,
) -> 'BashJob':
    """Train VQSR SNPs on the sites-only VCF."""

    local_res = utils.get_localised_resources_for_vqsr()
    siteonly_vcf = hail_batch.get_batch().read_input(sites_only_vcf)

    job = hail_batch.get_batch().new_job(
        'TrainVqsrSnpModelOnCombinerData',
        job_attrs | {'tool': 'gatk VariantRecalibrator'},
    )
    job.image(config.config_retrieve(['images', 'gatk']))

    # We run it for the entire dataset in one job, so can take an entire instance.
    res = resources.HIGHMEM.set_resources(
        j=job,
        fraction=1,
        storage_gb=utils.SNPS_RECAL_DISC_SIZE,
    )

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils.SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [f'-an {v}' for v in utils.SNP_ALLELE_SPECIFIC_FEATURES],
    )
    job.command(
        f"""set -euo pipefail
        gatk --java-options \
          "{res.java_mem_options()} {res.java_gc_thread_options()}" \\
          VariantRecalibrator \\
          -V {siteonly_vcf} \\
          -O {job.recalibration} \\
          --tranches-file {job.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          --use-allele-specific-annotations \\
          --sample-every-Nth-variant 10 \\
          --output-model {job.model_file} \\
          --max-gaussians 6 \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {local_res['hapmap'].base} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {local_res['omni'].base} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {local_res['one_thousand_genomes'].base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {local_res['dbsnp'].base}
          """,
    )

    hail_batch.get_batch().write_output(job.model_file, snp_model)
    return job
