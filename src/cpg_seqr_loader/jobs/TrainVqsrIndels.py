from typing import TYPE_CHECKING

from cpg_flow import resources
from cpg_utils import config, hail_batch

from cpg_seqr_loader import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def train_vqsr_indel_model(
    sites_only_vcf: str,
    output_prefix: str,
    job_attrs: dict,
) -> 'BashJob':
    """Train VQSR indels on the sites-only VCF."""
    local_resources = utils.get_localised_resources_for_vqsr()
    siteonly_vcf = hail_batch.get_batch().read_input(sites_only_vcf)

    """
    Run VariantRecalibrator to calculate VQSLOD tranches for indels

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value. 4 is a
    reasonable default for indels, as their number is smaller than SNPs.
    """
    job = hail_batch.get_batch().new_job(
        'TrainVqsrIndelModelOnCombinerData',
        job_attrs | {'tool': 'gatk VariantRecalibrator'},
    )
    job.image(config.config_retrieve(['images', 'gatk']))
    job.command('set -euo pipefail')

    # We run it for the entire dataset in one job, so can take an entire instance.
    res = resources.HIGHMEM.set_resources(j=job, fraction=1, storage_gb=utils.INDEL_RECAL_DISC_SIZE)

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils.INDEL_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [f'-an {v}' for v in utils.INDEL_ALLELE_SPECIFIC_FEATURES],
    )

    # delclare a resource group for the output
    job.declare_resource_group(
        output={
            'recal': '{root}.recal',
            'recal.idx': '{root}.recal.idx',
            'tranches': '{root}.tranches',
        },
    )
    job.command(
        f"""
        tabix {siteonly_vcf}
        gatk --java-options \
          "{res.java_mem_options()} {res.java_gc_thread_options()}" \\
          VariantRecalibrator \\
          -V {siteonly_vcf} \\
          -O {job.output.recal} \\
          --tranches-file {job.output.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          --use-allele-specific-annotations \\
          --max-gaussians 4 \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {local_resources['mills'].base} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {local_resources['axiom_poly'].base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {local_resources['dbsnp'].base}
        """,
    )
    hail_batch.get_batch().write_output(job.output, output_prefix)
    return job
