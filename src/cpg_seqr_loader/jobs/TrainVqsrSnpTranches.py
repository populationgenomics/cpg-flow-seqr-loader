from typing import TYPE_CHECKING

import loguru

from cpg_utils import hail_batch, config, Path
from cpg_flow import resources, utils as cpg_flow_utils
from cpg_seqr_loader import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from hailtop.batch.resource import ResourceFile


def train_vqsr_snp_tranches(
    manifest_file: Path,
    snp_model_path: str,
    output_path: str,
    temp_path: Path,
    job_attrs: dict,
) -> list['BashJob']:
    """
    Train vqsr tranches for SNPs, in a scattered manner.
    """

    vcf_resources = utils.get_all_fragments_from_manifest(manifest_file)

    fragment_count = len(vcf_resources)
    snps_recal_paths = [temp_path / f'snp_{i}.recal' for i in range(fragment_count)]
    snps_tranches_paths = [temp_path / f'snp_{i}.tranches' for i in range(fragment_count)]

    snp_model_in_batch = hail_batch.get_batch().read_input(snp_model_path)

    local_resources = utils.get_localised_resources_for_vqsr()

    # the list of all jobs (per-fragment, and the summary)
    scatter_jobs: list[BashJob] = []

    # to hold the resulting parts, as Resources
    snp_tranche_fragments: list[ResourceFile] = []

    # if we start this as -1, we can increment at the start of the loop, making the index counting easier to track
    vcf_counter = -1

    # iterate over all fragments, but in chunks of FRAGMENTS_PER_JOB
    for chunk_counter, fragment_chunk in enumerate(cpg_flow_utils.chunks(vcf_resources, utils.TRAINING_PER_JOB)):
        # NB this VQSR training stage is scattered, and we have a high number of very small VCF fragments
        # 99.5% of the time and cost of this task was pulling the docker image and loading reference data
        # the actual work took 5 seconds at negligible cost. Instead of running one job per VCF fragment,
        # we can batch fragments into fewer jobs, each stacking multiple fragments together but only requiring
        # one block of reference data to be loaded.
        # candidate splitting is 100-fragments-per-job, for a ~99% cost saving

        chunk_job = hail_batch.get_batch().new_job(f'TrainVqsrSnpTranches, Chunk {chunk_counter}', job_attrs)
        chunk_job.image(config.config_retrieve(['images', 'gatk']))
        chunk_job.command('set -euo pipefail')

        # add this job to the list of scatter jobs
        scatter_jobs.append(chunk_job)

        res = resources.STANDARD.set_resources(chunk_job, ncpu=4, storage_gb=utils.SNPS_GATHER_DISC_SIZE)

        # iterate over the fragment VCF resource groups
        for vcf_resource in fragment_chunk:
            vcf_counter += 1

            snps_recal_path = snps_recal_paths[vcf_counter]
            snps_tranche_path = snps_tranches_paths[vcf_counter]

            if cpg_flow_utils.can_reuse(snps_recal_path) and cpg_flow_utils.can_reuse(snps_tranche_path):
                loguru.logger.info(f'Reusing {snps_recal_path} and {snps_tranche_path}')
                snp_tranche_fragments.append(hail_batch.get_batch().read_input(str(snps_tranche_path)))
                continue

            tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils.SNP_RECALIBRATION_TRANCHE_VALUES])
            an_cmdl = ' '.join(
                [f'-an {v}' for v in utils.SNP_ALLELE_SPECIFIC_FEATURES],
            )

            # create a counter string to uniquely reference all job outputs
            counter_string = str(vcf_counter)

            # create a resource group for the recalibration output and its index
            chunk_job.declare_resource_group(
                **{
                    counter_string: {
                        'recal': '{root}.recal',
                        'recal.idx': '{root}.recal.idx',
                        'tranches': '{root}.tranches',
                    },
                },
            )

            # the mv command here is because the input VCF is a .bgz instead of .vcf.bgz
            # so GATK can't tell what type of file it is
            chunk_job.command(
                f"""
                MODEL_REPORT={snp_model_in_batch}
                mv {vcf_resource[utils.VCF_GZ]} input.vcf.bgz
                mv {vcf_resource[utils.VCF_GZ_TBI]} input.vcf.bgz.tbi
                gatk --java-options \
                  "{res.java_mem_options()} {res.java_gc_thread_options()}" \\
                  VariantRecalibrator \\
                  --verbosity WARNING \\
                  --QUIET \\
                  -V input.vcf.bgz \\
                  -O {chunk_job[counter_string]['recal']} \\
                  --tranches-file {chunk_job[counter_string]['tranches']} \\
                  --trust-all-polymorphic \\
                  {tranche_cmdl} \\
                  {an_cmdl} \\
                  -mode SNP \\
                  --use-allele-specific-annotations \\
                  --input-model {snp_model_in_batch} \\
                  --output-tranches-for-scatter \\
                  --max-gaussians 6 \\
                  -resource:hapmap,known=false,training=true,truth=true,prior=15 {local_resources['hapmap'].base} \\
                  -resource:omni,known=false,training=true,truth=true,prior=12 {local_resources['omni'].base} \\
                  -resource:1000G,known=false,training=true,truth=false,prior=10 {local_resources['one_thousand_genomes'].base} \\
                  -resource:dbsnp,known=true,training=false,truth=false,prior=7 {local_resources['dbsnp'].base}
                touch {chunk_job[counter_string]['recal.idx']}
                """,  # noqa: E501
            )

            # write the results out to GCP
            # TODO check this is OK without a str(output)
            hail_batch.get_batch().write_output(chunk_job[counter_string], temp_path / f'snp_{vcf_counter}')
            snp_tranche_fragments.append(chunk_job[counter_string].tranches)

    # one final job to write the success indicator
    final_job = hail_batch.get_batch().new_bash_job('Completion message', job_attrs)
    final_job.image(config.config_retrieve(['workflow', 'driver_image']))
    final_job.command(f'echo "All tranches trained" > {final_job.output}')
    final_job.depends_on(*scatter_jobs)
    scatter_jobs.append(final_job)
    hail_batch.get_batch().write_output(final_job.output, output_path)
    return scatter_jobs
