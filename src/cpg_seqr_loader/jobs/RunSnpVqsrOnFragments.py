from typing import TYPE_CHECKING

from cpg_flow import resources
from cpg_flow import utils as cpg_flow_utils
from cpg_utils import Path, config, hail_batch

from cpg_seqr_loader import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from hailtop.batch.resource import Resource, ResourceGroup


def quick_and_easy_bcftools_concat(
    input_vcfs: list['ResourceGroup | Resource'],
    storage_gb: int = 10,
    job_attrs: dict | None = None,
) -> 'BashJob':
    """
    A quick and easy way to concatenate VCFs
    Args:
        input_vcfs (list[hb.ResourceGroup]): list of VCFs to concatenate, requires the 'vcf.gz' attribute
        storage_gb (int): storage size for the job
        job_attrs (dict): job attributes, or None
    Returns:
        a ResourceGroup with the concatenated VCF, with the 'vcf.gz', and 'vcf.gz.tbi' attributes
    """
    job = hail_batch.get_batch().new_job(
        f'Concat {len(input_vcfs)} VCFs', (job_attrs or {}) | {'tool': 'bcftools concat'}
    )
    job.image(config.config_retrieve(['images', 'bcftools']))
    res = resources.STANDARD.set_resources(j=job, storage_gb=storage_gb)

    # declare a resource group for the concatenated output
    job.declare_resource_group(
        output={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
            'vcf.gz.csi': '{root}.vcf.gz.csi',
        },
    )
    job.command(
        f"""
    bcftools concat \\
        --threads {res.get_nthreads() - 1} \\
        -a {' '.join(vcf['vcf.gz'] for vcf in input_vcfs)} \\
        -Oz -o {job.output['vcf.gz']}
    tabix -p vcf {job.output['vcf.gz']}
    tabix -p vcf -C {job.output['vcf.gz']}
    """,
    )
    return job


def apply_snp_vqsr_to_fragments(
    manifest_file: Path,
    tranche_file: str,
    temp_path: Path,
    output_path: str,
    job_attrs: dict,
) -> list['BashJob']:
    """
    Apply SNP VQSR to the tranches
    I'm going to retry the stacking approach again to reduce job count

    Gather results into a single file

    Args:
        manifest_file (Path): path to the manifest file, locating all VCF fragments
        tranche_file ():
        temp_path (): Path to the temp from TrainVqsrSnpTranches
        output_path (str):
        job_attrs ():
    """

    vcf_resources = utils.get_all_fragments_from_manifest(manifest_file)
    fragment_count = len(vcf_resources)

    # read all the recal fragments into the batch as ResourceGroups
    # we're creating these paths in expectation that they were written by the tranches stage
    snps_recal_resources = [
        hail_batch.get_batch().read_input_group(
            recal=str(temp_path / f'snp_{i}.recal'),
            idx=str(temp_path / f'snp_{i}.recal.idx'),
        )
        for i in range(fragment_count)
    ]

    tranches_in_batch = hail_batch.get_batch().read_input(tranche_file)

    jobs = []
    recalibrated_vcfs = []

    vcf_counter = -1
    snp_filter_level = config.config_retrieve(['vqsr', 'snp_filter_level'])

    for chunk_counter, vcfs_recals in enumerate(
        cpg_flow_utils.generator_chunks(
            zip(vcf_resources, snps_recal_resources, strict=False), utils.RECALIBRATION_PER_JOB
        ),
    ):
        chunk_job = hail_batch.get_batch().new_bash_job(
            f'RunTrainedSnpVqsrOnCombinerFragments, Chunk {chunk_counter}', job_attrs
        )
        chunk_job.image(config.config_retrieve(['images', 'gatk']))

        # stores all the annotated VCFs in this chunk
        chunk_vcfs = []

        res = resources.STANDARD.set_resources(j=chunk_job, ncpu=1, storage_gb=10)

        # iterate over the zipped resource groups
        for vcf_resource, recal_resource in vcfs_recals:
            vcf_counter += 1
            # used in namespacing the outputs
            counter_string = str(vcf_counter)

            # create a resource group for the recalibration output and its index
            chunk_job.declare_resource_group(
                **{
                    counter_string: {
                        utils.VCF_GZ: '{root}.vcf.gz',
                        utils.VCF_GZ_TBI: '{root}.vcf.gz.tbi',
                    },
                },
            )

            chunk_job.command(
                f"""
            gatk --java-options "{res.java_mem_options()}" \\
                ApplyVQSR \\
                -O {chunk_job[counter_string]['vcf.gz']} \\
                -V {vcf_resource['vcf.gz']} \\
                --recal-file {recal_resource.recal} \\
                --tranches-file {tranches_in_batch} \\
                --truth-sensitivity-filter-level {snp_filter_level} \\
                --use-allele-specific-annotations \\
                -mode SNP
            tabix -p vcf -f {chunk_job[counter_string]['vcf.gz']}
            """,
            )
            chunk_vcfs.append(chunk_job[counter_string])

        # concatenates all VCFs in this chunk
        chunk_concat_job = quick_and_easy_bcftools_concat(
            chunk_vcfs,
            storage_gb=utils.SNPS_GATHER_DISC_SIZE,
            job_attrs=job_attrs,
        )
        chunk_concat_job.depends_on(chunk_job)
        jobs.extend([chunk_job, chunk_concat_job])
        recalibrated_vcfs.append(chunk_concat_job.output)

    # now we've got all the recalibrated VCFs, we need to gather them into a single VCF
    final_gather_job = quick_and_easy_bcftools_concat(
        recalibrated_vcfs,
        storage_gb=utils.SNPS_GATHER_DISC_SIZE,
        job_attrs=job_attrs,
    )
    final_gather_job.depends_on(*jobs)
    jobs.append(final_gather_job)

    hail_batch.get_batch().write_output(final_gather_job.output, output_path.removesuffix('.vcf.gz'))
    return jobs
