"""
All Stages relating to the seqr_loader pipeline, reimplemented from scratch to
use the gVCF combiner instead of joint-calling.
"""

from google.api_core.exceptions import PermissionDenied

import loguru

import cpg_flow

from cpg_seqr_loader import jobs, utils

from cpg_utils import Path, config, cloud

from metamist.graphql import gql, query

LATEST_ANALYSIS_QUERY = gql(
    """
    query LatestAnalysisEntry($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                meta
                output
                sequencingGroups {
                    id
                }
                timestampCompleted
            }
        }
    }
""",
)

SPECIFIC_VDS_QUERY = gql(
    """
    query getVDSByAnalysisId($vds_id: Int!) {
        analyses(id: {eq: $vds_id}) {
            output
            sequencingGroups {
                id
            }
        }
    }
""",
)
SHARD_MANIFEST = 'shard-manifest.txt'


def query_for_specific_vds(vds_id: int) -> tuple[str, set[str]] | None:
    """
    query for a specific analysis of type entry_type for a dataset
    if found, return the set of SG IDs in the VDS (using the metadata)

    - stolen from the cpg_workflows.large_cohort.combiner Stage, but duplicated here so we can split pipelines without
      further code changes

    Args:
        vds_id (int): analysis id to query for

    Returns:
        either None if the analysis wasn't found, or a set of SG IDs in the VDS
    """

    # query for the exact, single analysis entry
    query_results: dict[str, dict] = query(SPECIFIC_VDS_QUERY, variables={'vds_id': vds_id})

    if not query_results['analyses']:
        return None
    vds_path: str = query_results['analyses'][0]['output']
    sg_ids = {sg['id'] for sg in query_results['analyses'][0]['sequencingGroups']}
    return vds_path, sg_ids


def query_for_latest_vds(dataset: str, entry_type: str = 'combiner') -> dict | None:
    """
    query for the latest analysis of type entry_type for a dataset
    Args:
        dataset (str): project to query for
        entry_type (str): type of analysis entry to query for
    Returns:
        str, the path to the latest analysis
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(LATEST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': entry_type})

    analyses_by_date = {}

    for analysis in result['project']['analyses']:
        if analysis['output'] and (
            analysis['meta']['sequencing_type'] == config.config_retrieve(['workflow', 'sequencing_type'])
        ):
            analyses_by_date[analysis['timestampCompleted']] = analysis

    if not analyses_by_date:
        loguru.logger.warning(f'No analysis of type {entry_type} found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort as strings
    return analyses_by_date[sorted(analyses_by_date)[-1]]


@cpg_flow.stage.stage(analysis_type='combiner')
class CreateVdsFromGvcfsWithHailCombinerStage(cpg_flow.stage.MultiCohortStage):
    """undecided if this will be reimplemented"""

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> dict[str, Path]:
        return self.prefix / f'{multicohort.name}.vds'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        #
        # outputs = self.expected_outputs(multicohort)
        #
        # # we only ever build on top of a single VDS, or start from scratch
        # vds_path: str | None = None
        #
        # # create these as empty iterables
        # sg_ids_in_vds: set[str] = set()
        # sgs_to_remove: list[str] = []
        #
        # # check for a VDS by ID
        # if vds_id := config_retrieve(['workflow', 'use_specific_vds'], None):
        #     vds_result_or_none = query_for_specific_vds(vds_id)
        #     if vds_result_or_none is None:
        #         raise ValueError(f'Specified VDS ID {vds_id} not found in Metamist')
        #
        #     # if not none, unpack the result
        #     vds_path, sg_ids_in_vds = vds_result_or_none
        #
        # # check for existing VDS by getting all and fetching latest
        # elif config_retrieve(['workflow', 'check_for_existing_vds']):
        #     loguru.logger.info('Checking for existing VDS')
        #     if existing_vds_analysis_entry := query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
        #         vds_path = existing_vds_analysis_entry['output']
        #         sg_ids_in_vds = {sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']}
        #
        # else:
        #     loguru.logger.info('Not continuing from any previous VDS, creating new Combiner from gVCFs only')
        #
        # # quick check - if we found a VDS, guarantee it exists
        # if vds_path and not exists(vds_path):
        #     raise ValueError(f'VDS {vds_path} does not exist, but has an Analysis Entry')
        #
        # # this is a quick and confident check on current VDS contents, but does require a direct connection to the VDS
        # # by default this is True, and can be deactivated in config
        # if vds_path and config_retrieve(['workflow', 'manually_check_vds_sg_ids']):
        #     sg_ids_in_vds = utils.manually_find_ids_from_vds(vds_path)
        #
        # # technicality; this can be empty - in a situation where we have a VDS already and the current MCohort has FEWER
        # # SG IDs in it, this stage will re-run because the specific hash will be different to the previous VDS
        # # See https://github.com/populationgenomics/production-pipelines/issues/1126
        # new_sg_gvcfs: list[str] = [
        #     str(sg.gvcf)
        #     for sg in multicohort.get_sequencing_groups()
        #     if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
        # ]
        #
        # # final check - if we have a VDS, and we have a current MultiCohort
        # # detect any samples which should be _removed_ from the current VDS prior to further combining taking place
        # if sg_ids_in_vds:
        #     sgs_in_mc: list[str] = multicohort.get_sequencing_group_ids()
        #     loguru.logger.info(f'Found {len(sg_ids_in_vds)} SG IDs in VDS {vds_path}')
        #     loguru.logger.info(f'Total {len(sgs_in_mc)} SGs in this MultiCohort')
        #
        #     sgs_to_remove = sorted(set(sg_ids_in_vds) - set(sgs_in_mc))
        #
        #     if sgs_to_remove:
        #         loguru.logger.info(f'Removing {len(sgs_to_remove)} SGs from VDS {vds_path}')
        #         loguru.logger.info(f'SGs to remove: {sgs_to_remove}')
        #
        # if not (new_sg_gvcfs or sgs_to_remove):
        #     loguru.logger.info('No GVCFs to add to, or remove from, existing VDS')
        #     loguru.logger.info(f'Checking if VDS exists: {outputs["vds"]}: {outputs["vds"].exists()}')  # type: ignore
        #     return self.make_outputs(multicohort, outputs)
        #
        # combiner_job = get_batch().new_python_job('CreateVdsFromGvcfsWithHailCombiner', {'stage': self.name})
        # combiner_job.image(config_retrieve(['workflow', 'driver_image']))
        # combiner_job.memory(config_retrieve(['combiner', 'driver_memory']))
        # combiner_job.storage(config_retrieve(['combiner', 'driver_storage']))
        # combiner_job.cpu(config_retrieve(['combiner', 'driver_cores']))
        #
        # # set this job to be non-spot (i.e. non-preemptible)
        # # previous issues with preemptible VMs led to multiple simultaneous QOB groups processing the same data
        # combiner_job.spot(config_retrieve(['combiner', 'preemptible_vms']))
        #
        # # Default to GRCh38 for reference if not specified
        # combiner_job.call(
        #     combiner.run,
        #     output_vds_path=str(outputs['vds']),
        #     save_path=str(self.tmp_prefix / 'combiner_plan.json'),
        #     sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
        #     tmp_prefix=str(self.tmp_prefix / 'temp_dir'),
        #     genome_build=genome_build(),
        #     gvcf_paths=new_sg_gvcfs,
        #     vds_path=vds_path,
        #     force_new_combiner=config_retrieve(['combiner', 'force_new_combiner'], False),
        #     sgs_to_remove=sgs_to_remove,
        # )

        return self.make_outputs(multicohort, data=None, jobs=None)


@cpg_flow.stage.stage(
    required_stages=[CreateVdsFromGvcfsWithHailCombinerStage],
    analysis_type='matrixtable',
    analysis_keys=['mt'],
)
class CreateDenseMtFromVdsWithHail(cpg_flow.stage.MultiCohortStage):
    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> dict:
        """
        the MT and both shard_manifest files are Paths, so this stage will rerun if any of those are missing
        the VCFs are written as a directory, rather than a single VCF, so we can't check its existence well

        Needs a range of INFO fields to be present in the VCF
        """
        prefix = self.prefix
        temp_prefix = self.tmp_prefix

        return {
            'mt': prefix / f'{multicohort.name}.mt',
            # this will be the write path for fragments of sites-only VCF (header-per-shard)
            'hps_vcf_dir': str(prefix / f'{multicohort.name}.vcf.bgz'),
            # this will be the file which contains the name of all fragments (header-per-shard)
            'hps_shard_manifest': prefix / f'{multicohort.name}.vcf.bgz' / SHARD_MANIFEST,
            # this will be the write path for fragments of sites-only VCF (separate header)
            'separate_header_vcf_dir': str(temp_prefix / f'{multicohort.name}_separate.vcf.bgz'),
            # this will be the file which contains the name of all fragments (separate header)
            'separate_header_manifest': temp_prefix / f'{multicohort.name}_separate.vcf.bgz' / SHARD_MANIFEST,
        }

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        outputs = self.expected_outputs(multicohort)

        job = jobs.CreateDenseMtFromVdsWithHail.generate_densify_jobs(
            input_vds=inputs.as_str(multicohort, CreateVdsFromGvcfsWithHailCombinerStage),
            output_mt=outputs['mt'],
            output_sites_only=outputs['hps_vcf_dir'],
            output_separate_header=outputs['separate_header_vcf_dir'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(target=multicohort, data=outputs, jobs=job)


@cpg_flow.stage.stage(required_stages=[CreateDenseMtFromVdsWithHail])
class ConcatenateVcfsWithGcloud(cpg_flow.stage.MultiCohortStage):
    """
    Takes a manifest of VCF fragments, and produces a single VCF file
    This is disconnected from the previous stage, but requires it to be run first
    So we check for the exact same output, and fail if we're not ready to start
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'gcloud_composed_sitesonly.vcf.bgz'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        """
        Submit jobs to take a manifest of VCF fragments, and produce a single VCF file
        The VCF being composed here has a single header in a separate file, the first entry in the manifest
        This means we can combine the VCF header and data fragments through concatenation
        and the result will be a spec-compliant VCF
        """
        dense_inputs = inputs.as_dict(
            target=multicohort,
            stage=CreateDenseMtFromVdsWithHail,
        )

        if not dense_inputs['separate_header_manifest'].exists():
            raise ValueError(
                f'Manifest file {dense_inputs["separate_header_manifest"]!s} does not exist, '
                f're-run the combiner workflow with workflows.last_stages=["CreateDenseMtFromVdsWithHailStage"]',
            )

        output = self.expected_outputs(multicohort)

        job = jobs.ConcatenateVcfFragmentsWithGcloud.create_and_run_compose_script(
            multicohort=multicohort,
            manifest_file=dense_inputs['separate_header_manifest'],
            manifest_dir=dense_inputs['separate_header_vcf_dir'],
            output=output,
            tmp_dir=self.tmp_prefix,
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@cpg_flow.stage.stage(required_stages=[ConcatenateVcfsWithGcloud])
class TrainVqsrIndelModel(cpg_flow.stage.MultiCohortStage):
    """
    Train VQSR Indel model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> dict[str, Path | str]:
        prefix = self.prefix
        return {
            'indel_recalibrations': prefix / 'indel.recal',
            'indel_tranches': prefix / 'indel.tranches',
            'indel_prefix': str(prefix / 'indel'),
        }

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        """
        Submit jobs to train VQSR on the combiner data.
        """

        outputs = self.expected_outputs(multicohort)

        composed_sitesonly_vcf = inputs.as_str(multicohort, ConcatenateVcfsWithGcloud)

        job = jobs.TrainVqsrIndels.train_vqsr_indel_model(
            sites_only_vcf=composed_sitesonly_vcf,
            output_prefix=outputs['indel_prefix'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@cpg_flow.stage.stage(required_stages=[ConcatenateVcfsWithGcloud])
class TrainVqsrSnpModel(cpg_flow.stage.MultiCohortStage):
    """
    Train VQSR SNP model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        return self.prefix / 'snp_model'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        """
        Submit jobs to train VQSR on the combiner data
        """

        composed_sitesonly_vcf = inputs.as_str(multicohort, ConcatenateVcfsWithGcloud)
        outputs = self.expected_outputs(multicohort)
        job = jobs.TrainVqsrSnps.train_vqsr_snp_model(
            sites_only_vcf=composed_sitesonly_vcf,
            snp_model=str(outputs),
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@cpg_flow.stage.stage(required_stages=[CreateDenseMtFromVdsWithHail, TrainVqsrSnpModel])
class TrainVqsrSnpTranches(cpg_flow.stage.MultiCohortStage):
    """
    Scattered training of VQSR tranches for SNPs
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'tranches_trained'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        output = self.expected_outputs(multicohort)

        manifest_file = inputs.as_str(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')

        snp_model_path = inputs.as_str(target=multicohort, stage=TrainVqsrSnpModel)

        job_list = jobs.TrainVqsrSnpTranches.train_vqsr_snp_tranches(
            manifest_file=manifest_file,
            snp_model_path=snp_model_path,
            output_path=output,
            temp_path=self.tmp_prefix / 'vqsr_snp_tranches',
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=output, jobs=job_list)


@cpg_flow.stage.stage(required_stages=[CreateDenseMtFromVdsWithHail, TrainVqsrSnpTranches])
class GatherTrainedVqsrSnpTranches(cpg_flow.stage.MultiCohortStage):
    """Scattered training of VQSR tranches for SNPs."""

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        return self.prefix / 'snp_tranches'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        manifest_file = inputs.as_path(
            target=multicohort,
            stage=CreateDenseMtFromVdsWithHail,
            key='hps_shard_manifest',
        )

        output = self.expected_outputs(multicohort)

        job = jobs.GatherTrainedVqsrSnpTranches.gather_tranches(
            manifest_file=manifest_file,
            temp_path=self.tmp_prefix / 'vqsr_snp_tranches',
            output_path=output,
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=output, jobs=job)


@cpg_flow.stage.stage(
    required_stages=[
        CreateDenseMtFromVdsWithHail,
        GatherTrainedVqsrSnpTranches,
        TrainVqsrSnpTranches,
    ],
)
class RunSnpVqsrOnFragments(cpg_flow.stage.MultiCohortStage):
    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'vqsr.vcf.gz'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        output = self.expected_outputs(multicohort)

        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')

        tranche_file = inputs.as_str(target=multicohort, stage=GatherTrainedVqsrSnpTranches)

        job_list = jobs.RunSnpVqsrOnFragments.apply_snp_vqsr_to_fragments(
            manifest_file=manifest_file,
            tranche_file=tranche_file,
            temp_path=self.tmp_prefix / 'vqsr_snp_tranches',
            output_path=str(output),
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=output, jobs=job_list)


@cpg_flow.stage.stage(
    analysis_type='qc',
    required_stages=[
        RunSnpVqsrOnFragments,
        TrainVqsrIndelModel,
    ],
)
class RunIndelVqsr(cpg_flow.stage.MultiCohortStage):
    """
    Run Indel VQSR on the reconstituted, SNP-annotated, VCF
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        return self.prefix / 'vqsr_snps_and_indels.vcf.gz'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        output = self.expected_outputs(multicohort)
        annotated_vcf = inputs.as_str(target=multicohort, stage=RunSnpVqsrOnFragments)
        indel_recalibrations = inputs.as_str(
            target=multicohort,
            stage=TrainVqsrIndelModel,
            key='indel_recalibrations',
        )
        indel_tranches = inputs.as_str(
            target=multicohort,
            stage=TrainVqsrIndelModel,
            key='indel_tranches',
        )

        job = jobs.RunIndelVqsr.apply_recalibration_indels(
            snp_annotated_vcf=annotated_vcf,
            indel_recalibration=indel_recalibrations,
            indel_tranches=indel_tranches,
            output_path=str(output),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(
            multicohort,
            data=output,
            jobs=job,
        )


# TODO use a proper analysis type
@cpg_flow.stage.stage(
    analysis_type='custom',
    required_stages=[CreateDenseMtFromVdsWithHail],
)
class AnnotateVcfsWithVep(cpg_flow.stage.MultiCohortStage):
    """
    Annotate VCF with VEP.
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        """
        Should this be in tmp? We'll never use it again maybe?
        """
        return self.prefix / 'vep.ht'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        outputs = self.expected_outputs(multicohort)
        manifest_file = inputs.as_path(
            target=multicohort,
            stage=CreateDenseMtFromVdsWithHail,
            key='hps_shard_manifest',
        )

        vep_jobs = jobs.AnnotateVcfsWithVep.add_vep_jobs(
            manifest_file=manifest_file,
            final_out_path=outputs,
            tmp_prefix=self.tmp_prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=vep_jobs)


@cpg_flow.stage.stage(
    required_stages=[
        CreateDenseMtFromVdsWithHail,
        AnnotateVcfsWithVep,
        RunIndelVqsr,
    ],
)
class AnnotateCohort(cpg_flow.stage.MultiCohortStage):
    """
    Annotate SNP/Indel MT with VEP and VQSR.
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return self.tmp_prefix / 'annotate_cohort.mt'

    def queue_jobs(
        self, multicohort: cpg_flow.targets.MultiCohort, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        """

        Args:
            multicohort ():
            inputs ():
        """

        outputs = self.expected_outputs(multicohort)
        vep_ht_path = inputs.as_str(target=multicohort, stage=AnnotateVcfsWithVep)
        vqsr_vcf = inputs.as_str(target=multicohort, stage=RunIndelVqsr)
        variant_mt = inputs.as_str(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='mt')

        job = jobs.AnnotateCohort.create_annotate_cohort_job(
            output_mt=outputs,
            vep_ht=vep_ht_path,
            vqsr_vcf=vqsr_vcf,
            variant_mt=variant_mt,
            checkpoint_prefix=self.tmp_prefix / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@cpg_flow.stage.stage(required_stages=[AnnotateCohort])
class SubsetMtToDatasetWithHail(cpg_flow.stage.DatasetStage):
    """
    Subset the MT to a single dataset - or a subset of families within a dataset
    Skips this stage if the MultiCohort has only one dataset
    """

    def expected_outputs(self, dataset: cpg_flow.targets.Dataset) -> dict[str, Path] | None:
        """
        Expected to generate a MatrixTable
        This is kinda transient, so shove it in tmp.

        If subsetting to certain families, the output will be named accordingly
        """
        output_prefix = dataset.tmp_prefix() / 'mt' / self.name
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            return {
                'mt': output_prefix
                / f'{cpg_flow.workflow.get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.mt',
                'id_file': dataset.tmp_prefix()
                / f'{cpg_flow.workflow.get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}-SG-ids.txt',
            }
        elif len(cpg_flow.workflow.get_multicohort().get_datasets()) == 1:
            loguru.logger.info(f'Skipping SubsetMatrixTableToDataset for single Dataset {dataset}')
            return None
        else:
            return {
                'mt': output_prefix / f'{cpg_flow.workflow.get_workflow().output_version}-{dataset.name}.mt',
                'id_file': dataset.tmp_prefix()
                / f'{cpg_flow.workflow.get_workflow().output_version}-{dataset.name}-SG-ids.txt',
            }

    def queue_jobs(
        self, dataset: cpg_flow.targets.Dataset, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        # only create dataset MTs for datasets specified in the config
        # and only run this stage if the callset has multiple datasets
        if (outputs is None) or (
            dataset.name not in config.config_retrieve(['workflow', 'write_mt_for_datasets'], default=[])
        ):
            loguru.logger.info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        variant_mt = inputs.as_path(target=cpg_flow.workflow.get_multicohort(), stage=AnnotateCohort)

        job = jobs.SubsetMtToDatasetWithHail.create_subset_mt_job(
            dataset=dataset,
            input_mt=variant_mt,
            id_file=outputs['id_file'],
            output_mt=outputs['mt'],
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=job)


@cpg_flow.stage.stage(
    required_stages=[
        AnnotateCohort,
        SubsetMtToDatasetWithHail,
    ],
    analysis_type='matrixtable',
)
class AnnotateDataset(cpg_flow.stage.DatasetStage):
    def expected_outputs(self, dataset: cpg_flow.targets.Dataset) -> Path:
        """
        Expected to generate a matrix table
        """
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            return (
                dataset.prefix()
                / 'mt'
                / self.name
                / f'{cpg_flow.workflow.get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.mt'
            )
        else:
            return (
                dataset.prefix()
                / 'mt'
                / self.name
                / f'{cpg_flow.workflow.get_workflow().output_version}-{dataset.name}.mt'
            )

    def queue_jobs(
        self, dataset: cpg_flow.targets.Dataset, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput:
        # only create final MTs for datasets specified in the config
        if dataset.name not in config.config_retrieve(['workflow', 'write_mt_for_datasets'], default=[]):
            loguru.logger.info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        output = self.expected_outputs(dataset)

        # choose whether to run directly from AnnotateCohort, or use the Dataset/Family Subset MT from previous Stage
        family_sgs = utils.get_family_sequencing_groups(dataset)
        # choose the input MT based on the number of datasets in the MultiCohort and the presence of family SGs
        if len(cpg_flow.workflow.get_multicohort().get_datasets()) == 1 and family_sgs is None:
            input_mt = inputs.as_path(target=cpg_flow.workflow.get_multicohort(), stage=AnnotateCohort)
        else:
            input_mt = inputs.as_path(target=dataset, stage=SubsetMtToDatasetWithHail, key='mt')

        job = jobs.AnnotateDataset.create_annotate_cohort_job(
            dataset=dataset,
            input_mt=input_mt,
            output_mt=output,
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=output, jobs=job)


@cpg_flow.stage.stage(
    required_stages=[AnnotateDataset],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'VARIANTS'},
)
class ExportMtAsEsIndex(cpg_flow.stage.DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: cpg_flow.targets.Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            index_name = f'{dataset.name}-{sequencing_type}-{family_sgs["name_suffix"]}-{cpg_flow.workflow.get_workflow().run_timestamp}'.lower()
        else:
            index_name = f'{dataset.name}-{sequencing_type}-{cpg_flow.workflow.get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(
        self, dataset: cpg_flow.targets.Dataset, inputs: cpg_flow.stage.StageInput
    ) -> cpg_flow.stage.StageOutput | None:
        """
        Transforms the MT into a Seqr index, no DataProc
        """
        # only create the elasticsearch index for the datasets specified in the config
        eligible_datasets = config.config_retrieve(['workflow', 'create_es_index_for_datasets'], default=[])
        if dataset.name not in eligible_datasets:
            loguru.logger.info(f'Skipping ES index creation for {dataset}')
            return None

        # try to generate a password here - we'll find out inside the script anyway, but
        # by that point we'd already have localised the MT, wasting time and money
        try:
            _es_password_string = cloud.read_secret(
                project_id=config.config_retrieve(['elasticsearch', 'password_project_id']),
                secret_name=config.config_retrieve(['elasticsearch', 'password_secret_id']),
                fail_gracefully=False,
            )
        except PermissionDenied:
            loguru.logger.warning(f'No permission to access ES password, skipping for {dataset}')
            return self.make_outputs(dataset)
        except KeyError:
            loguru.logger.warning(f'ES section not in config, skipping for {dataset}')
            return self.make_outputs(dataset)

        # get the absolute path to the MT
        mt_path = inputs.as_str(target=dataset, stage=AnnotateDataset)

        outputs = self.expected_outputs(dataset)

        job = jobs.ExportMtAsEsIndex.create_annotate_cohort_job(
            index_name=outputs['index_name'],
            done_flag=outputs['done_flag'],
            mt_path=mt_path,
            sample_size=len(dataset.get_sequencing_group_ids()),
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
