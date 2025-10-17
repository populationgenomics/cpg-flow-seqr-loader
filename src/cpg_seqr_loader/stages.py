"""
All Stages relating to the seqr_loader pipeline, reimplemented from scratch to
use the gVCF combiner instead of joint-calling.
"""

import loguru
from cpg_flow import stage, targets, workflow
from cpg_flow.stage import StageOutput
from cpg_utils import Path, cloud, config
from google.api_core.exceptions import PermissionDenied

from cpg_seqr_loader import utils
from cpg_seqr_loader.jobs.AnnotateCohort import create_annotate_cohort_job
from cpg_seqr_loader.jobs.AnnotateDataset import create_annotate_dataset_job
from cpg_seqr_loader.jobs.AnnotatedDatasetMtToVcf import cohort_to_vcf_job
from cpg_seqr_loader.jobs.AnnotateVcfsWithVep import add_vep_jobs
from cpg_seqr_loader.jobs.CombineGvcfsIntoVds import create_combiner_jobs
from cpg_seqr_loader.jobs.ConcatenateVcfFragmentsWithGcloud import create_and_run_compose_script
from cpg_seqr_loader.jobs.CreateDenseMtFromVdsWithHail import generate_densify_jobs
from cpg_seqr_loader.jobs.DeleteCombinerTemp import delete_temp_data_recursive
from cpg_seqr_loader.jobs.ExportMtAsEsIndex import create_es_export_job
from cpg_seqr_loader.jobs.GatherTrainedVqsrSnpTranches import gather_tranches
from cpg_seqr_loader.jobs.RunIndelVqsr import apply_recalibration_indels
from cpg_seqr_loader.jobs.RunSnpVqsrOnFragments import apply_snp_vqsr_to_fragments
from cpg_seqr_loader.jobs.SubsetMtToDatasetWithHail import create_subset_mt_job
from cpg_seqr_loader.jobs.TrainVqsrIndels import train_vqsr_indel_model
from cpg_seqr_loader.jobs.TrainVqsrSnpModel import train_vqsr_snp_model
from cpg_seqr_loader.jobs.TrainVqsrSnpTranches import train_vqsr_snp_tranches

SHARD_MANIFEST = 'shard-manifest.txt'


@stage.stage(analysis_type='combiner', analysis_keys=['vds'])
class CombineGvcfsIntoVds(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path | str]:
        return {
            'vds': self.prefix / f'{multicohort.name}.vds',
            'tmp': str(self.tmp_prefix / 'temp_dir'),
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)

        job = create_combiner_jobs(
            multicohort=multicohort,
            output_vds=outputs['vds'],
            combiner_plan=self.tmp_prefix / 'combiner_plan.json',
            temp_dir_string=outputs['tmp'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=CombineGvcfsIntoVds)
class DeleteCombinerTemp(stage.MultiCohortStage):
    """
    new stage for purging the temp data created during VDS combining

    This could be an additional job in the previous stage, but it is a separate job
    VDSs can take a LONG time to delete given the high level of fragmentation, and the previous step can generate
    many VDSs. As such, if this was in the previous stage, the start of the next step (Creating a MT from the VDS)
    would be on hold until all VDSs were deleted, which could take a long time.

    Splitting this out into a separate dependent stage allows the next step to start immediately, whilst deletions are
    going on in parallel. It also allows this Stage to be skipped through configuration, if the user does not want to
    explicitly delete the temporary data created by the combiner.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'CombinerIntermediatesDeleted.txt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> StageOutput:
        combiner_outputs = inputs.as_dict(multicohort, CombineGvcfsIntoVds)

        output = self.expected_outputs(multicohort)

        job = delete_temp_data_recursive(
            temp_path=combiner_outputs['tmp'],
            logfile=output,
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(required_stages=CombineGvcfsIntoVds)
class CreateDenseMtFromVdsWithHail(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        """
        the MT and both shard_manifest files are Paths, so this stage will rerun if any of those are missing
        the VCFs are written as a directory, rather than a single VCF, so we can't check its existence well

        Needs a range of INFO fields to be present in the VCF
        """
        temp_prefix = self.tmp_prefix

        return {
            'mt': temp_prefix / f'{multicohort.name}.mt',
            # this will be the write path for fragments of sites-only VCF (header-per-shard)
            'hps_vcf_dir': str(temp_prefix / f'{multicohort.name}.vcf.bgz'),
            # this will be the file which contains the name of all fragments (header-per-shard)
            'hps_shard_manifest': temp_prefix / f'{multicohort.name}.vcf.bgz' / SHARD_MANIFEST,
            # this will be the write path for fragments of sites-only VCF (separate header)
            'separate_header_vcf_dir': str(temp_prefix / f'{multicohort.name}_separate.vcf.bgz'),
            # this will be the file which contains the name of all fragments (separate header)
            'separate_header_manifest': temp_prefix / f'{multicohort.name}_separate.vcf.bgz' / SHARD_MANIFEST,
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)

        job = generate_densify_jobs(
            input_vds=inputs.as_str(multicohort, CombineGvcfsIntoVds, 'vds'),
            output_mt=outputs['mt'],
            output_sites_only=outputs['hps_vcf_dir'],
            output_separate_header=outputs['separate_header_vcf_dir'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(target=multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=[CreateDenseMtFromVdsWithHail])
class ConcatenateVcfsWithGcloud(stage.MultiCohortStage):
    """
    Takes a manifest of VCF fragments, and produces a single VCF file
    This is disconnected from the previous stage, but requires it to be run first
    So we check for the exact same output, and fail if we're not ready to start
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'gcloud_composed_sitesonly.vcf.bgz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
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

        job = create_and_run_compose_script(
            multicohort=multicohort,
            manifest_file=dense_inputs['separate_header_manifest'],
            manifest_dir=dense_inputs['separate_header_vcf_dir'],
            output=output,
            tmp_dir=self.tmp_prefix,
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(required_stages=[ConcatenateVcfsWithGcloud])
class TrainVqsrIndelModel(stage.MultiCohortStage):
    """
    Train VQSR Indel model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path | str]:
        prefix = self.tmp_prefix
        return {
            'indel_recalibrations': prefix / 'indel.recal',
            'indel_tranches': prefix / 'indel.tranches',
            'indel_prefix': str(prefix / 'indel'),
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Submit jobs to train VQSR on the combiner data.
        """

        outputs = self.expected_outputs(multicohort)

        composed_sitesonly_vcf = inputs.as_str(multicohort, ConcatenateVcfsWithGcloud)

        job = train_vqsr_indel_model(
            sites_only_vcf=composed_sitesonly_vcf,
            output_prefix=outputs['indel_prefix'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=[ConcatenateVcfsWithGcloud])
class TrainVqsrSnpModel(stage.MultiCohortStage):
    """
    Train VQSR SNP model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'snp_model'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Submit jobs to train VQSR on the combiner data
        """

        composed_sitesonly_vcf = inputs.as_str(multicohort, ConcatenateVcfsWithGcloud)
        outputs = self.expected_outputs(multicohort)
        job = train_vqsr_snp_model(
            sites_only_vcf=composed_sitesonly_vcf,
            snp_model=str(outputs),
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=[CreateDenseMtFromVdsWithHail, TrainVqsrSnpModel])
class TrainVqsrSnpTranches(stage.MultiCohortStage):
    """
    Scattered training of VQSR tranches for SNPs
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path | str]:
        return {
            'tranches': self.tmp_prefix / 'tranches_trained',
            'temp': self.tmp_prefix / 'vqsr_snp_tranches',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)

        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')

        snp_model_path = inputs.as_str(target=multicohort, stage=TrainVqsrSnpModel)

        job_list = train_vqsr_snp_tranches(
            manifest_file=manifest_file,
            snp_model_path=snp_model_path,
            output_path=outputs['tranches'],
            temp_path=outputs['temp'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job_list)


@stage.stage(required_stages=[CreateDenseMtFromVdsWithHail, TrainVqsrSnpTranches])
class GatherTrainedVqsrSnpTranches(stage.MultiCohortStage):
    """Scattered training of VQSR tranches for SNPs."""

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'snp_tranches'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        manifest_file = inputs.as_path(
            target=multicohort,
            stage=CreateDenseMtFromVdsWithHail,
            key='hps_shard_manifest',
        )

        output = self.expected_outputs(multicohort)

        # temp dir for outputs from the previous stage
        temp_dir = inputs.as_path(multicohort, TrainVqsrSnpTranches, key='temp')

        job = gather_tranches(
            manifest_file=manifest_file,
            temp_path=temp_dir,
            output_path=output,
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(
    required_stages=[
        CreateDenseMtFromVdsWithHail,
        GatherTrainedVqsrSnpTranches,
        TrainVqsrSnpTranches,
    ],
)
class RunSnpVqsrOnFragments(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'vqsr.vcf.gz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)

        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')

        tranche_file = inputs.as_str(target=multicohort, stage=GatherTrainedVqsrSnpTranches)

        # temp dir for outputs from the previous stage
        temp_dir = inputs.as_path(multicohort, TrainVqsrSnpTranches, key='temp')

        job_list = apply_snp_vqsr_to_fragments(
            manifest_file=manifest_file,
            tranche_file=tranche_file,
            temp_path=temp_dir,
            output_path=str(output),
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=output, jobs=job_list)


@stage.stage(
    required_stages=[
        RunSnpVqsrOnFragments,
        TrainVqsrIndelModel,
    ],
)
class RunIndelVqsr(stage.MultiCohortStage):
    """
    Run Indel VQSR on the reconstituted, SNP-annotated, VCF
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.tmp_prefix / 'vqsr_snps_and_indels.vcf.gz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
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

        job = apply_recalibration_indels(
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


@stage.stage(required_stages=[CreateDenseMtFromVdsWithHail])
class AnnotateVcfsWithVep(stage.MultiCohortStage):
    """
    Annotate VCF with VEP.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        """
        Should this be in tmp? We'll never use it again maybe?
        """
        return self.tmp_prefix / 'vep.ht'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)
        manifest_file = inputs.as_path(
            target=multicohort,
            stage=CreateDenseMtFromVdsWithHail,
            key='hps_shard_manifest',
        )

        vep_jobs = add_vep_jobs(
            manifest_file=manifest_file,
            final_out_path=outputs,
            tmp_prefix=self.tmp_prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=vep_jobs)


@stage.stage(
    required_stages=[
        CreateDenseMtFromVdsWithHail,
        AnnotateVcfsWithVep,
        RunIndelVqsr,
    ],
)
class AnnotateCohort(stage.MultiCohortStage):
    """
    Annotate SNP/Indel MT with VEP and VQSR.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return self.prefix / 'annotate_cohort.mt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:

        outputs = self.expected_outputs(multicohort)
        vep_ht_path = inputs.as_str(target=multicohort, stage=AnnotateVcfsWithVep)
        vqsr_vcf = inputs.as_str(target=multicohort, stage=RunIndelVqsr)
        variant_mt = inputs.as_str(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='mt')

        job = create_annotate_cohort_job(
            output_mt=outputs,
            vep_ht=vep_ht_path,
            vqsr_vcf=vqsr_vcf,
            variant_mt=variant_mt,
            checkpoint_prefix=self.tmp_prefix / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=AnnotateCohort)
class SubsetMtToDatasetWithHail(stage.DatasetStage):
    """
    Subset the MT to a single dataset - or a subset of families within a dataset
    Skips this stage if the MultiCohort has only one dataset
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path] | None:
        """
        Expected to generate a MatrixTable
        This is kinda transient, so shove it in tmp.

        If subsetting to certain families, the output will be named accordingly
        """
        output_prefix = dataset.tmp_prefix() / 'mt' / self.name
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            return {
                'mt': output_prefix
                / f'{workflow.get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.mt',
                'id_file': dataset.tmp_prefix()
                / f'{workflow.get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}-SG-ids.txt',
            }
        if len(workflow.get_multicohort().get_datasets()) == 1:
            loguru.logger.info(f'Skipping SubsetMatrixTableToDataset for single Dataset {dataset}')
            return None

        return {
            'mt': output_prefix / f'{workflow.get_workflow().output_version}-{dataset.name}.mt',
            'id_file': dataset.tmp_prefix() / f'{workflow.get_workflow().output_version}-{dataset.name}-SG-ids.txt',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        # only create dataset MTs for datasets specified in the config
        # and only run this stage if the callset has multiple datasets
        if (outputs is None) or (
            dataset.name not in config.config_retrieve(['workflow', 'write_mt_for_datasets'], default=[])
        ):
            loguru.logger.info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        variant_mt = inputs.as_path(target=workflow.get_multicohort(), stage=AnnotateCohort)

        job = create_subset_mt_job(
            dataset=dataset,
            input_mt=variant_mt,
            id_file=outputs['id_file'],
            output_mt=outputs['mt'],
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=job)


@stage.stage(
    required_stages=[
        AnnotateCohort,
        SubsetMtToDatasetWithHail,
    ],
    # analyses recorded to pass through to Metamist/other workflows
    analysis_type='matrixtable',
)
class AnnotateDataset(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        """
        Expected to generate a matrix table
        """
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            mt_name = f'{dataset.name}-{family_sgs["name_suffix"]}.mt'
        else:
            mt_name = f'{dataset.name}.mt'

        return (
            dataset.prefix()
            / workflow.get_workflow().name
            / workflow.get_workflow().output_version
            / self.name
            / mt_name
        )

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # only create final MTs for datasets specified in the config
        # catches specific instructions to just run this Stage, and for the two stages which can follow this
        if not utils.run_annotate_dataset(dataset.name):
            loguru.logger.info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        output = self.expected_outputs(dataset)

        # choose whether to run directly from AnnotateCohort, or use the Dataset/Family Subset MT from previous Stage
        family_sgs = utils.get_family_sequencing_groups(dataset)
        # choose the input MT based on the number of datasets in the MultiCohort and the presence of family SGs
        if len(workflow.get_multicohort().get_datasets()) == 1 and family_sgs is None:
            input_mt = inputs.as_path(target=workflow.get_multicohort(), stage=AnnotateCohort)
        else:
            input_mt = inputs.as_path(target=dataset, stage=SubsetMtToDatasetWithHail, key='mt')

        job = create_annotate_dataset_job(
            dataset=dataset,
            input_mt=input_mt,
            output_mt=output,
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=output, jobs=job)


@stage.stage(required_stages=[AnnotateDataset], analysis_type='custom', analysis_keys=['vcf'])
class AnnotatedDatasetMtToVcf(stage.DatasetStage):
    """
    Take the per-dataset annotated MT and write out as a VCF
    Optional stage set by dataset name in the config file
    """

    def expected_outputs(self, dataset: targets.Dataset):
        """
        Expected to generate a VCF from the single-dataset MT
        """
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            return (
                dataset.prefix()
                / 'vcf'
                / f'{workflow.get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.vcf.bgz'
            )
        return dataset.prefix() / 'vcf' / f'{workflow.get_workflow().output_version}-{dataset.name}.vcf.bgz'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """Run a MT -> VCF extraction on selected cohorts - only run this on manually defined list of Datasets."""

        # only run this selectively, most datasets it's not required
        eligible_datasets = config.config_retrieve(['workflow', 'write_vcf'])
        if dataset.name not in eligible_datasets:
            return None

        output = self.expected_outputs(dataset)

        mt_path = inputs.as_str(target=dataset, stage=AnnotateDataset)

        job = cohort_to_vcf_job(
            input_mt=mt_path,
            output_vcf=output,
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=output, jobs=job)


@stage.stage(
    required_stages=[AnnotateDataset],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'VARIANTS'},  # noqa: ARG005
)
class ExportMtAsEsIndex(stage.DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
        prelude = f'{dataset.name}-{sequencing_type}'
        if family_sgs := utils.get_family_sequencing_groups(dataset):
            index_name = f'{prelude}-{family_sgs["name_suffix"]}-{workflow.get_workflow().output_version}'.lower()
        else:
            index_name = f'{prelude}-{workflow.get_workflow().output_version}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
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

        job = create_es_export_job(
            index_name=outputs['index_name'],
            done_flag=outputs['done_flag'],
            mt_path=mt_path,
            sample_size=len(dataset.get_sequencing_group_ids()),
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
