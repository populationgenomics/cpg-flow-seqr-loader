"""
Stages for the synthetic-proband workflow.

The synthetic-proband workflow runs in parallel with the main seqr-loader combiner
work: it invents proband gVCFs so a cohort of parental duos can be shaped into a
trio for the combiner, and its outputs feed a separate combiner run whose final
annotations are joined against the main callset (rather than re-annotated from
scratch, which would let the synthetic data pollute AC/AN/AF and friends).

This module is kept separate from `stages.py` so the two workflows don't share a
file - see PR review discussion. Only the synthetic-proband workflow entry point
(`synthetic_proband_gvcf_workflow.py`) imports from here.
"""

import hashlib

from cpg_flow import stage, targets, workflow
from cpg_utils import Path, config

from cpg_seqr_loader import utils
from cpg_seqr_loader.jobs.AnnotateFromGlobalCallset import create_annotate_from_global_callset_job
from cpg_seqr_loader.jobs.CombineGvcfsIntoVdsFromManifest import create_combiner_from_manifest_jobs
from cpg_seqr_loader.jobs.GenerateSyntheticProbandCombinerInputs import write_combiner_inputs
from cpg_seqr_loader.jobs.GenerateSyntheticProbandGvcfs import create_synthetic_gvcf_jobs
from cpg_seqr_loader.stages import (
    AnnotateDataset,
    CreateDenseMtFromVdsWithHail,
    ExportMtAsEsIndex,
    SubsetMtToDatasetWithHail,
)


@stage.stage
class GenerateSyntheticProbandGvcfs(stage.MultiCohortStage):
    """
    Generate synthetic proband gVCFs for each duo family in the multicohort, so the combiner can
    build a trio-shaped VDS for cohorts that only contain unaffected parental duos.

    Metamist analysis registration for each output is queued separately (see Task 4) rather than
    via the @stage.stage(analysis_type=...) decorator, because the stage produces N outputs (one
    per family) and the decorator only supports one analysis registration per stage run.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        """Two paths per qualifying family: the gVCF, and a sentinel written by the registration
        script on success. cpg_flow only skips this stage as [REUSE] when *all* keys point at
        existing files, so both artifacts must be present. Tracking the sentinel is what stops
        the framework from silently skipping metamist registration when the gVCF is already on
        disk from a previous run (see the debugging trail in git history).
        """
        families = utils.get_families_for_synthetic_probands(multicohort)
        outputs: dict[str, Path] = {}
        for family in families:
            # Filenames use the external family ID so collaborators recognise the artefacts.
            # Dict keys use the same, so downstream lookups stay symmetrical.
            outputs[f'{family.external_family_id}_gvcf'] = (
                self.prefix / f'{family.external_family_id}_synthetic_proband.g.vcf.gz'
            )
            outputs[f'{family.external_family_id}_registered'] = (
                self.prefix / f'{family.external_family_id}_registered.txt'
            )
        return outputs

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        families = utils.get_families_for_synthetic_probands(multicohort)
        outputs = self.expected_outputs(multicohort)

        gvcf_paths = {f.family_id: outputs[f'{f.external_family_id}_gvcf'] for f in families}
        marker_paths = {f.family_id: outputs[f'{f.external_family_id}_registered'] for f in families}

        jobs = create_synthetic_gvcf_jobs(
            families=families,
            gvcf_paths=gvcf_paths,
            marker_paths=marker_paths,
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[GenerateSyntheticProbandGvcfs])
class GenerateSyntheticProbandCombinerInputs(stage.MultiCohortStage):
    """
    Build the combiner inputs for the separate synthetic-trio combiner run: a pedigree with the
    synthetic probands inserted, and a gVCF manifest listing every real parental gVCF in the
    multicohort plus every synthetic gVCF from Stage 1.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        """PED + manifest, both scoped by a hash of the SG set and qualifying family set.

        The hash lives in the directory prefix, so any change to who's in the multicohort or
        which families qualify writes to a new directory rather than overwriting the previous
        run's outputs. cpg_flow's REUSE check keys off the hash-scoped paths themselves - if the
        family set changes, the new paths don't exist yet, and Stage 2 re-runs; if the family set
        is unchanged, both paths exist and the stage is skipped.

        Downstream stages (combiner, seqr sync) should pick the paths up via
        `inputs.as_path(stage=GenerateSyntheticProbandCombinerInputs, key=...)` rather than
        hard-coding the location, so they resolve to the current-config paths automatically.
        """
        families = utils.get_families_for_synthetic_probands(multicohort)
        sg_ids = sorted(sg.id for sg in multicohort.get_sequencing_groups())
        family_ids = sorted(family.family_id for family in families)
        family_set_hash = hashlib.sha256(('\n'.join([*sg_ids, '--', *family_ids])).encode()).hexdigest()[:12]
        scoped = self.prefix / f'families-{family_set_hash}'
        return {
            'gvcfs_list': scoped / 'all_gvcf_paths_including_synthetic_gvcfs.txt',
            'pedigree': scoped / 'synthetic_pedigree.ped',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        families = utils.get_families_for_synthetic_probands(multicohort)
        outputs = self.expected_outputs(multicohort)

        # Stage 1's outputs contain both `_gvcf` and `_registered` keys per family; we only want
        # the gVCF paths here (indexed by internal family_id, as the manifest builder expects).
        stage_1_outputs = inputs.as_dict(
            target=multicohort,
            stage=GenerateSyntheticProbandGvcfs,
        )
        synthetic_gvcf_paths = {
            family.family_id: stage_1_outputs[f'{family.external_family_id}_gvcf'] for family in families
        }

        # Both output files are written directly from the driver (see write_combiner_inputs) so
        # there are no Batch jobs to queue. cpg_flow's REUSE check keys off the output paths
        # existing on disk, which is fine with jobs=[].
        write_combiner_inputs(
            families=families,
            multicohort=multicohort,
            synthetic_gvcf_paths=synthetic_gvcf_paths,
            output_ped=outputs['pedigree'],
            output_manifest=outputs['gvcfs_list'],
        )
        return self.make_outputs(multicohort, data=outputs, jobs=[])


@stage.stage(
    analysis_type='combiner',
    analysis_keys=['vds'],
    required_stages=[GenerateSyntheticProbandCombinerInputs],
)
class CombineGvcfsIntoVdsFromManifest(stage.MultiCohortStage):
    """Combine the synthetic-trio gVCF manifest into an isolated VDS.

    Rebuilds from scratch every run using the manifest produced by
    GenerateSyntheticProbandCombinerInputs. The VDS is scoped to this cohort so synthetic
    proband data never contaminates AC/AN/AF in the global callset.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path | str]:
        return {
            'vds': self.prefix / f'{multicohort.name}.vds',
            'tmp': str(self.tmp_prefix / 'temp_dir'),
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)
        manifest_path = inputs.as_path(
            target=multicohort,
            stage=GenerateSyntheticProbandCombinerInputs,
            key='gvcfs_list',
        )
        job = create_combiner_from_manifest_jobs(
            manifest_path=manifest_path,
            output_vds=outputs['vds'],
            combiner_plan=self.tmp_prefix / 'combiner_plan.json',
            temp_dir_string=outputs['tmp'],
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=CombineGvcfsIntoVdsFromManifest)
class CreateDenseMtFromVdsWithHailNoFragments(CreateDenseMtFromVdsWithHail):
    """Densify variant for workflows that skip VQSR / VEP.

    Overrides the base class on two axes:
      - Reads the input VDS from the manifest-driven combiner instead of the standard one.
      - Omits the four sites-only VCF-fragment outputs, which the base class emits to feed
        VQSR / VEP downstream. This workflow joins pre-computed annotations from the global
        callset instead, so the fragments are dead weight.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        return {'mt': self.tmp_prefix / f'{multicohort.name}.mt'}

    def _get_input_vds(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> str:
        return inputs.as_str(multicohort, CombineGvcfsIntoVdsFromManifest, 'vds')


@stage.stage(
    required_stages=[CreateDenseMtFromVdsWithHailNoFragments],
    analysis_type='matrixtable',
)
class AnnotateFromGlobalCallset(stage.MultiCohortStage):
    """Attach row annotations from the standard seqr-loader's annotate_cohort.mt to this cohort.

    Replaces the standard AnnotateCohort stage for workflows that must skip the annotation
    stack (VEP / VQSR / gnomAD / clinvar) because the cohort contains synthetic samples that
    would pollute the resulting population stats. AC/AN/AF, VQSR fields, VEP consequences,
    reference joins, and clinvar are all copied from the global MT by (locus, alleles) join.

    The global MT path is resolved from metamist via query_for_latest_annotate_cohort_mt,
    unless overridden by config key `annotate_from_global_callset.source_mt`. Both the
    metamist query and the resolved path must succeed at DAG-planning time - the stage
    fails loudly rather than let the pipeline start with a missing or stale source.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'annotate_cohort.mt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)
        input_mt = inputs.as_str(target=multicohort, stage=CreateDenseMtFromVdsWithHailNoFragments, key='mt')

        global_mt_override = config.config_retrieve(['annotate_from_global_callset', 'source_mt'], None)
        if global_mt_override:
            global_mt = str(global_mt_override)
        else:
            source_dataset = config.config_retrieve(['annotate_from_global_callset', 'source_dataset'], 'seqr')
            global_mt = utils.query_for_latest_annotate_cohort_mt(source_dataset)

        job = create_annotate_from_global_callset_job(
            input_mt=input_mt,
            global_mt=global_mt,
            output_mt=outputs,
            job_attrs=self.get_job_attrs(multicohort),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=AnnotateFromGlobalCallset)
class SubsetMtToDatasetFromGlobalCallset(SubsetMtToDatasetWithHail):
    """SubsetMtToDatasetWithHail variant that subsets from AnnotateFromGlobalCallset.

    Only kicks in when the synthetic workflow runs a multi-dataset multicohort or has
    only_families set on a dataset (see AnnotateDataset._get_cohort_mt for the branch).
    Otherwise the downstream AnnotateDatasetFromGlobalCallset reads directly from
    AnnotateFromGlobalCallset and this stage is skipped.
    """

    def _get_cohort_mt(self, inputs: stage.StageInput) -> Path:
        return inputs.as_path(target=workflow.get_multicohort(), stage=AnnotateFromGlobalCallset)


@stage.stage(
    required_stages=[AnnotateFromGlobalCallset, SubsetMtToDatasetFromGlobalCallset],
    analysis_type='matrixtable',
)
class AnnotateDatasetFromGlobalCallset(AnnotateDataset):
    """AnnotateDataset variant that reads from the global-join stack.

    Single-dataset multicohorts with no only_families config read directly from
    AnnotateFromGlobalCallset. Multi-dataset or family-filtered runs fall through to
    SubsetMtToDatasetFromGlobalCallset, mirroring the base class's branch structure.
    """

    def _get_cohort_mt(self, dataset: targets.Dataset, inputs: stage.StageInput) -> Path:
        family_sgs = utils.get_family_sequencing_groups(dataset)
        if len(workflow.get_multicohort().get_datasets()) == 1 and family_sgs is None:
            return inputs.as_path(target=workflow.get_multicohort(), stage=AnnotateFromGlobalCallset)
        return inputs.as_path(target=dataset, stage=SubsetMtToDatasetFromGlobalCallset, key='mt')


@stage.stage(
    required_stages=[AnnotateDatasetFromGlobalCallset],
    analysis_type='es-index',
    analysis_keys=['done_flag'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'VARIANTS'},  # noqa: ARG005
)
class ExportMtAsEsIndexFromGlobalCallset(ExportMtAsEsIndex):
    """ExportMtAsEsIndex variant sourced from the global-join AnnotateDataset variant."""

    def _get_annotated_mt_path(self, dataset: targets.Dataset, inputs: stage.StageInput) -> str:
        return inputs.as_str(target=dataset, stage=AnnotateDatasetFromGlobalCallset)
