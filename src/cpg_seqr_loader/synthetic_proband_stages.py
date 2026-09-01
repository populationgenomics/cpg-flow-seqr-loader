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

from cpg_flow import stage, targets
from cpg_utils import Path

from cpg_seqr_loader import utils
from cpg_seqr_loader.jobs.GenerateSyntheticProbandCombinerInputs import write_combiner_inputs
from cpg_seqr_loader.jobs.GenerateSyntheticProbandGvcfs import create_synthetic_gvcf_jobs


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
