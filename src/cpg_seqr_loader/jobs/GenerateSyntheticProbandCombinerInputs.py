"""
Driver-side writer for the synthetic proband combiner-inputs stage.

The single public entry point `write_combiner_inputs` produces both of the stage's output
artifacts from the driver, no Batch jobs required:

  - Pedigree PED file: 6-column TSV with three rows per qualifying duo family (mother, father,
    synthetic proband), plus a row for every real SG in the multicohort.

  - gVCF manifest: newline-separated list of every gVCF that will go into the synthetic-trio
    combiner run. It includes every real SG in the multicohort that has a gVCF (qualifying or
    not, so we don't silently drop samples from the seqr load), plus every synthetic gVCF
    from Stage 1's outputs.

Both files are tiny text blobs, so the previous approach of spinning up a Batch job per file
(only to run `cat > $output <<HEREDOC`) added VM + container startup cost and a heredoc-quoting
hazard for no benefit. Writing directly to `gs://` via `cpg_utils.Path.open('w')` is the same
pattern used elsewhere in the repo (see jobs/SubsetMtToDatasetWithHail.py).

The write is skipped in dry-run mode - unlike Batch jobs, which cpg_flow describes rather than
executes under `workflow.dry_run`, a driver-side `.open('w')` executes unconditionally, so the
callers must gate on that config themselves.
"""

from cpg_flow import targets
from cpg_utils import Path, config

from cpg_seqr_loader.utils import (
    SyntheticProbandFamily,
    build_gvcf_manifest_content,
    build_synthetic_pedigree_content,
)


def write_combiner_inputs(
    families: list[SyntheticProbandFamily],
    multicohort: targets.MultiCohort,
    synthetic_gvcf_paths: dict[str, Path],
    output_ped: Path,
    output_manifest: Path,
) -> None:
    """Write the PED and gVCF manifest for the synthetic-trio combiner run.

    No Batch jobs are queued - both files are written directly from the driver. The stage's
    `queue_jobs` still returns an empty job list so cpg_flow's REUSE check keys off the output
    paths existing on disk, exactly as it would with jobs.
    """
    if config.config_retrieve(['workflow', 'dry_run'], False):
        return

    pedigree_content = build_synthetic_pedigree_content(multicohort, families)
    with output_ped.open('w') as f:
        f.write(pedigree_content)

    real_paths = [str(sg.gvcf) for sg in multicohort.get_sequencing_groups() if sg.gvcf is not None]
    synthetic_paths = [str(synthetic_gvcf_paths[family.family_id]) for family in families]
    manifest_content = build_gvcf_manifest_content(real_paths + synthetic_paths)
    with output_manifest.open('w') as f:
        f.write(manifest_content)
