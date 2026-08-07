"""
Job factories for the synthetic proband combiner-inputs stage.

Two artifacts, one Batch job each:

  - Pedigree PED file: 6-column TSV with three rows per qualifying duo family (mother, father,
    synthetic proband). Persistent so seqr sync can reuse it across loads.

  - gVCF manifest: newline-separated list of every gVCF that will go into the synthetic-trio
    combiner run (the ravenscroft-rpl invocation, in the initial rollout) - every real SG in the
    multicohort that has a gVCF (qualifying or not, so we don't silently drop samples from the
    seqr load), plus every synthetic gVCF from Stage 1's outputs.

Content for each file is computed in the driver via the pure builders in `utils`, then dropped
into a Batch job via a quoted heredoc; `write_output` uploads the produced file to the durable
gs:// path once the job succeeds. This keeps the writes gated on Batch success and gives the
framework a real job to depend on rather than a driver-side side effect.
"""

from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

from cpg_seqr_loader.utils import (
    SyntheticProbandFamily,
    build_gvcf_manifest_content,
    build_synthetic_pedigree_content,
)

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def _write_content_job(name: str, content: str, out_path: Path, job_attrs: dict) -> 'BashJob':
    """Return a Batch job that writes a fixed string to a durable output path.

    Uses a quoted heredoc (`<<'EOF'`) so tabs / newlines pass through literally with no shell
    expansion. Safe for our alphanumeric family_id / SG_id content; would need escaping for
    arbitrary user input.
    """
    job = hail_batch.get_batch().new_bash_job(name, attributes=job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f"cat > {job.output} <<'STAGE2_EOF'\n{content}STAGE2_EOF\n")
    hail_batch.get_batch().write_output(job.output, str(out_path))
    return job


def create_pedigree_job(
    families: list[SyntheticProbandFamily],
    output_ped: Path,
    job_attrs: dict,
) -> 'BashJob':
    """One Batch job that writes the synthetic-proband PED for all qualifying duo families."""
    content = build_synthetic_pedigree_content(families)
    return _write_content_job(
        name='WriteSyntheticPedigree',
        content=content,
        out_path=output_ped,
        job_attrs=job_attrs | {'tool': 'python'},
    )


def create_manifest_job(
    families: list[SyntheticProbandFamily],
    synthetic_gvcf_paths: dict[str, Path],
    multicohort: targets.MultiCohort,
    output_manifest: Path,
    job_attrs: dict,
) -> 'BashJob':
    """One Batch job that writes the combiner gVCF manifest.

    Manifest composition (deliberately inclusive so samples aren't silently dropped from seqr):
      - every real SG in the multicohort whose `.gvcf` is set, whether or not their family
        qualified for synthetic proband synthesis;
      - every synthetic gVCF from Stage 1's outputs (indexed by family_id).

    Real SGs whose `.gvcf` is None are logged upstream by cpg_flow and simply won't appear here
    - there's no useful path to write for them.
    """
    real_paths: list[str] = []
    for sg in multicohort.get_sequencing_groups():
        if sg.gvcf is None:
            continue
        real_paths.append(str(sg.gvcf))

    synthetic_paths = [str(synthetic_gvcf_paths[f.family_id]) for f in families]

    content = build_gvcf_manifest_content(real_paths + synthetic_paths)
    return _write_content_job(
        name='WriteGvcfManifestWithSynthetics',
        content=content,
        out_path=output_manifest,
        job_attrs=job_attrs | {'tool': 'python'},
    )
