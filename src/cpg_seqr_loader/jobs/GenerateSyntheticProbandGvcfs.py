"""
Job factories for the synthetic proband gVCF stage.

Two factories work together:

  - create_synthesis_jobs: one BashJob per duo family that invokes
    scripts/create_synthetic_proband_gvcf.py with the family's parental gVCFs, writing the
    output gVCF to the durable per-family path supplied by the stage.

  - create_analysis_registration_jobs: one BashJob per duo family that invokes
    scripts/register_synthetic_gvcf_analysis.py to record (or refresh) the synthetic gVCF as a
    metamist Analysis of type SYNTHETIC_GVCF_ANALYSIS_TYPE. Each registration job depends_on
    the matching synthesis job so registration only fires once the gVCF actually exists.

The synthesis script accepts gs:// paths directly (it localises inputs itself), so we don't
need Batch's read_input_group here.
"""

from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

from cpg_seqr_loader.utils import SyntheticProbandFamily

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_synthesis_jobs(
    families: list[SyntheticProbandFamily],
    output_paths: dict[str, Path],
    job_attrs: dict,
) -> dict[str, 'BashJob']:
    """Build one Batch synthesis job per qualifying duo family.

    Args:
        families: enumerated by utils.get_families_for_synthetic_probands.
        output_paths: keyed by family_id; the durable output gVCF path per family.
        job_attrs: base job attributes from the calling stage (dataset, workflow, etc.).

    Returns:
        Dict keyed by family_id -> BashJob, so downstream registration jobs can depend on the
        matching synthesis job. Empty dict if no families qualify.
    """
    jobs: dict[str, 'BashJob'] = {}
    for family in families:
        output = output_paths[family.family_id]

        job = hail_batch.get_batch().new_bash_job(
            f'SynthesizeProband_{family.family_id}',
            attributes=job_attrs | {'tool': 'pysam', 'family_id': family.family_id},
        )
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        # Sized for a WGS trio synthesis: 2 parental gVCFs (~20 GB each) + one output (~20 GB)
        # + slack. Overridable via config so exome / smaller-genome runs can shrink this.
        job.storage(config.config_retrieve(['synthetic_proband', 'storage'], '80G'))
        job.memory(config.config_retrieve(['synthetic_proband', 'memory'], 'standard'))

        job.command(
            f"""
            python -m cpg_seqr_loader.scripts.create_synthetic_proband_gvcf \\
                --mother_gvcf {family.mother_sg.gvcf!s} \\
                --father_gvcf {family.father_sg.gvcf!s} \\
                --sample_name {family.synthetic_sample_name} \\
                --output {output!s}
            """,
        )

        jobs[family.family_id] = job

    return jobs


def create_analysis_registration_jobs(
    families: list[SyntheticProbandFamily],
    output_paths: dict[str, Path],
    synthesis_jobs: dict[str, 'BashJob'],
    project: str,
    job_attrs: dict,
) -> list['BashJob']:
    """Build one Batch job per family that registers (or refreshes) the metamist Analysis.

    Each returned job depends on the matching synthesis job so registration only fires once the
    synthetic gVCF actually exists. If the framework has skipped a synthesis job because the
    output already exists, the corresponding key won't be present in synthesis_jobs; in that
    case the registration job simply has no upstream dependency and runs immediately (the
    script itself is idempotent).

    Args:
        families: same list passed to create_synthesis_jobs.
        output_paths: same dict passed to create_synthesis_jobs.
        synthesis_jobs: return value of create_synthesis_jobs (family_id -> BashJob).
        project: metamist project name (the dataset the parental SGs live in).
        job_attrs: base job attributes from the calling stage.
    """
    jobs: list['BashJob'] = []
    for family in families:
        output = output_paths[family.family_id]

        job = hail_batch.get_batch().new_bash_job(
            f'RegisterSyntheticProband_{family.family_id}',
            attributes=job_attrs | {'tool': 'metamist', 'family_id': family.family_id},
        )
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        upstream = synthesis_jobs.get(family.family_id)
        if upstream is not None:
            job.depends_on(upstream)

        job.command(
            f"""
            python -m cpg_seqr_loader.scripts.register_synthetic_gvcf_analysis \\
                --project {project} \\
                --gvcf_path {output!s} \\
                --sample_name {family.synthetic_sample_name} \\
                --family_id {family.family_id} \\
                --mother_sg_id {family.mother_sg.id} \\
                --father_sg_id {family.father_sg.id} \\
                --mother_source_gvcf {family.mother_sg.gvcf!s} \\
                --father_source_gvcf {family.father_sg.gvcf!s}
            """,
        )

        jobs.append(job)

    return jobs
