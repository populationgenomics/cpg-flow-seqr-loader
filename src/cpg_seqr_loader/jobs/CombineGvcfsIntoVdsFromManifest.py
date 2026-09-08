from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_combiner_from_manifest_jobs(
    manifest_path: Path,
    output_vds: Path,
    combiner_plan: Path,
    temp_dir_string: str,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """Build a fresh VDS from a pre-baked gVCF manifest.

    Unlike the metamist-driven combiner, this job takes the manifest path directly and skips
    all incremental-VDS machinery (existing-VDS discovery, refresh handling, sample removal).
    Each run rebuilds from scratch, which is fine at synthetic-cohort scale.
    """
    localised_manifest = hail_batch.get_batch().read_input(str(manifest_path))

    job = hail_batch.get_batch().new_bash_job(
        'CombineGvcfsIntoVdsFromManifest',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    # Non-spot for the same reason as the standard combiner: preemption has previously caused
    # multiple simultaneous QOB groups to race on the same data.
    job.spot(False)

    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.run_combiner \\
            --output_vds {output_vds!s} \\
            --plan {combiner_plan!s} \\
            --tmp {temp_dir_string!s} \\
            --gvcf_add_file {localised_manifest}
        """
    )

    return job
