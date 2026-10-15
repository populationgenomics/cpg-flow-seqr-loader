import json
from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_seqr_vcf_job(
    input_mt: str,
    vqsr_vcf: str,
    sg_to_participant: dict[str, str],
    id_map_file: Path,
    output_vcf: Path,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """Write the SG -> participant ID mapping, and create a job to write a seqr VCF from the dense MT."""

    if not config.config_retrieve(['workflow', 'dry_run'], False):
        with id_map_file.open('w') as f:
            json.dump(sg_to_participant, f, indent=2)

    job = hail_batch.get_batch().new_bash_job('Seqr VCF from dense MT', attributes=job_attrs | {'tool': 'hail'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.spot(False)
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.seqr_vcf_from_dense_mt \\
            --input {input_mt} \\
            --vqsr {vqsr_vcf} \\
            --id_map {id_map_file!s} \\
            --output {output_vcf!s}
        """
    )
    return job
