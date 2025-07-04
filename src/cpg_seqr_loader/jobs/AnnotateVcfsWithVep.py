"""
Creates a Hail Batch job to run the command line VEP tool.
"""

from typing import TYPE_CHECKING

import hailtop.batch as hb

from cpg_utils import Path, to_path, config, hail_batch
from cpg_seqr_loader import utils
from cpg_flow import utils as cpg_flow_utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def add_vep_jobs(
    manifest_file: Path,
    tmp_prefix: Path,
    final_out_path: Path,
    job_attrs: dict,
) -> list['BashJob']:
    """
    Runs VEP on provided VCF fragments. Writes annotations per-fragment as JSON output.
    """

    input_vcfs = utils.get_all_fragments_from_manifest(manifest_file)

    jobs = []

    fragment_count = len(input_vcfs)

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'
    result_part_paths = []
    for idx, resource in enumerate(input_vcfs):
        result_part_path = result_parts_bucket / f'part{idx + 1}.jsonl'

        result_part_paths.append(result_part_path)

        if cpg_flow_utils.can_reuse(result_part_path):
            continue

        jobs.append(
            vep_one(
                vcf=resource[utils.VCF_GZ],
                out_path=str(result_part_paths[idx]),
                job_attrs=job_attrs | {'part': f'{idx + 1}/{fragment_count}'},
            ),
        )

    j = gather_vep_json_to_ht(
        vep_results_paths=result_part_paths,
        out_path=final_out_path,
        job_attrs=job_attrs,
    )
    j.depends_on(*jobs)
    jobs.append(j)
    return jobs


def gather_vep_json_to_ht(
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict,
) -> 'BashJob':
    """
    Parse results from VEP with annotations formatted in JSON, and write into a Hail Table using a Batch job.
    """

    job = hail_batch.get_batch().new_job('VEP', attributes=job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    json_paths = ' '.join(str(p) for p in vep_results_paths)
    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.vep_json_to_ht \\
            --input {json_paths} \\
            --output {out_path}
        """
    )
    return job


def vep_one(
    vcf: hb.ResourceFile,
    out_path: str,
    job_attrs: dict,
) -> 'BashJob':
    """Run a single VEP job."""

    # check that the cache and image for this version exist
    vep_image = config.config_retrieve(['images', 'vep'])
    vep_mount_path = to_path(config.config_retrieve(['references', 'vep_mount']))

    job = hail_batch.get_batch().new_bash_job('AnnotateFragmentedVcfWithVep', job_attrs | {'tool': 'VEP'})
    job.image(vep_image)

    # vep is single threaded, with a middling memory requirement
    # during test it can exceed 8GB, so we'll give it 16GB
    job.memory('16Gi').storage('15Gi').cpu(1)

    # gcsfuse works only with the root bucket, without prefix:
    data_mount = to_path(f'/{vep_mount_path.drive}')
    job.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

    loftee_conf = {
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
        'loftee_path': '$VEP_DIR_PLUGINS',
    }

    # sexy new plugin - only present in 110 build
    alpha_missense_plugin = f'--plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz '

    job.command(
        f"""
    set -x

    vep \\
        --format vcf \\
        --json \\
        -o {job.output} \\
        -i {vcf} \\
        --everything \\
        --mane_select \\
        --allele_number \\
        --minimal \\
        --species homo_sapiens \\
        --cache --offline --assembly GRCh38 \\
        --dir_cache {vep_dir}/vep/ \\
        --fasta /cpg-common-main/references/vep/110/mount/vep/homo_sapiens/110/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \\
        {alpha_missense_plugin} \\
        --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())} \\
        --plugin UTRAnnotator,file=$UTR38
    """,
    )

    hail_batch.get_batch().write_output(job.output, out_path)

    return job
