from typing import TYPE_CHECKING

import loguru

from cpg_utils import hail_batch, config, Path
from cpg_seqr_loader import utils

from cpg_flow import targets, utils as cpg_flow_utils


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_combiner_jobs(
    multicohort: targets.MultiCohort,
    output_vds: Path,
    combiner_plan: Path,
    temp_dir: Path,
    job_attrs: dict[str, str],
) -> 'BashJob | None':
    vds_path: str | None = None
    sg_ids_in_vds: set[str] = set()
    sgs_to_remove: list[str] = []

    # check for a VDS by ID - this is not the typical RD process
    if vds_id := config.config_retrieve(['workflow', 'use_specific_vds'], None):
        vds_result_or_none = utils.query_for_specific_vds(vds_id)
        if vds_result_or_none is None:
            raise ValueError(f'Specified VDS ID {vds_id} not found in Metamist')

        # if not none, unpack the result
        vds_path, sg_ids_in_vds = vds_result_or_none

    # check for existing VDS by getting all and fetching latest
    elif config.config_retrieve(['workflow', 'check_for_existing_vds']):
        loguru.logger.info('Checking for existing VDS')
        if existing_vds_analysis_entry := utils.query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
            vds_path = existing_vds_analysis_entry['output']
            sg_ids_in_vds = {sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']}

    else:
        loguru.logger.info('Not continuing from any previous VDS, creating new Combiner from gVCFs only')

    # quick check - if we found a VDS, guarantee it exists
    if vds_path and not cpg_flow_utils.exists(vds_path):
        raise ValueError(f'VDS {vds_path} does not exist, but has an Analysis Entry')

    # this is a quick and confident check on current VDS contents, but does require a direct connection to the VDS
    # by default this is True, and can be deactivated in config
    if vds_path and config.config_retrieve(['workflow', 'manually_check_vds_sg_ids']):
        sg_ids_in_vds = utils.manually_find_ids_from_vds(vds_path)

    new_sg_gvcfs = [
        str(sg.gvcf)
        for sg in multicohort.get_sequencing_groups()
        if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
    ]

    # final check - if we have a VDS, and we have a current MultiCohort
    # detect any samples which should be _removed_ from the current VDS prior to further combining taking place
    sg_remove_arg = ''
    if sg_ids_in_vds:
        sgs_in_mc: list[str] = multicohort.get_sequencing_group_ids()
        loguru.logger.info(f'Found {len(sg_ids_in_vds)} SG IDs in VDS {vds_path}')
        loguru.logger.info(f'Total {len(sgs_in_mc)} SGs in this MultiCohort')

        sgs_to_remove = sorted(set(sg_ids_in_vds) - set(sgs_in_mc))

        if sgs_to_remove:
            loguru.logger.info(f'Removing {len(sgs_to_remove)} SGs from VDS {vds_path}')
            loguru.logger.info(f'SGs to remove: {sgs_to_remove}')

            # make a temp file containing these IDs, and write them to it

            # write the gVCF paths into a temporary file
            sg_remove_file = temp_dir / 'sgs_to_remove.txt'
            with sg_remove_file.open('w') as write_handle:
                for sgid in sgs_to_remove:
                    write_handle.write(f'{sgid!s}\n')
            sg_remove_arg = f'--sg_remove_file {sg_remove_file!s}'

    if not (new_sg_gvcfs or sgs_to_remove):
        loguru.logger.info('No GVCFs to add to, or remove from, existing VDS')
        loguru.logger.info(f'Checking if VDS exists: {output_vds}: {output_vds.exists()}')
        return None

    gvcf_add_arg = ''
    if new_sg_gvcfs:
        # write the gVCF paths into a temporary file
        gvcf_path_file = temp_dir / 'gvcfs_to_combine.txt'
        with gvcf_path_file.open('w') as write_handle:
            for gvcf_path in new_sg_gvcfs:
                write_handle.write(f'{gvcf_path!s}\n')
        gvcf_add_arg = f'--gvcf_add_file {gvcf_path_file!s}'

    job = hail_batch.get_batch().new_bash_job(
        'CombineGvcfsIntoVds',
        attributes=job_attrs | {'tool': 'hail'},
    )
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    # set this job to be non-spot (i.e. non-preemptible)
    # previous issues with preemptible VMs led to multiple simultaneous QOB groups processing the same data
    job.spot(config.config_retrieve(['workflow', 'preemptible_vms']))

    input_vds_arg = f'--input_vds {vds_path!s}' if vds_path else ''

    job.command(
        f"""
        python -m cpg_seqr_loader.scripts.run_combiner \\
            --output_vds {output_vds!s} \\
            --plan {combiner_plan!s} \\
            --tmp {temp_dir / 'temp_dir'!s} {input_vds_arg} {gvcf_add_arg} {sg_remove_arg}
        """
    )

    return job
