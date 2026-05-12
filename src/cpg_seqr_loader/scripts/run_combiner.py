"""
Runs the hail combiner
This is copied from the large_cohort implementation, as RD needs to alter some defaults
Specifically we need to drop/modify the Force argument to resume from previous Combiner plan/runs
"""

import argparse
import logging

import loguru
from cpg_flow import utils
from cpg_utils import config, hail_batch, to_path

import hail as hl


def main(
    output_vds_path: str,
    tmp_prefix: str,
    combiner_plan: str | None,
    gvcfs_to_combine_file: str | None = None,
    vds_path: str | None = None,
    sgs_to_remove_file: str | None = None,
) -> None:
    """
    Runs the combiner -

    1. this method receives an existing VDS path or None - if None, combine from scratch. If populated, the combiner
       run will use that as a base, and add additional samples to it
    2. if this method receives sgs_to_remove and a vds_path, the first step is to create a new VDS from temp with the
       samples removed. This is then used as the base for the combiner run
    3. if this method receives gvcf_paths, the combiner will run with those as the input paths to add
    4. if there are no gVCFs, the VDS with samples removed is written to the output path directly

    Args:
        output_vds_path (str): eventual output path for the VDS
        tmp_prefix (str): where to store temporary combiner intermediates
        combiner_plan (str | None): where to store the combiner plan, or where to resume from
        gvcfs_to_combine_file (str | None): A file containing all paths to GVCFs we're combining, or None
        vds_path (str | None): a single VDS, or None - this is where a combiner can continue on from
        sgs_to_remove_file (str | None): A file containing AG IDs to remove from the combiner, or None
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),
        worker_cores=config.config_retrieve(['combiner', 'worker_cores']),
    )

    # allow force-new from config. This would force a fresh combining of the input gVCFs and prior VDSs as appropriate
    # it would not restart a completely fresh VDS from component gVCFs
    # this is relevant for instances where a previous combiner attempt started, but failed due to resourcing issues, and
    # we want to start fresh with new intervals/partition strategy
    force_new_combiner = config.config_retrieve(['combiner', 'force_new_combiner'])

    # Load from save, if supplied (log correctly depending on force_new_combiner)
    if combiner_plan and force_new_combiner:
        loguru.logger.info(f'Combiner plan {combiner_plan} will be ignored/written new')

    # combiner plan as an argument comes from the Stage, so it always has a real value, but the file may not exist.
    elif combiner_plan and to_path(combiner_plan).exists():
        loguru.logger.info(f'Resuming combiner plan from {combiner_plan}.')

    # do we have new gVCFs to add?
    if gvcfs_to_combine_file is not None:
        # if gvcfs_to_combine is a file, read it and split into a list
        with open(gvcfs_to_combine_file) as f:
            gvcf_paths = [line.strip() for line in f if line.strip()]
    else:
        gvcf_paths = None

    # do we have any SGs to remove? i.e. scenarios where we have to remove SGs due to retraction
    if sgs_to_remove_file is not None:
        # if gvcfs_to_combine is a file, read it and split into a list
        with open(sgs_to_remove_file) as f:
            sgs_to_remove = [line.strip() for line in f if line.strip()]
            num_sgs_to_remove = len(sgs_to_remove)

            # quick check - if we're removing many samples, confirm intent
            if num_sgs_to_remove >= config.config_retrieve(['combiner', 'sg_remove_threshold']):
                if not config.config_retrieve(['combiner', 'sg_remove_confirm']):
                    raise ValueError(
                        f'Attempting to remove {num_sgs_to_remove} samples from {gvcf_paths}, '
                        'if this is intentional set the config parameter combiner.sg_remove_confirm to proceed'
                    )
                loguru.logger.warning(f'Removing {num_sgs_to_remove} samples from {vds_path}')
    else:
        sgs_to_remove = None

    # if neither of these are true, we don't need to be here - the Stage should not be running.
    if not (sgs_to_remove or gvcf_paths):
        raise ValueError('No samples to remove or gVCFs to add - please provide at least one of these')

    # logical steps -
    # 1. if there are samples to remove, do that first
    #   a. if there are no samples to add, just remove the samples and write to eventual output path
    #   b. if there are samples to add, write this to a temporary path, set force to true, and then run the combiner
    # 2. if there are samples to add, run the combiner with the final path as the output path

    # 1 - do we need to do removal?
    if vds_path and sgs_to_remove:
        loguru.logger.info(f'Removing sample groups {sgs_to_remove} from {vds_path}')

        temp_path = f'{tmp_prefix}/combiner_removal_temp.vds'

        # this will only exist if the previous removal was successful, AND we have additional gVCFs to add
        if utils.exists(temp_path):
            loguru.logger.info(f'Found existing VDS at {temp_path}, skipping removal step')
            vds_path = temp_path

        else:
            vds = hl.vds.read_vds(vds_path)
            vds = hl.vds.filter_samples(vds, samples=sgs_to_remove, keep=False, remove_dead_alleles=True)

            # 1a. if there are no samples to add, just remove the samples and write to eventual output path
            if not gvcf_paths:
                loguru.logger.info(f'Writing to {output_vds_path}')
                vds.write(output_vds_path)
                return

            # 1b. if there are samples to add, write this to a temporary path, set force to true, then run the combiner
            # we've changed the VDS base, so we need to force a new combiner plan
            force_new_combiner = True

            loguru.logger.info(f'Writing with removed SGs to {temp_path}')
            vds.write(temp_path)
            vds_path = temp_path

    # final round of combiner arguments from config
    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])

    # for Combining phase, aim for a high number of partitions, ~5k
    vds_intervals_path = config.config_retrieve(['combiner', 'vds_intervals', sequencing_type])

    if not vds_intervals_path:
        raise ValueError(f'Provided path for VDS intervals: {vds_intervals_path} - please provide a real path.')

    vds_intervals = hl.import_bed(vds_intervals_path, reference_genome='GRCh38').interval.collect()

    # 2 - do we need to run the combiner?
    hl.vds.new_combiner(
        output_path=output_vds_path,
        save_path=combiner_plan,
        gvcf_paths=gvcf_paths,
        vds_paths=[vds_path] if vds_path else None,
        reference_genome='GRCh38',
        temp_path=tmp_prefix,
        force=force_new_combiner,
        intervals=vds_intervals,
        branch_factor=config.config_retrieve(['combiner', 'branch_factor']),
        target_records=config.config_retrieve(
            ['combiner', 'target_records']
        ),  # keep this high to prevent hail overriding preset partitions
        gvcf_batch_size=config.config_retrieve(['combiner', 'gvcf_batch_size']),
    ).run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_vds', required=True, help='Location for the new combiner product')
    parser.add_argument('--input_vds', help='Path to a VDS this run will build on top of')
    parser.add_argument('--plan', help='Path to a Combiner plan')
    parser.add_argument('--tmp', help='Temporary directory to use', required=True)
    parser.add_argument('--gvcf_add_file', help='File containing gVCF paths to combine')
    parser.add_argument('--sg_remove_file', help='File containing SG IDs to remove from the input VDS')
    args = parser.parse_args()

    # boot up the logging module so Hail will report on combiner progress
    logging.basicConfig(level=logging.INFO)

    main(
        output_vds_path=args.output_vds,
        vds_path=args.input_vds,
        combiner_plan=args.plan,
        tmp_prefix=args.tmp,
        gvcfs_to_combine_file=args.gvcf_add_file,
        sgs_to_remove_file=args.sg_remove_file,
    )
