"""
pull out a single Dataset's samples from a MatrixTable
"""

from argparse import ArgumentParser

import loguru
from cpg_utils import config, hail_batch, to_path

import hail as hl


def subset_mt_to_samples(input_mt: str, sg_id_file: str, output: str):
    """
    Subset the MatrixTable to the provided list of samples and to variants present in those samples

    Args:
        input_mt (str): cohort-level matrix table from VCF.
        sg_id_file (list[str]): samples to take from the matrix table.
        output (str): path to write the result.
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),
    )

    mt = hl.read_matrix_table(input_mt)

    id_file = to_path(sg_id_file)
    if not id_file.exists():
        raise ValueError(f'Sample ID file {id_file} does not exist')

    id_list = {each.strip() for each in id_file.read_text().splitlines()}

    mt_sample_ids = set(mt.s.collect())

    if sample_ids_not_in_mt := id_list - mt_sample_ids:
        raise ValueError(
            f'Found {len(sample_ids_not_in_mt)}/{len(id_list)} IDs in the requested subset not in the callset.\n'
            f"IDs that aren't in the callset: {sample_ids_not_in_mt}\n"
            f'All callset sample IDs: {mt_sample_ids}',
        )

    loguru.logger.info(f'Found {len(mt_sample_ids)} samples in mt, subsetting to {len(id_list)} samples.')

    mt = mt.filter_cols(hl.literal(id_list).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt.write(output, overwrite=True)
    loguru.logger.info(f'Finished subsetting to {len(id_list)} samples, written to {output}.')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input MatrixTable to subset')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    parser.add_argument('--sg_id_file', required=True, help='Samples to include in the subset')
    args = parser.parse_args()
    subset_mt_to_samples(
        input_mt=args.input,
        sg_id_file=args.sg_id_file,
        output=args.output,
    )


if __name__ == '__main__':
    cli_main()
