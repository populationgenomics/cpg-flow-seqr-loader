from argparse import ArgumentParser

from cpg_utils import config, hail_batch
from loguru import logger

import hail as hl


def vcf_from_mt_subset(input_mt: str, output: str):
    """
    Read the MT in, and write out to a VCF
    If we wanted to translate sample IDs to external samples
    then we could do that here, otherwise rely on VCF re-heading

    Args:
        input_mt (str): path of the single-dataset MT to read in
        output (str): path of the vcf.bgz to generate
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['vcf_from_mt', 'worker_memory']),
        driver_memory=config.config_retrieve(['vcf_from_mt', 'driver_memory']),
        driver_cores=config.config_retrieve(['vcf_from_mt', 'driver_cores']),
    )

    mt = hl.read_matrix_table(input_mt)
    logger.info(f'Dataset MT dimensions: {mt.count()}')

    # filter to this sample's non-ref calls
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    # drop gvcf_info - a dict, can't sit in FORMAT
    # only relevant when doing gVCF -> VDS -> MT
    if 'gvcf_info' in mt.row:
        mt = mt.drop('gvcf_info')
    mt = mt.drop('variant_qc')

    hl.export_vcf(mt, out_path, tabix=True)
    logger.info(f'Written {output}')


if __name__ == '__main__':
    parser = ArgumentParser(description='Write a VCF from a single-dataset MT')
    parser.add_argument('--input', required=True, help='Path to the single-dataset MT to read in')
    parser.add_argument('--output', required=True, help='Path to write the VCF to')
    args = parser.parse_args()

vcf_from_mt_subset(input_mt=args.input, output=args.output)
