"""
takes a path to a VDS, and an output path
densifies the VDS into a MatrixTable
Writes the MT to the output path

Additional arguments:
    --sites_only: optional, if used write a sites-only VCF directory to this location
                  it's far more efficient for Hail to write out a VCF per-partition, rather than a single VCF
                  this also prevents the need for writing out a single VCF, then fragmenting that to run VEP, VQSR
    --separate_header: optional, if used write a sites-only VCF directory with a separate header to this location
                  this makes the VCF write more efficient, and is easily rejoined into a single VCF using compose
                  for situations where we need a single whole-genome/exome VCF (e.g. VQSR training)
"""

import argparse

import hail as hl
import loguru
from cpg_flow import utils
from cpg_utils import config, hail_batch
from cpg_seqr_loader.hail_scripts.sparse_mt import default_compute_info
from cpg_seqr_loader.hail_scripts.vcf import adjust_vcf_incompatible_types


def main(
    vds_in: str,
    dense_mt_out: str,
    sites_only: str | None = None,
    separate_header: str | None = None,
) -> None:
    """
    Load a sparse VariantDataset
    Write a dense MatrixTable with split multiallelics
    optionally, coalesce the data into a predetermined number of partitions
    - defaulting to the most computationally efficient coalesce, instead of a full repartition
    - option to override this with a naive coalesce, or a shuffle repartition

    Args:
        vds_in (str):
        dense_mt_out (str):
        sites_only (str): optional, if used write a sites-only VCF directory to this location
        separate_header (str): optional, write a sites-only VCF directory with a separate header to this location
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory'], 'highmem'),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory'], 'highmem'),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores'], 2),
    )

    partitions = config.config_retrieve(['workflow', 'densify_partitions'], 2000)

    # check here to see if we can reuse the dense MT
    if not utils.can_reuse(dense_mt_out):
        loguru.logger.info(f'Densifying data, using {partitions} partitions')

        # providing n_partitions here gets Hail to calculate the intervals per partition on the VDS var and ref data
        vds = hl.vds.read_vds(vds_in, n_partitions=partitions)

        mt = hl.vds.to_dense_mt(vds)

        # taken from _filter_rows_and_add_tags in large_cohort/site_only_vcf.py
        # remove any monoallelic or non-ref-in-any-sample sites
        mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

        # annotate site-level DP to avoid name collision
        mt = mt.annotate_rows(
            site_dp=hl.agg.sum(mt.DP),
            ANS=hl.agg.count_where(hl.is_defined(mt.LGT)) * 2,
        )

        # content shared with large_cohort.site_only_vcf.py
        info_ht = default_compute_info(mt, site_annotations=True, n_partitions=mt.n_partitions())
        info_ht = info_ht.annotate(info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp))

        info_ht = adjust_vcf_incompatible_types(
            info_ht,
            # with default INFO_VCF_AS_PIPE_DELIMITED_FIELDS, AS_VarDP will be converted
            # into a pipe-delimited value e.g.: VarDP=|132.1|140.2
            # which breaks VQSR parser (it doesn't recognise the delimiter and treats
            # it as an array with a single string value "|132.1|140.2", leading to
            # an IndexOutOfBound exception when trying to access value for second allele)
            pipe_delimited_annotations=[],
        )

        # annotate this info back into the main MatrixTable
        mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

        # unpack mt.info.info back into mt.info. Must be better syntax for this?
        mt = mt.drop('gvcf_info')

        loguru.logger.info('Splitting multiallelics, in a sparse way')
        mt = hl.experimental.sparse_split_multi(mt)

        loguru.logger.info(f'Writing fresh data into {dense_mt_out}')
        mt.write(dense_mt_out, overwrite=True)

    else:
        loguru.logger.info(f'Accepting existing data in {dense_mt_out}')

    if not (sites_only or separate_header):
        return

    # read the dense MT and obtain the sites-only HT
    mt = hl.read_matrix_table(dense_mt_out)

    sites_only_ht = mt.rows()

    # write a directory containing all the per-partition VCF fragments, each with a VCF header
    if sites_only:
        loguru.logger.info('Writing sites-only VCF, header-per-shard')
        hl.export_vcf(sites_only_ht, sites_only, tabix=True, parallel='header_per_shard')

    # write a directory containing all the per-partition VCF fragments, with a separate VCF header file
    if separate_header:
        loguru.logger.info('Writing sites-only VCF, separate-header')
        hl.export_vcf(sites_only_ht, separate_header, tabix=True, parallel='separate_header')


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        help='Path to the input VDS',
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Path to write the result',
        required=True,
    )
    parser.add_argument(
        '--sites_only',
        help='Specify an output path for a sites-only VCF, or None',
    )
    parser.add_argument(
        '--separate_header',
        help='Specify an output path for a sites-only VCF, with a separate header file, or None',
    )
    args = parser.parse_args()
    main(
        vds_in=args.input,
        dense_mt_out=args.output,
        sites_only=args.sites_only,
        separate_header=args.separate_header,
    )


if __name__ == '__main__':
    cli_main()
