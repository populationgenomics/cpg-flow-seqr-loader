"""
Read in a MT, and re-jig the annotations ready for Seqr Export
"""

from argparse import ArgumentParser

from cpg_utils import config, hail_batch
from loguru import logger

import hail as hl


def annotate_dataset_mt(mt_path: str, out_mt_path: str) -> None:
    """
    Add dataset-level annotations.
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['annotate_dataset', 'worker_memory']),
        driver_memory=config.config_retrieve(['annotate_dataset', 'driver_memory']),
        driver_cores=config.config_retrieve(['annotate_dataset', 'driver_cores']),
    )

    mt = hl.read_matrix_table(mt_path)

    # Convert the mt genotype entries into num_alt, gq, ab, dp, and sample_id.
    is_called = hl.is_defined(mt.GT)
    genotype_fields = {
        'num_alt': hl.if_else(is_called, mt.GT.n_alt_alleles(), -1),
        'gq': hl.if_else(is_called, mt.GQ, hl.null(hl.tint)),
        'ab': hl.bind(
            lambda total: hl.if_else(
                is_called & (total != 0) & (hl.len(mt.AD) > 1),
                hl.float(mt.AD[1] / total),
                hl.missing(hl.tfloat),
            ),
            hl.sum(mt.AD),
        ),
        'ad': hl.if_else(is_called, hl.delimit(mt.AD, ','), hl.missing(hl.tstr)),
        'dp': hl.if_else(is_called, hl.int(hl.min(mt.DP, 32000)), hl.missing(hl.tfloat)),
        'sample_id': mt.s,
    }
    logger.info('Annotating genotypes')
    mt = mt.annotate_rows(genotypes=hl.agg.collect(hl.struct(**genotype_fields)))

    def _genotype_filter_samples(fn) -> hl.expr.SetExpression:
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    # top level - decorator
    def _capture_i_decorator(func):  # noqa: ANN202
        # call the returned_function(i) which locks in the value of i
        def _inner_filter(i):  # noqa: ANN202
            # the _genotype_filter_samples will call this _func with g
            def _func(g):  # noqa: ANN202
                return func(i, g)

            return _func

        return _inner_filter

    @_capture_i_decorator
    def _filter_num_alt(i, g) -> bool:
        return i == g.num_alt

    @_capture_i_decorator
    def _filter_samples_gq(i, g) -> hl.expr.BooleanExpression:
        return (g.gq >= i) & (g.gq < i + 5)

    @_capture_i_decorator
    def _filter_samples_ab(i, g) -> hl.expr.BooleanExpression:
        return (g.num_alt == 1) & ((g.ab * 100) >= i) & ((g.ab * 100) < i + 5)

    mt = mt.annotate_rows(
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{str(i): _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq=hl.struct(
            **{f'{i}_to_{i + 5}': _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 95, 5)},
        ),
        samples_ab=hl.struct(
            **{f'{i}_to_{i + 5}': _genotype_filter_samples(_filter_samples_ab(i)) for i in range(0, 45, 5)},
        ),
    )
    mt.write(out_mt_path, overwrite=True)
    logger.info(f'Written {out_mt_path}')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input MatrixTable to subset')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    args = parser.parse_args()
    annotate_dataset_mt(
        mt_path=args.input,
        out_mt_path=args.output,
    )


if __name__ == '__main__':
    cli_main()
