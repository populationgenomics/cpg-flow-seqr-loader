"""
Read in a MT, and write it out as a single long-format Parquet file, one row per variant x sample

Query-on-Batch has no Spark, so we can't use to_spark().write.parquet(). Instead, chunks of partitions are collected
to this process and streamed into a single Parquet file with pyarrow, one or more row groups per chunk. The MT is
sorted by locus, and contig/position are top-level columns, so the row group min/max statistics act as an index of
the genomic interval covered by each row group.

Row fields are only collected once per variant (via localize_entries), and expanded to one row per sample in arrow.
"""

from argparse import ArgumentParser
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from loguru import logger

import hail as hl

from cpg_utils import config, hail_batch

# row annotations added in AnnotateDataset which summarise all samples per variant. These are redundant in a
# long-format table, and would be repeated for every sample on every row
REDUNDANT_ROW_FIELDS = ['genotypes', 'samples_no_call', 'samples_num_alt', 'samples_gq', 'samples_ab']

ENTRIES_FIELD = 'mt_entries'
COLS_FIELD = 'mt_cols'

PRIMITIVE_TYPES = {
    hl.tstr: pa.string(),
    hl.tint32: pa.int32(),
    hl.tint64: pa.int64(),
    hl.tfloat32: pa.float32(),
    hl.tfloat64: pa.float64(),
    hl.tbool: pa.bool_(),
}


def hail_to_arrow_type(dtype: hl.HailType) -> pa.DataType:
    """Map a Hail type to an equivalent Arrow type. Calls, Loci, and Intervals are represented as strings."""
    if dtype in PRIMITIVE_TYPES:
        return PRIMITIVE_TYPES[dtype]
    if dtype == hl.tcall or isinstance(dtype, hl.tlocus | hl.tinterval):
        return pa.string()
    if isinstance(dtype, hl.tarray | hl.tset):
        return pa.list_(hail_to_arrow_type(dtype.element_type))
    if isinstance(dtype, hl.tdict):
        return pa.map_(hail_to_arrow_type(dtype.key_type), hail_to_arrow_type(dtype.value_type))
    if isinstance(dtype, hl.tstruct):
        return pa.struct([pa.field(name, hail_to_arrow_type(field_type)) for name, field_type in dtype.items()])
    if isinstance(dtype, hl.ttuple):
        return pa.struct(
            [pa.field(f'f{i}', hail_to_arrow_type(field_type)) for i, field_type in enumerate(dtype.types)]
        )
    raise TypeError(f'No Arrow equivalent configured for Hail type {dtype}')


def struct_to_arrow_schema(dtype: hl.tstruct) -> pa.Schema:
    return pa.schema([pa.field(name, hail_to_arrow_type(field_type)) for name, field_type in dtype.items()])


def to_arrow_value(value: Any, dtype: hl.HailType) -> Any:  # noqa: PLR0911
    """Convert a collected Hail value into plain python objects matching the type from hail_to_arrow_type."""
    if value is None:
        return None
    if isinstance(dtype, hl.tstruct):
        return {name: to_arrow_value(value[name], field_type) for name, field_type in dtype.items()}
    if isinstance(dtype, hl.ttuple):
        return {
            f'f{i}': to_arrow_value(each_value, field_type)
            for i, (each_value, field_type) in enumerate(zip(value, dtype.types, strict=True))
        }
    if isinstance(dtype, hl.tarray | hl.tset):
        return [to_arrow_value(each_value, dtype.element_type) for each_value in value]
    if isinstance(dtype, hl.tdict):
        return [(to_arrow_value(k, dtype.key_type), to_arrow_value(v, dtype.value_type)) for k, v in value.items()]
    if dtype == hl.tcall or isinstance(dtype, hl.tlocus | hl.tinterval):
        return str(value)
    return value


def structs_to_arrow_table(values: list, dtype: hl.tstruct, schema: pa.Schema) -> pa.Table:
    # missing structs (e.g. filtered entries) become a row of nulls
    return pa.Table.from_pylist([to_arrow_value(value, dtype) or {} for value in values], schema=schema)


def mt_to_parquet(input_mt: str, output: str) -> None:
    """
    Read the MT in, and write out to a single long-format Parquet file

    Args:
        input_mt (str): path of the MT to read in
        output (str): path of the Parquet file to generate
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['mt_to_parquet', 'worker_memory']),
        driver_memory=config.config_retrieve(['mt_to_parquet', 'driver_memory']),
        driver_cores=config.config_retrieve(['mt_to_parquet', 'driver_cores']),
    )
    partitions_per_chunk = config.config_retrieve(['mt_to_parquet', 'partitions_per_chunk'])

    mt = hl.read_matrix_table(input_mt)
    mt = mt.drop(*[field for field in REDUNDANT_ROW_FIELDS if field in mt.row])
    n_partitions = mt.n_partitions()

    # one row per variant, with all entries in an array, and the columns as an ordered array in the globals
    ht = mt.localize_entries(entries_array_field_name=ENTRIES_FIELD, columns_array_field_name=COLS_FIELD)
    col_type = ht[COLS_FIELD].dtype.element_type
    entry_type = ht[ENTRIES_FIELD].dtype.element_type
    col_values = hl.eval(ht.index_globals()[COLS_FIELD])
    n_cols = len(col_values)

    # split the locus into top-level columns, so Parquet row group statistics can be used to find genomic regions
    ht = ht.key_by()
    ht = ht.select_globals()
    ht = ht.select(
        chrom=ht.locus.contig,
        pos=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
        **{field: ht[field] for field in ht.row if field not in ('locus', 'alleles', ENTRIES_FIELD)},
        **{ENTRIES_FIELD: ht[ENTRIES_FIELD]},
    )
    row_type = hl.tstruct(**{field: field_type for field, field_type in ht.row.dtype.items() if field != ENTRIES_FIELD})

    row_schema = struct_to_arrow_schema(row_type)
    col_schema = struct_to_arrow_schema(col_type)
    entry_schema = struct_to_arrow_schema(entry_type)
    schema = pa.schema([*row_schema, *col_schema, *entry_schema])

    col_table = structs_to_arrow_table(col_values, col_type, col_schema)

    logger.info(f'Writing {n_partitions} partitions x {n_cols} columns to {output}, {partitions_per_chunk} at a time')

    total_variants = 0
    writer = pq.ParquetWriter(output, schema)
    try:
        for start in range(0, n_partitions, partitions_per_chunk):
            end = min(start + partitions_per_chunk, n_partitions)
            rows = ht._filter_partitions(range(start, end)).collect()  # noqa: SLF001
            if not rows:
                continue

            row_table = structs_to_arrow_table(rows, row_type, row_schema)
            entry_table = structs_to_arrow_table(
                [entry for row in rows for entry in row[ENTRIES_FIELD]],
                entry_type,
                entry_schema,
            )

            # expand to one row per variant x sample, matching the variant-major order of the flattened entries
            row_indices = pa.array(np.repeat(np.arange(len(rows)), n_cols))
            col_indices = pa.array(np.tile(np.arange(n_cols), len(rows)))
            chunk = pa.Table.from_arrays(
                [*row_table.take(row_indices).columns, *col_table.take(col_indices).columns, *entry_table.columns],
                schema=schema,
            )
            writer.write_table(chunk)

            total_variants += len(rows)
            logger.info(f'Written partitions {start}-{end - 1} of {n_partitions}, {total_variants} variants so far')
    finally:
        writer.close()

    logger.info(f'Written {total_variants} variants x {n_cols} columns = {total_variants * n_cols} rows to {output}')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser(description='Write a long-format Parquet file from a MT')
    parser.add_argument('--input', required=True, help='Path to the MT to read in')
    parser.add_argument('--output', required=True, help='Path to write the Parquet file to')
    args = parser.parse_args()

    mt_to_parquet(input_mt=args.input, output=args.output)


if __name__ == '__main__':
    cli_main()
