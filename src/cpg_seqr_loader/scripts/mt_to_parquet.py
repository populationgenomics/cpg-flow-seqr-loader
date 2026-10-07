"""
Read in a MT, and write it out as a single long-format Parquet file, one row per variant x sample

Query-on-Batch has no Spark, so we can't use to_spark().write.parquet(), and collecting the MT through the driver
doesn't scale. Instead this runs in two phases, neither of which holds more than a bounded amount of data in memory:

1. Hail localises the entries (one row per variant, with an array of all the entries), converts all fields to
   JSON/Arrow-compatible types, and the workers export the variant rows as block-gzipped newline-delimited JSON shards
2. This process streams the shards in order through pyarrow's JSON reader in fixed-size blocks, expands each block to
   one row per variant x sample, and writes those as row groups to a single Parquet file

The MT is sorted by locus, and contig/position are top-level columns, so the row group min/max statistics act as an
index of the genomic interval covered by each row group.
"""

import json
from argparse import ArgumentParser
from collections.abc import Iterator

import numpy as np
import pyarrow as pa
import pyarrow.fs as pafs
import pyarrow.json as pj
import pyarrow.parquet as pq
from loguru import logger

import hail as hl

from cpg_utils import config, hail_batch

# row annotations added in AnnotateDataset which summarise all samples per variant. These are redundant in a
# long-format table, and would be repeated for every sample on every row
REDUNDANT_ROW_FIELDS = ['genotypes', 'samples_no_call', 'samples_num_alt', 'samples_gq', 'samples_ab']

ENTRIES_FIELD = 'mt_entries'
COLS_FIELD = 'mt_cols'

SHARD_MANIFEST = 'shard-manifest.txt'

PRIMITIVE_TYPES = {
    hl.tstr: pa.string(),
    hl.tint32: pa.int32(),
    hl.tint64: pa.int64(),
    hl.tfloat32: pa.float32(),
    hl.tfloat64: pa.float64(),
    hl.tbool: pa.bool_(),
}


def to_arrow_compatible(expr: hl.Expression) -> hl.Expression:  # noqa: PLR0911
    """
    Recursively convert a Hail expression into types which have a direct JSON & Arrow representation

    - Calls, Loci, and Intervals become strings, e.g. '0/1', 'chr1:12345'
    - Sets become arrays, Tuples become structs with fields f0, f1, ...
    - Dicts become arrays of key/value structs (pyarrow's JSON reader doesn't support map types)
    - Non-finite floats become missing, as Hail writes them to JSON as strings
    """
    dtype = expr.dtype
    if isinstance(dtype, hl.tstruct):
        return hl.struct(**{name: to_arrow_compatible(expr[name]) for name in dtype})
    if isinstance(dtype, hl.ttuple):
        return hl.struct(**{f'f{i}': to_arrow_compatible(expr[i]) for i in range(len(dtype.types))})
    if isinstance(dtype, hl.tarray | hl.tset):
        return hl.array(expr).map(to_arrow_compatible)
    if isinstance(dtype, hl.tdict):
        return hl.array(expr).map(
            lambda kv: hl.struct(key=to_arrow_compatible(kv[0]), value=to_arrow_compatible(kv[1])),
        )
    if dtype == hl.tcall or isinstance(dtype, hl.tlocus | hl.tinterval):
        return hl.str(expr)
    if dtype in (hl.tfloat32, hl.tfloat64):
        return hl.or_missing(hl.is_finite(expr), expr)
    if dtype in PRIMITIVE_TYPES:
        return expr
    raise TypeError(f'No Arrow equivalent configured for Hail type {dtype}')


def hail_to_arrow_type(dtype: hl.HailType) -> pa.DataType:
    """Map a Hail type (already passed through to_arrow_compatible) to the equivalent Arrow type."""
    if dtype in PRIMITIVE_TYPES:
        return PRIMITIVE_TYPES[dtype]
    if isinstance(dtype, hl.tarray):
        return pa.list_(hail_to_arrow_type(dtype.element_type))
    if isinstance(dtype, hl.tstruct):
        return pa.struct([pa.field(name, hail_to_arrow_type(field_type)) for name, field_type in dtype.items()])
    raise TypeError(f'No Arrow equivalent configured for Hail type {dtype}')


def struct_to_arrow_schema(dtype: hl.tstruct) -> pa.Schema:
    return pa.schema([pa.field(name, hail_to_arrow_type(field_type)) for name, field_type in dtype.items()])


def export_json_shards(input_mt: str, json_dir: str) -> tuple[pa.Schema, pa.Table]:
    """
    Export the MT as variant-level NDJSON shards, each line holding the row fields and an array of all entries

    Returns:
        the Arrow schema of the JSON lines, and an Arrow table of the column fields in entry-array order
    """
    mt = hl.read_matrix_table(input_mt)
    mt = mt.drop(*[field for field in REDUNDANT_ROW_FIELDS if field in mt.row])

    # one row per variant, with all entries in an array, and the columns as an ordered array in the globals
    ht = mt.localize_entries(entries_array_field_name=ENTRIES_FIELD, columns_array_field_name=COLS_FIELD)

    # the column fields are small (one struct per sample), so evaluate them here
    cols = to_arrow_compatible(ht.index_globals()[COLS_FIELD])
    col_table = pa.Table.from_pylist(
        json.loads(hl.eval(hl.json(cols))),
        schema=struct_to_arrow_schema(cols.dtype.element_type),
    )

    # split the locus into top-level columns, so Parquet row group statistics can be used to find genomic regions
    ht = ht.key_by()
    ht = ht.select_globals()
    ht = ht.select(
        contig=ht.locus.contig,
        position=ht.locus.position,
        **{field: to_arrow_compatible(ht[field]) for field in ht.row if field not in ('locus', ENTRIES_FIELD)},
        **{ENTRIES_FIELD: to_arrow_compatible(ht[ENTRIES_FIELD])},
    )
    json_schema = struct_to_arrow_schema(ht.row.dtype)

    logger.info(f'Exporting {mt.n_partitions()} partitions x {col_table.num_rows} columns as JSON to {json_dir}')
    ht.select(json=hl.json(ht.row)).export(json_dir, header=False, parallel='separate_header')

    return json_schema, col_table


def get_ordered_shards(json_dir: str) -> tuple[pafs.FileSystem, list[str]]:
    """Use the manifest Hail writes alongside the shards to find them in partition order."""
    filesystem, base_path = pafs.FileSystem.from_uri(json_dir)
    with filesystem.open_input_stream(f'{base_path}/{SHARD_MANIFEST}') as handle:
        manifest = handle.read().decode().split()
    # the manifest may include the (empty) header file, only keep the data shards
    return filesystem, [f'{base_path}/{shard}' for shard in manifest if shard.startswith('part-')]


def expand_to_long_format(
    batch: pa.RecordBatch,
    col_table: pa.Table,
    schema: pa.Schema,
    max_rows: int,
) -> Iterator[pa.Table]:
    """
    Expand a batch of variant rows into one row per variant x sample, in slices of at most max_rows rows

    The row fields are repeated for each sample, so the expansion is done in slices to bound memory
    """
    n_cols = col_table.num_rows
    variants_per_slice = max(1, max_rows // n_cols)
    row_fields = [name for name in batch.schema.names if name != ENTRIES_FIELD]

    for offset in range(0, batch.num_rows, variants_per_slice):
        variants = batch.slice(offset, variants_per_slice)
        n_variants = variants.num_rows

        # flatten the per-variant entry arrays (variant-major), then split the entry structs into columns
        # StructArray.flatten() carries the null-ness of missing entries into each of the entry field columns
        entries = variants.column(ENTRIES_FIELD).flatten()
        if len(entries) != n_variants * n_cols:
            raise ValueError(f'Expected {n_variants * n_cols} entries, found {len(entries)}')

        row_indices = pa.array(np.repeat(np.arange(n_variants), n_cols))
        col_indices = pa.array(np.tile(np.arange(n_cols), n_variants))

        yield pa.Table.from_arrays(
            [
                *[variants.column(name).take(row_indices) for name in row_fields],
                *col_table.take(col_indices).columns,
                *entries.flatten(),
            ],
            schema=schema,
        )


def mt_to_parquet(input_mt: str, json_dir: str, output: str) -> None:
    """
    Read the MT in, and write out to a single long-format Parquet file

    Args:
        input_mt (str): path of the MT to read in
        json_dir (str): temporary directory to write the intermediate NDJSON shards to
        output (str): path of the Parquet file to generate
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['mt_to_parquet', 'worker_memory']),
        driver_memory=config.config_retrieve(['mt_to_parquet', 'driver_memory']),
        driver_cores=config.config_retrieve(['mt_to_parquet', 'driver_cores']),
    )
    block_size = config.config_retrieve(['mt_to_parquet', 'json_block_size_mb']) * 1024 * 1024
    rows_per_row_group = config.config_retrieve(['mt_to_parquet', 'rows_per_row_group'])

    json_schema, col_table = export_json_shards(input_mt, json_dir)

    # final schema - row fields, then column fields, then entry fields
    entry_type = json_schema.field(ENTRIES_FIELD).type.value_type
    schema = pa.schema(
        [
            *[field for field in json_schema if field.name != ENTRIES_FIELD],
            *col_table.schema,
            *[entry_type.field(i) for i in range(entry_type.num_fields)],
        ],
    )

    filesystem, shards = get_ordered_shards(json_dir)
    logger.info(f'Streaming {len(shards)} JSON shards into {output}')

    read_options = pj.ReadOptions(block_size=block_size)
    parse_options = pj.ParseOptions(explicit_schema=json_schema, unexpected_field_behavior='error')

    total_variants = 0
    with pq.ParquetWriter(output, schema) as writer:
        for index, shard in enumerate(shards, start=1):
            with filesystem.open_input_stream(shard, compression='gzip') as stream:
                for batch in pj.open_json(stream, read_options=read_options, parse_options=parse_options):
                    for table in expand_to_long_format(batch, col_table, schema, rows_per_row_group):
                        writer.write_table(table, row_group_size=rows_per_row_group)
                    total_variants += batch.num_rows
            logger.info(f'Written shard {index} of {len(shards)}, {total_variants} variants so far')

    n_cols = col_table.num_rows
    logger.info(f'Written {total_variants} variants x {n_cols} columns = {total_variants * n_cols} rows to {output}')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser(description='Write a long-format Parquet file from a MT')
    parser.add_argument('--input', required=True, help='Path to the MT to read in')
    parser.add_argument('--json_dir', required=True, help='Temporary directory for intermediate JSON shards')
    parser.add_argument('--output', required=True, help='Path to write the Parquet file to')
    args = parser.parse_args()

    mt_to_parquet(input_mt=args.input, json_dir=args.json_dir, output=args.output)


if __name__ == '__main__':
    cli_main()
