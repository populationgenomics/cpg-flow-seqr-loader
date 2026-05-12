#!/usr/bin/env python3

"""
Script to read in the gnomAD 4.1 data table, use it to obtain evenly distributed (variant-balanced) intervals

n.b. requires an additional Mito interval... (does it? we should produce with/without, we have a mito calling workflow)

requires some additional re-splitting of the resulting intervals e.g. to have clean breaks at contig ends

target interval counts - 50, 100, 500, 1000, 2000

some inspo here
- https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/scripts/generate_new_intervals.py
- https://github.com/populationgenomics/production-pipelines/pull/883
- meres: gs://cpg-common-main/references/hg38/v0/hg38.telomeresAndMergedCentromeres.interval_list
"""

import logging
from argparse import ArgumentParser

from cpg_utils import hail_batch, to_path

import hail as hl


def get_naive_intervals(ht: hl.Table, intervals: int) -> list[tuple[str, int, int]]:
    """
    Get naive new intervals from a Table

    Args:
        ht (hl.Table): Input Table
        intervals (int): Number of intervals to generate

    Returns:
        list[tuple[str, int, int]]: List of intervals, in non-hail form
    """
    naive_interval_structs = ht._calculate_new_partitions(intervals)  # noqa: SLF001
    logging.info(f'Generated {len(naive_interval_structs)} intervals')

    # store both start and end positions separately
    new_intervals = [
        (
            interval.start.locus.contig,
            interval.start.locus.position,
            interval.end.locus.contig,
            interval.end.locus.position,
        )
        for interval in naive_interval_structs
    ]

    # check for cross-contig intervals, and split
    final_intervals: list[tuple[str, int, int]] = []
    parsed_chroms: set[str] = set()
    for index, (contig1, start, contig2, end) in enumerate(new_intervals, start=1):
        this_start = start
        # shift start position to 1, regardless of the variant positions used when generating the intervals
        if contig1 not in parsed_chroms:
            parsed_chroms.add(contig1)
            this_start = 1

        if contig1 != contig2:
            parsed_chroms.add(contig2)
            # get the end of the first
            first_end = hl.get_reference('GRCh38').lengths[contig1]
            final_intervals.append((contig1, this_start, first_end))
            final_intervals.append((contig2, 1, end))

        elif index == len(new_intervals) or new_intervals[index][0] != contig1:
            # if this is the final interval on this chromosome, shift the end position to the end of the chromosome
            # true if this is the latest interval on this chromosome
            final_intervals.append((contig1, this_start, hl.get_reference('GRCh38').lengths[contig1]))

        else:
            final_intervals.append((contig1, this_start, end))

    logging.info(f'Naive intervals post splitting: {len(final_intervals)}')

    return final_intervals


def overlaps(a: tuple[int, int], b: tuple[int, int]) -> int:
    """
    Return the amount of overlap, in bp
    between a and b.
    If >0, the number of bp of overlap
    If 0,  they are book-ended.
    If <0, the distance in bp between them
    """

    if (overlap := min(a[1], b[1]) - max(a[0], b[0])) > 0:
        return overlap
    return 0


def split_interval(interval: tuple[str, int, int], max_interval_size: int) -> list[tuple[str, int, int]]:
    """
    Split an interval into approximately evenly sized smaller intervals of at most max_interval_size bp
    """
    chrom, start, end = interval
    interval_size = end - start

    num_splits = (interval_size + max_interval_size - 1) // max_interval_size
    split_size = interval_size // num_splits

    new_intervals = []
    for i in range(num_splits):
        new_start = start + i * split_size
        new_end = new_start + split_size if i < num_splits - 1 else end
        new_intervals.append((chrom, new_start, new_end))

    return new_intervals


def read_meres(meres_path: str) -> tuple[dict[str, tuple[int, int]], dict[str, list[tuple[int, int]]]]:
    telomeres: dict[str, list[tuple[int, int]]] = {}
    centromeres: dict[str, tuple[int, int]] = {}

    logging.info(f'Reading centromere and telomere regions from {meres_path}')

    with to_path(meres_path).open('r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(('@SQ', '@HD')):
                continue

            contig, start_str, end_str, strand, region_type = line.strip().split('\t')
            if 'centromere' in region_type:
                centromeres[contig] = (int(start_str), int(end_str))
            else:
                telomeres.setdefault(contig, []).append((int(start_str), int(end_str)))

    return centromeres, telomeres


def polish_intervals(
    naive_intervals: list[tuple[str, int, int]],
    centromeres: dict[str, tuple[int, int]],
    telomeres: dict[str, list[tuple[int, int]]],
) -> list[tuple[str, int, int]]:
    """
    Polish intervals by splitting/removing centromere and telomere regions
    and setting a max interval length
    Args:
        naive_intervals (list): naive intervals, each limited to a single contig
        centromeres (dict): contig indexed tuple of centromere locations
        telomeres (dict): contig indexed tuple of telomere locations

    Returns:
        a polished version of the intervals file, where no centromere or telomere regions are present
    """

    new_intervals: list[tuple[str, int, int]] = []

    logging.info('Polishing intervals')
    for chrom, start, end in naive_intervals:
        this_start, this_end = start, end
        # does this overlap with a telomere?
        for telo in telomeres.get(chrom, []):
            # there's an overlap
            if overlaps(telo, (start, end)):
                # this is a 'start' telomere, shift the start coord
                if telo[0] == 1:
                    this_start = telo[1]
                else:
                    this_end = telo[0]
        # does it overlap with a centromere? If so split around it
        if (centro_region := centromeres.get(chrom)) and overlaps(centro_region, (start, end)):
            # check we have some interval left
            if centro_region[0] > start:
                new_intervals.append((chrom, this_start, centro_region[0]))
            if end > centro_region[1]:
                new_intervals.append((chrom, centro_region[1], this_end))
            continue

        new_intervals.append((chrom, this_start, this_end))

    logging.info(f'Final intervals: {len(new_intervals)}')

    return new_intervals


def main(input_path: str, count: list[int], meres: str, output: str):
    """Run all the things!"""

    hail_batch.init_batch()

    # read the input object using an appropriate method
    if input_path.endswith('.mt'):
        ht = hl.read_matrix_table(input_path).rows()
    elif input_path.endswith('.ht'):
        ht = hl.read_table(input_path)
    else:
        raise Exception(f'Unknown input format: {input_path}')

    # set up a path to do this once per interval count
    output_root = to_path(output)

    # parse the 'meres file
    centromeres, telomeres = read_meres(meres)

    for interval_count in count:
        intervals = get_naive_intervals(ht, interval_count)

        # do a little bit of interval grooming - hail will only divide the variants it sees, so there is a possibility
        # of as-yet-unseen variants occuring before the observed chr1 variants, or after the final observed variant

        # set the first base of the first interval on chr1 to 1, to maximise included regions
        first_chr, first_start, first_end = intervals[0]
        intervals[0] = (first_chr, 1, first_end)

        # remove a final interval on chrM if present
        if intervals[-1][0] == 'chrM':
            _chr_m_pop = intervals.pop(-1)

        # maximise the final non-chrM interval to the end of the contig
        last_chr, last_start, last_end = intervals[-1]
        intervals[-1] = (last_chr, last_start, hl.get_reference('GRCh38').lengths[last_chr])

        # and shove on a single max region for mitochondria, removing a previous mito interval if it was in the list
        intervals.append(('chrM', 1, 16569))

        with (output_root / f'{interval_count}_var_balanced_intervals.bed').open('w') as f:
            for contig, start, end in intervals:
                f.write(f'{contig}\t{start}\t{end}\n')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    parser = ArgumentParser()
    parser.add_argument('--input_path', required=True, type=str, help='Path to a Hail MT or HT.')
    parser.add_argument(
        '--count',
        nargs='+',
        type=int,
        required=True,
        help='Interval counts to generate.',
    )
    parser.add_argument(
        '--output',
        required=True,
        type=str,
        help='Dir to write output files.',
    )
    parser.add_argument(
        '--meres',
        default='gs://cpg-common-main/references/hg38/v0/hg38.telomeresAndMergedCentromeres.interval_list',
        help='Path to centro/telomere regions file.',
    )
    args = parser.parse_args()
    main(input_path=args.input_path, count=args.count, meres=args.meres, output=args.output)
