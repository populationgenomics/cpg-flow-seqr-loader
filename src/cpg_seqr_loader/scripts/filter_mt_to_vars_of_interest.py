#!/usr/bin/env python3
# ruff: noqa: PLR2004

"""
Takes a path to a MatrixTable and a name prefix
Writes data into the test bucket for the dataset

Optionally sample IDs and a locus can be provided to reduce the output
If sample IDs are specified, output subset contains only those
--format arg sets the output type (mt|vcf|both), default is MT

new behaviour: `--genes` as a CLI argument, a list of whitespace-delimited ENSG IDs. Filters MT to rows/variants which
    have at least one of the query genes present in the consequence annotations
"""

import logging
import sys
import argparse
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch
from cpg_utils import hail_batch, config


"""
The code reads a Hail MatrixTable, filters for rare variants (AF < 0.01), and attempts
to select variants with specific functional consequences (3' UTR, 5' UTR, splice region).
It exports deduplicated variant information to a TSV file.

Possible improvements:

Use more granular allele frequency thresholds (e.g., population-specific AF)
Select variants in specific genes or gene sets
Filter by additional functional consequences (e.g., missense, LoF)
Use external annotations (e.g., ClinVar, gnomAD)
Filter by genomic intervals (e.g., exons, promoters)
Apply sample-level filters (e.g., phenotype, ancestry)
"""




# Default values
POP_AF = 0.01
CALLSET_AF = 0.01
GQ_THRESHOLD = 35
CSQ = ['3_prime_UTR_variant', '5_prime_UTR_variant', 'splice_region_variant']
GENES = []


def filter_to_gene_ids(
    mt: hl.MatrixTable,
    gene_ids: set[str],
) -> hl.MatrixTable:
    # turn the python set of Strings into a Hail Set Expression
    hl_gene_set = hl.set(gene_ids)

    # return rows where at least one of the query gene IDs is in the row annotation
    return mt.filter_rows(hl.len(hl_gene_set.intersection(mt.geneIds)) > 0)


def get_thresholds_from_config() -> tuple[float, float, int]:
    """
    Quick method to allow fetching thresholds from the config file, instead of hard coding.

    Returns:
        A single tuple, containing Callset AF, Population AF, and GQ thresholds.
    """

    if thresholds := config.config_retrieve('alphagenome_params', None):
        return (
            thresholds.get('callset_af', CALLSET_AF),
            thresholds.get('pop_af', POP_AF),
            thresholds.get('gq_threshold', GQ_THRESHOLD),
            thresholds.get('genes', GENES)
        )
    else:
        return CALLSET_AF, POP_AF, GQ_THRESHOLD, GENES


def main(
    mt_path: str,
    output_path: str,
):
    """

    Parameters
    ----------
    mt_path : path to input MatrixTable
    prefix : prefix for file naming
    out_format : whether to write as a MT, VCF, or Both
    genes : optional, set of Ensembl ENSG gene IDs
    """

    mt = hl.read_matrix_table(mt_path)
    callset_af, pop_af, gq_threshold, genes = get_thresholds_from_config()

    if genes:
        mt = filter_to_gene_ids(mt, genes)

        # Filter for rare variants within the callset
        filtered_mt = mt.filter_rows(mt.info.AF[0] < callset_af)

        # Filter for rare variants in gnomAD
        filtered_mt = filtered_mt.filter_rows(filtered_mt.gnomad_genomes.AF < pop_af)

        filtered_mt = hl.variant_qc(filtered_mt)

        filtered_mt_GQ = filtered_mt.filter_rows(filtered_mt.variant_qc.gq_stats.mean > gq_threshold)
        exploded = filtered_mt_GQ.explode_rows(filtered_mt_GQ.vep.transcript_consequences)

        interesting_consequences = hl.set(config.config_retrieve('alphagenome_consequences', CSQ))
        results = exploded.filter_rows(
            exploded.vep.transcript_consequences.consequence_terms.any(
                lambda term: interesting_consequences.contains(term))
        )

        # only need the rows
        results = results.rows()

        # results = exploded.filter_rows(
        #     exploded.vep.transcript_consequences.consequence_terms.contains(
        #         '3_prime_UTR_variant' or '5_prime_UTR_variant' or 'splice_region_variant'
        #     )
        # )

        selected = results.annotate(
            CHROM=results.locus.contig,
            POS=results.locus.position,
            REF=results.alleles[0],
            ALT=results.alleles[1],
            ##    FILTER_CRITERIA=hl.delimit(
            ##        interesting_consequences.intersection(
            ##            hl.set(rows.vep.transcript_consequences.consequence_terms)
            ##        ))
        )
        deduped = selected.key_by('CHROM', 'POS', 'REF', 'ALT').select().distinct()

        # Re-key by the selected fields to remove duplicates
        deduped.export(output_path)



if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the input MatrixTable', required=True)
    parser.add_argument(
        '--out',
        help='Full prefix for MT/VCF name\n'
        '("output" will become output.vcf.bgz or output.mt)',
        required=True,
    )

    args = parser.parse_args()

    init_batch()

    main(
        mt_path=args.i,
        output_path=args.out,
    )
