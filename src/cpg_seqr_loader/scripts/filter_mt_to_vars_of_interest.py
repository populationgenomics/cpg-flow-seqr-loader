"""
Takes a path to a MatrixTable and a name prefix
Writes data into the test bucket for the dataset

Optionally sample IDs and a locus can be provided to reduce the output
If sample IDs are specified, output subset contains only those
--format arg sets the output type (mt|vcf|both), default is MT

new behaviour: `--genes` as a CLI argument, a list of whitespace-delimited ENSG IDs. Filters MT to rows/variants which
    have at least one of the query genes present in the consequence annotations
"""

from argparse import ArgumentParser

from cpg_utils import config
from cpg_utils.hail_batch import init_batch

import hail as hl

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


def filter_to_gene_ids(
    mt: hl.MatrixTable,
    gene_ids: set[str],
) -> hl.MatrixTable:
    # turn the python set of Strings into a Hail Set Expression
    hl_gene_set = hl.set(gene_ids)

    # return rows where at least one of the query gene IDs is in the row annotation
    return mt.filter_rows(hl.len(hl_gene_set.intersection(mt.geneIds)) > 0)


def main(
    mt_path: str,
    output_path: str,
):
    """
    Parameters
    ----------
    mt_path : path to input MatrixTable
    output_path : output location of a TSV file
    """

    mt = hl.read_matrix_table(mt_path)

    if gene_list := config.config_retrieve(['alphagenome_params', 'genes'], []):
        mt = filter_to_gene_ids(mt, gene_list)

        # Filter for rare variants within the callset
        filtered_mt = mt.filter_rows(
            mt.info.AF[0] < config.config_retrieve(['alphagenome_params', 'callset_af'], CALLSET_AF)
        )

        # Filter for rare variants in gnomAD
        filtered_mt = filtered_mt.filter_rows(
            config.config_retrieve(['alphagenome_params', 'pop_af'], POP_AF) > filtered_mt.gnomad_genomes.AF
        )

        filtered_mt = hl.variant_qc(filtered_mt)

        filtered_mt_gq = filtered_mt.filter_rows(
            filtered_mt.variant_qc.gq_stats.mean > config.config_retrieve(['alphagenome_params', 'gq_threshold'], GQ_THRESHOLD)
        )
        exploded = filtered_mt_gq.explode_rows(filtered_mt_gq.vep.transcript_consequences)

        interesting_consequences = hl.set(config.config_retrieve(['alphagenome_params', 'consequences'], CSQ))
        results = exploded.filter_rows(
            exploded.vep.transcript_consequences.consequence_terms.any(
                lambda term: interesting_consequences.contains(term)
            )
        )

        # only need the rows
        results = results.rows()

        selected = results.annotate(
            CHROM=results.locus.contig,
            POS=results.locus.position,
            REF=results.alleles[0],
            ALT=results.alleles[1],
            FILTER_CRITERIA=hl.delimit(
                interesting_consequences.intersection(hl.set(results.vep.transcript_consequences.consequence_terms))
            ),
        )
        deduped = selected.key_by('CHROM', 'POS', 'REF', 'ALT', 'FILTER_CRITERIA').select().distinct()

        # Re-key by the selected fields to remove duplicates
        deduped.export(output_path)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the input MatrixTable', required=True)
    parser.add_argument(
        '--out',
        help='Full prefix for MT/VCF name\n("output" will become output.vcf.bgz or output.mt)',
        required=True,
    )

    args = parser.parse_args()

    init_batch()

    main(
        mt_path=args.i,
        output_path=args.out,
    )
