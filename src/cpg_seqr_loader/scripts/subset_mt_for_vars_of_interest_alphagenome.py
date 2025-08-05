import argparse
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

from cpg_utils import hail_batch, config


# Default values
POP_AF = 0.01
CALLSET_AF = 0.01
GQ_THRESHOLD = 35
CSQ = ['3_prime_UTR_variant', '5_prime_UTR_variant', 'splice_region_variant']


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    main(input_path=args.i, output_path=args.out)


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
        )
    else:
        return CALLSET_AF, POP_AF, GQ_THRESHOLD


def main(input_path: str, output_path: str):
    """Main function to run the script."""

    # start a batch with a service backend
    hail_batch.init_batch()

    mt = hl.read_matrix_table(input_path)

    callset_af, pop_af, gq_threshold = get_thresholds_from_config()

    # Filter for rare variants within the callset
    filtered_mt = mt.filter_rows(mt.info.AF[0] < callset_af)

    # Filter for rare variants in gnomAD
    filtered_mt = filtered_mt.filter_rows(filtered_mt.gnomad_genomes.AF < pop_af)

    # Todo remove, unused?
    # filtered_mt = hl.sample_qc(filtered_mt)

    filtered_mt = hl.variant_qc(filtered_mt)

    filtered_mt_GQ = filtered_mt.filter_rows(filtered_mt.variant_qc.gq_stats.mean > gq_threshold)
    exploded = filtered_mt_GQ.explode_rows(filtered_mt_GQ.vep.transcript_consequences)

    interesting_consequences = hl.set(config.config_retrieve('alphagenome_consequences', CSQ))
    results = exploded.filter_rows(
        exploded.vep.transcript_consequences.consequence_terms.any(lambda term: interesting_consequences.contains(term))
    )

    # only need the rows
    results = results.rows()

    # results = exploded.filter_rows(
    #     exploded.vep.transcript_consequences.consequence_terms.contains(
    #         '3_prime_UTR_variant' or '5_prime_UTR_variant' or 'splice_region_variant'
    #     )
    # )

    selected = (
        results.select(
            CHROM=results.locus.contig,
            POS=results.locus.position,
            REF=results.alleles[0],
            ALT=results.alleles[1],
            ##    FILTER_CRITERIA=hl.delimit(
            ##        interesting_consequences.intersection(
            ##            hl.set(rows.vep.transcript_consequences.consequence_terms)
            ##        ))
        )
        .key_by('CHROM', 'POS', 'REF', 'ALT')
        .distinct()
    )

    # Re-key by the selected fields to remove duplicates
    selected.export(output_path)


if __name__ == '__main__':
    cli_main()
