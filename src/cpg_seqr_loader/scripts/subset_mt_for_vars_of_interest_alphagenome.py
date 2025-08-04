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

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_mt', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()


hl.init()
mt=hl.read_matrix_table(args.input_mt)

##mt.aggregate_rows(hl.agg.count())         # Count of variants
##mt.aggregate_cols(hl.agg.collect(mt.s))   # List of sample IDs
##mt.aggregate_entries(hl.agg.mean(mt.DP))  # Mean depth across all entries


# Filter for rare variants
filtered_mt = mt.filter_rows(mt.info.AF[0] < 0.01)
filtered_mt = hl.sample_qc(filtered_mt)
filtered_mt = hl.variant_qc(filtered_mt)

filtered_mt_GQ=filtered_mt.filter_rows(filtered_mt.row.variant_qc.gq_stats.mean > 35)
exploded = filtered_mt_GQ.explode_rows(filtered_mt_GQ.vep.transcript_consequences)
#exploded = exploded.checkpoint('exploded.mt')

interesting_consequences = hl.set(['3_prime_UTR_variant', '5_prime_UTR_variant', 'splice_region_variant'])
results=exploded.filter_rows(exploded.vep.transcript_consequences.consequence_terms.contains("3_prime_UTR_variant" or '5_prime_UTR_variant' or'splice_region_variant'))


rows = results.rows().key_by()
selected = rows.select(
    CHROM = rows.locus.contig,
    POS = rows.locus.position,
    REF = rows.alleles[0],
    ALT = rows.alleles[1],
##    FILTER_CRITERIA=hl.delimit(
##        interesting_consequences.intersection(
##            hl.set(rows.vep.transcript_consequences.consequence_terms)
##        ))
)

# Re-key by the selected fields to remove duplicates
deduped = selected.key_by('CHROM', 'POS', 'REF', 'ALT').distinct()
deduped.export(args.out)