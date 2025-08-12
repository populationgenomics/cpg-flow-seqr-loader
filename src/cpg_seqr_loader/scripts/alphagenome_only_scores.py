
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm
import gcsfs as gcsfs
from cloudpathlib.anypath import to_anypath

from cpg_utils import config

from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
#!/usr/bin/env python3
"""
alpha_variant_scan.py — Configurable & Stable AlphaGenome Variant Expression Scanner

Requires: alphagenome, pandas, numpy, matplotlib

This script scans predicted expression effects of genetic variants across organs using the AlphaGenome model.
It loads variant data and compares reference and alternate allele and applies the AlphaGenome model to score variants.
The script supports configurable parameters for input files, organs, thresholds and output formats.
"""




# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    P = argparse.ArgumentParser(
        prog='alpha_variant_scan',
        description='Scan predicted expression effects of variants across organs using AlphaGenome.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    P.add_argument('--variants', required=True, help='TSV/CSV/VCF of variants.')
    P.add_argument(
        '--output-table',
            default='gs://cpg-kidgen-test/alphagenome_scan_results.csv',
        help='Output summary table path (csv/tsv/xlsx).',
    )
    P.add_argument(
        '--output-table-sum',
        default='gs://cpg-kidgen-test/alphagenome_scan_results_variant_organ_summary.csv',
        help='Output organ summary table path (csv/tsv/xlsx).',
    )
    P.add_argument('--output-dir', default='gs://cpg-kidgen-test/alphagenome_scan_plots', help='Dir for plot images.')

    return P


# ------------------------------------------------------------------
# I/O
# --------------------------------------------


# ------------------------------------------------------------------
# Model / Annotation
# ------------------------------------------------------------------



api_key = config.config_retrieve(['alphagenome_params', 'api_key'], default=None)
organ = config.config_retrieve(
    ['alphagenome_params', 'organs'],
    default=[
        'UBERON:0000992',
        'UBERON:0002371',
        'UBERON:0000948',
        'UBERON:0000955',
        'UBERON:0001134',
        'UBERON:0001264',
    ],
)


def get_dna_model(api_key: str | None = None):
    if api_key is None:
        try:
            api_key = colab_utils.get_api_key()
        except Exception as e:
            # If not Colab
            print(f'Error getting API key from Colab: {e}')
            api_key = ''  # Insert your API key here if not working
    return dna_client.create(api_key)





def write_table(df: pd.DataFrame, path: str):
    ext = Path(path).suffix.lower()
    if ext in {'.tsv', '.txt'}:
        df.to_csv(path, sep='\t', index=False)
    elif ext in {'.xlsx', '.xls'}:
        df.to_excel(path, index=False)
    else:
        # default csv
        df.to_csv(path, index=False)

def check_gcs_path(gcs_path: str) -> str:
    """
    Check if the provided GCS path starts with 'b/' and convert it to 'gs://' format.
    """
    if gcs_path.startswith("b/"):
        print('Detected GCS path starting with "b/". Converting to gs:// format.')
        # convert from b/.../o/... to gs://bucket/object
        parts = gcs_path.split("/")
        bucket = parts[1]
        obj = "/".join(parts[3:])
        return f"gs://{bucket}/{obj}"
    return gcs_path


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    gcs_path = args.variants
    gcs_path = check_gcs_path(gcs_path)
    fs = gcsfs.GCSFileSystem()
    with fs.open(gcs_path, "rb") as f:
        variants_df = pd.read_csv(f, sep="\t")

    #transcript_extractor = #load_transcript_extractor(gtf)
    dna_model = get_dna_model(api_key)
    organism = 'human'  # @param ["human", "mouse"] {type:"string"}
    organs = organ or [
        'UBERON:0000992',
        'UBERON:0002371',
        'UBERON:0000948',
        'UBERON:0000955',
        'UBERON:0001134',
        'UBERON:0001264',
    ]
    #out_dir = to_anypath(args.output_dir)
    sequence_length = '1MB'  # @param ["2KB", "16KB", "100KB", "500KB", "1MB"] { type:"string" }
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]

    # @markdown Specify which scorers to use to score your variants:
    score_rna_seq = True  # @param { type: "boolean"}
    score_cage = True  # @param { type: "boolean" }
    score_procap = True  # @param { type: "boolean" }
    score_atac = True  # @param { type: "boolean" }
    score_dnase = True  # @param { type: "boolean" }
    score_chip_histone = True  # @param { type: "boolean" }
    score_chip_tf = True  # @param { type: "boolean" }
    score_polyadenylation = True  # @param { type: "boolean" }
    score_splice_sites = True  # @param { type: "boolean" }
    score_splice_site_usage = True  # @param { type: "boolean" }
    score_splice_junctions = True  # @param { type: "boolean" }

    # @markdown Other settings:
    download_predictions = False  # @param { type: "boolean" }

    # Parse organism specification.
    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism = organism_map[organism]

    # Parse scorer specification.
    scorer_selections = {
        'rna_seq': score_rna_seq,
        'cage': score_cage,
        'procap': score_procap,
        'atac': score_atac,
        'dnase': score_dnase,
        'chip_histone': score_chip_histone,
        'chip_tf': score_chip_tf,
        'polyadenylation': score_polyadenylation,
        'splice_sites': score_splice_sites,
        'splice_site_usage': score_splice_site_usage,
        'splice_junctions': score_splice_junctions,
    }

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key]
        for key in all_scorers
        if scorer_selections.get(key.lower(), False)
    ]
    # Remove any scorers or output types that are not supported for the chosen organism.
    unsupported_scorers = [
        scorer
        for scorer in selected_scorers
        if (
                   organism.value
                   not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
           )
           | (
                   (scorer.requested_output == dna_client.OutputType.PROCAP)
                   & (organism == dna_client.Organism.MUS_MUSCULUS)
           )
    ]
    if len(unsupported_scorers) > 0:
        print(
            f'Excluding {unsupported_scorers} scorers as they are not supported for'
            f' {organism}.'
        )
        for unsupported_scorer in unsupported_scorers:
            selected_scorers.remove(unsupported_scorer)

    ## revise this print statement to include all parameters
    print(
        f"Input variants path: {args.variants} | "
        f"Output table path: {args.output_table} | "
        #f"Output summary table path: {args.output_table_sum} | "
        #f"Output directory for plots: {args.output_dir} | "
    )

    # Score variants in the VCF file.
    results = []
    for ontology in organs:
        number_rank = 0
        for i, vcf_row in enumerate(tqdm(variants_df.itertuples(index=False), total=len(variants_df))):
            number_rank += 1
            variant = genome.Variant(
                chromosome=str(vcf_row.CHROM),
                position=int(vcf_row.POS),
                reference_bases=vcf_row.REF,
                alternate_bases=vcf_row.ALT,
                name=f'{ontology}_{number_rank}_{vcf_row.CHROM}_{vcf_row.POS}_{vcf_row.REF}_{vcf_row.ALT}',
            )
            interval = variant.reference_interval.resize(sequence_length)

            variant_scores = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=selected_scorers,
                organism=organism,
            )
            results.append(variant_scores)

    df_scores = variant_scorers.tidy_scores(results)
    write_table(df_scores, args.output_table)

if __name__ == '__main__':
    raise SystemExit(main())


