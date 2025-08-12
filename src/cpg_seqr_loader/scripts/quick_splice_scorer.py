import io
from argparse import ArgumentParser
from csv import DictReader

import matplotlib.pyplot as plt
import pandas as pd
from cloudpathlib.anypath import to_anypath
from cpg_utils import config

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

ARBITRARY_THRESHOLD = config.config_retrieve(['alphagenome_params', 'sig_threshold'], default=0.98)


def save_figure(path_logdir, fig):
    buf = io.BytesIO()
    fig.set_size_inches(12, 8)
    fig.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    with to_anypath(path_logdir).open('wb') as handle:
        handle.write(buf.read())
    buf.close()
    plt.close()


def load_variants_table(path: str):
    """
    Load variants from a TSV file.
    Args:
        path (str): Path to the TSV file containing variants.
    Returns:
        list[genome.Variant]: List of genome.Variant objects.
    """
    variants = []
    with to_anypath(path).open() as tsv_reader:
        reader = DictReader(tsv_reader, delimiter='\t')

        for row in reader:
            variants.append(
                genome.Variant(
                    chromosome=row['CHROM'],
                    position=int(row['POS']),
                    reference_bases=row['REF'],
                    alternate_bases=row['ALT'],
                )
            )
    return variants


def main(input_variants: str, output_root: str, ontology: list[str], api_key: str):
    """
    Main function to read variants from a TSV file and score them.
    For variants passing a significance threshold, it plots the splice site predictions.

    Args:
        input_variants (str): Path to the input TSV file containing variants.
        output_root (str): Prefix for the output files.
        ontology (list[str]): List of ontology terms to use for scoring.
        api_key (str): API key for authentication.
    """

    # as generated from the filter_mt_to_vars_of_interest.py script
    variants = load_variants_table(input_variants)

    significant_results: pd.DataFrame | None = None

    model = dna_client.create(api_key)

    for var in variants:
        interval = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

        variant_scores = model.score_variant(
            interval=interval,
            variant=var,
            variant_scorers=[
                variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITES'],
            ],
        )

        tidied_scores = variant_scorers.tidy_scores(
            [variant_scores],
            match_gene_strand=True,
        )

        filtered_scores = tidied_scores[tidied_scores['quantile_score'] >= ARBITRARY_THRESHOLD]

        # If there are no scores above the threshold, skip to the next variant
        if filtered_scores.empty:
            print(f'No scores above threshold for variant {var}.')
            continue

        # either extend the significant results DataFrame or initialize it
        if significant_results is None:
            significant_results = filtered_scores
        else:
            significant_results = pd.concat([significant_results, filtered_scores])

        # if there's significance, go on to plot the variant
        variant_prediction = model.predict_variant(
            interval=interval,
            variant=var,
            ontology_terms=ontology,
            requested_outputs=[dna_client.OutputType.SPLICE_SITES],
        )

        plot = plot_components.plot(
            [
                plot_components.OverlaidTracks(
                    tdata={
                        'REF': variant_prediction.reference.splice_sites,
                        'ALT': variant_prediction.alternate.splice_sites,
                    },
                    colors={'REF': 'dimgrey', 'ALT': 'red'},
                ),
            ],
            interval=variant_prediction.reference.splice_sites.interval.resize(2**14),
            # Annotate the location of the variant as a vertical line.
            annotations=[plot_components.VariantAnnotation([var], alpha=0.5)],
        )

        # use the buffered graph plotter
        save_figure(f'{var!s}.png', plot)
        break

    if significant_results is None:
        print('No significant results found.')
        return

    # Write the significant results to a TSV file via a StringIO buffer
    buffer = io.StringIO()
    significant_results.to_string(buffer)
    buffer.seek(0)
    with to_anypath(f'{output_root}.tsv').open('w') as handle:
        handle.write(buffer.read())
    buffer.close()


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the input MatrixTable', required=True)
    parser.add_argument('--out', help='prefix to write output files to)', required=True)
    parser.add_argument('--ontology', help='Ontology term to use for scoring', default=['UBERON:0001134'], nargs='+')
    parser.add_argument('--api_key', help='API key to authenticate', required=True)

    args = parser.parse_args()
    main(input_variants=args.i, output_root=args.out, ontology=args.ontology, api_key=args.api_key)
