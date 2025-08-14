import io
from argparse import ArgumentParser
from csv import DictReader

import matplotlib.pyplot as plt
import pandas as pd
from cloudpathlib.anypath import to_anypath
from cpg_utils import config

from alphagenome.data import genome,gene_annotation
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

ARBITRARY_THRESHOLD = config.config_retrieve(['alphagenome_params', 'sig_threshold'], default=0.98)
gtf = config.config_retrieve(
    ['alphagenome_params', 'gtf'],
    default='gs://cpg-common-main/references/alphagenome/gencode.v46.annotation.gtf.gz.feather',
)

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


def align_reference_for_indel(variant, interval, vout, length_alter: int):
    """Shift REF track in-place to align with ALT for indels (original logic)."""
    if length_alter > 0:  # deletion
        vout.reference.splice_sites.values[
            (variant.position - interval.start) : (interval.end - interval.start - length_alter)
        ] = vout.reference.splice_sites.values[
            (variant.position - interval.start + length_alter) : (interval.end - interval.start)
        ]
        vout.reference.splice_sites.values[
            (interval.end - interval.start - length_alter) : (interval.end - interval.start)
        ] = np.nan
    elif length_alter < 0:  # insertion
        vout.reference.splice_sites.values[
            (variant.position - interval.start - length_alter) : (interval.end - interval.start)
        ] = vout.reference.splice_sites.values[
            (variant.position - interval.start) : (interval.end - interval.start + length_alter)
        ]
        vout.reference.splice_sites.values[
            (variant.position - interval.start) : (variant.position - interval.start - length_alter)
        ] = np.nan
    # SNV => no shift


def load_transcript_extractor(gtf_path: str):
    fs = gcsfs.GCSFileSystem()
    with fs.open(gtf_path, 'rb') as f:
        gtf = pd.read_feather(f)
        gtf_t = gene_annotation.filter_protein_coding(gtf)
        gtf_t = gene_annotation.filter_to_longest_transcript(gtf_t)
        return transcript_utils.TranscriptExtractor(gtf_t)


def plot_variant_tracks(variant, interval, vout, transcript_extractor, plot_size: int, outpath: Path):
    transcripts = transcript_extractor.extract(interval)
    ref_output = vout.reference
    alt_output = vout.alternate
    ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}
    plot_elements = []
    variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']

    variant_scores = dna_model.score_variant(
        interval=interval, variant=variant, variant_scorers=[variant_scorer]
    )
    # Transcript annotation
    if transcripts:
        plot_elements.append(plot_components.TranscriptAnnotation(transcripts))

    # Sashimi plots
    ref_sashimi = ref_output.splice_junctions.filter_to_strand('-')
    alt_sashimi = alt_output.splice_junctions.filter_to_strand('-')
    if ref_sashimi and len(ref_sashimi) > 0:
        plot_elements.append(
            plot_components.Sashimi(ref_sashimi, ylabel_template='Reference {biosample_name} ({strand})\n{name}')
        )
    if alt_sashimi and len(alt_sashimi) > 0:
        plot_elements.append(
            plot_components.Sashimi(alt_sashimi, ylabel_template='Alternate {biosample_name} ({strand})\n{name}')
        )

    # Overlaid tracks — RNA Seq
    ref_rna = ref_output.rna_seq.filter_to_nonpositive_strand()
    alt_rna = alt_output.rna_seq.filter_to_nonpositive_strand()
    if ref_rna:
        plot_elements.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ref_rna, 'ALT': alt_rna},
                colors=ref_alt_colors,
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            )
        )

    # Overlaid tracks — Splice Sites
    ref_sites = ref_output.splice_sites.filter_to_nonpositive_strand()
    alt_sites = alt_output.splice_sites.filter_to_nonpositive_strand()
    if ref_sites:
        plot_elements.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ref_sites, 'ALT': alt_sites},
                colors=ref_alt_colors,
                ylabel_template='SPLICE SITES: {name} ({strand})',
            )
        )

    # Overlaid tracks — Splice Site Usage
    ref_usage = ref_output.splice_site_usage.filter_to_nonpositive_strand()
    alt_usage = alt_output.splice_site_usage.filter_to_nonpositive_strand()
    if ref_usage:
        plot_elements.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ref_usage, 'ALT': alt_usage},
                colors=ref_alt_colors,
                ylabel_template='SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}',
            )
        )

        # Sashimi plots
    ref_sashimi = ref_output.splice_junctions.filter_to_strand('+')
    alt_sashimi = alt_output.splice_junctions.filter_to_strand('+')
    if ref_sashimi and len(ref_sashimi) > 0:
        plot_elements.append(
            plot_components.Sashimi(ref_sashimi, ylabel_template='Reference {biosample_name} ({strand})\n{name}')
        )
    if alt_sashimi and len(alt_sashimi) > 0:
        plot_elements.append(
            plot_components.Sashimi(alt_sashimi, ylabel_template='Alternate {biosample_name} ({strand})\n{name}')
        )

    # Overlaid tracks — RNA Seq
    ref_rna = ref_output.rna_seq.filter_to_positive_strand()
    alt_rna = alt_output.rna_seq.filter_to_positive_strand()
    if ref_rna:
        plot_elements.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ref_rna, 'ALT': alt_rna},
                colors=ref_alt_colors,
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            )
        )

    # Overlaid tracks — Splice Sites
    ref_sites = ref_output.splice_sites.filter_to_positive_strand()
    alt_sites = alt_output.splice_sites.filter_to_positive_strand()
    if ref_sites:
        plot_elements.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ref_sites, 'ALT': alt_sites},
                colors=ref_alt_colors,
                ylabel_template='SPLICE SITES: {name} ({strand})',
            )
        )

    # Overlaid tracks — Splice Site Usage
    ref_usage = ref_output.splice_site_usage.filter_to_positive_strand()
    alt_usage = alt_output.splice_site_usage.filter_to_positive_strand()
    if ref_usage:
        plot_elements.append(
            plot_components.OverlaidTracks(
                tdata={'REF': ref_usage, 'ALT': alt_usage},
                colors=ref_alt_colors,
                ylabel_template='SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}',
            )
        )
    # Only plot if there's something to show
    if plot_elements:
        try:
            plot = plot_components.plot(
                plot_elements,
                interval=vout.reference.splice_sites.interval.resize(plot_size),
                annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
                title='Predicted REF vs. ALT effects of variant in Kidney tissue',
            )
            save_figure(outpath, plot)
            plt.close()
        except ValueError:  # raised if `y` is empty.
            print('No plot elements found. Skipping plot rendering.')


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
    transcript_extractor = load_transcript_extractor(gtf)
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
            requested_outputs={
                dna_client.OutputType.RNA_SEQ,
                dna_client.OutputType.SPLICE_SITES,
                dna_client.OutputType.SPLICE_SITE_USAGE,
                dna_client.OutputType.SPLICE_JUNCTIONS,
            }
        )

        # track count
        alt_vals = variant_prediction.alternate.splice_sites.values
        ref_vals = variant_prediction.reference.splice_sites.values
        n_tracks = alt_vals.shape[1]
        if n_tracks == 0:
            warnings.warn(f'No tracks available for {ontology}; skipping variant {var}.')
            continue

        # indel alignment
        length_alter = len(var.reference_bases) - len(var.alternate_bases)
        if length_alter != 0:
            align_reference_for_indel(var, interval, variant_prediction, length_alter)
            ref_vals = variant_prediction.reference.splice_sites.values  # updated

        plot_variant_tracks(var,interval,variant_prediction,transcript_extractor,plot_size=2**15, f'{var!s}.png')


    if significant_results is None:
        print('No significant results found.')
        return

    # Write the significant results to a TSV file via a StringIO buffer
    buffer = io.StringIO()
    significant_results.to_csv(buffer,sep='\t', index=False)
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
