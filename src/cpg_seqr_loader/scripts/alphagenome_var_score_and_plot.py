# ruff: noqa: PLR0915
import io
import warnings
from argparse import ArgumentParser
from csv import DictReader

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
from cloudpathlib.anypath import to_anypath
from cpg_utils import config
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas

"""
Notes:
This was originally written to run with a config file, and churn through variants of interest to produce
a pdf summary and a tsv file of significant results. It will be kept in the script directory in case we want to use
or adapt it in future.


max raw score is also associated with a position, i.e the most affected site in the 1MB window, would have been
helpful to track that position and see if multiple variants hit the same position
"""
# This is a configurable parameter that can be set in the config file.
# It determines the threshold for significance in variant scoring.
THRESHOLD = config.config_retrieve(['alphagenome_params', 'sig_threshold'], default=0.99)

# This is the path to the GTF file used for gene annotation.
gtf = config.config_retrieve(
    ['alphagenome_params', 'gtf'],
    default='gs://cpg-common-main/references/alphagenome/gencode.v46.annotation.gtf.gz.feather',
)

# This is a list of variant scorers that will be used to score the variants.
SCORER_CHOICES = [
    'RNA_SEQ',
    'SPLICE_JUNCTIONS',
    'SPLICE_SITES',
    'SPLICE_SITE_USAGE',
]


def pngs_to_pdf_streaming(directory, summary_text, variant_max_scores):
    """
    Save PNG files from cloud storage directly to a PDF in streaming mode.
    Images are loaded one at a time, written immediately to PDF page by page.
    """
    dir_path = to_anypath(directory)
    png_files = sorted([f for f in dir_path.iterdir() if f.suffix.lower() == '.png'])
    if not png_files:
        print('No PNG files found.')
        return

    pdf_path = to_anypath(f'{dir_path!s}/all_significant_vars.pdf')

    # Open PDF for writing in cloud storage
    with pdf_path.open('wb') as out_file:
        c = canvas.Canvas(out_file, pagesize=landscape(letter))
        width, height = landscape(letter)

        # --- Page 1: summary text ---
        c.setFont('Helvetica', 16)
        lines = summary_text.split('\n')
        y = height / 2 + (len(lines) // 2) * 20  # start a bit above center
        for line in lines:
            c.drawCentredString(width / 2, y, line)
            y -= 20  # move down for next line
        c.showPage()  # move to next page

        # --- Remaining pages: PNG images ---
        for png_file in png_files:
            with png_file.open('rb') as f:
                img = ImageReader(f)  # load only this image

            # Fit image to page while preserving aspect ratio
            img_width, img_height = img.getSize()
            scale = min(width / img_width, height / img_height) * 0.95  # leave margin
            draw_width = img_width * scale
            draw_height = img_height * scale
            x = (width - draw_width) / 2
            y = (height - draw_height) / 2

            c.drawImage(img, x, y, draw_width, draw_height)
            # Draw the filename as title above the image
            c.setFont('Helvetica-Bold', 14)
            c.drawCentredString(width / 2, height - 40, png_file.stem)

            # Draw max_raw_score just below the image
            variant_id = png_file.stem  # assuming the stem is a valid variant string like "chr1-123456-A-T"
            try:
                # Try to match the variant from the keys
                variant_obj = next(v for v in variant_max_scores if str(v) == variant_id)
                max_score = variant_max_scores[variant_obj]
                c.setFont('Helvetica', 12)
                c.drawCentredString(width / 2, y - 20, f'Max Raw Score: {max_score:.4f}')
            except StopIteration:
                c.setFont('Helvetica', 12)
                c.drawCentredString(width / 2, y - 20, 'Max Raw Score: N/A')
            c.showPage()  # finalize page

        c.save()


def save_figure(path_logdir, fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    # reset buffer position to the beginning
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
        vout.reference.splice_sites.values[  # noqa: PD011
            (variant.position - interval.start) : (interval.end - interval.start - length_alter)
        ] = vout.reference.splice_sites.values[  # noqa: PD011
            (variant.position - interval.start + length_alter) : (interval.end - interval.start)
        ]
        vout.reference.splice_sites.values[  # noqa: PD011
            (interval.end - interval.start - length_alter) : (interval.end - interval.start)
        ] = np.nan
    elif length_alter < 0:  # insertion
        vout.reference.splice_sites.values[  # noqa: PD011
            (variant.position - interval.start - length_alter) : (interval.end - interval.start)
        ] = vout.reference.splice_sites.values[  # noqa: PD011
            (variant.position - interval.start) : (interval.end - interval.start + length_alter)
        ]
        vout.reference.splice_sites.values[  # noqa: PD011
            (variant.position - interval.start) : (variant.position - interval.start - length_alter)
        ] = np.nan
    # SNV => no shift


def load_transcript_extractor(gtf_path: str):
    gtf = pd.read_feather(to_anypath(gtf_path))
    gtf_t = gene_annotation.filter_protein_coding(gtf)
    gtf_t = gene_annotation.filter_to_longest_transcript(gtf_t)
    return transcript_utils.TranscriptExtractor(gtf_t)


def plot_variant_tracks(variant, vout, transcript_extractor, outpath: str, significant_types: set[str]):
    """Plot the variant tracks for a given variant and save the figure.
    Args:
        variant (genome.Variant): The variant to plot.
        vout (dna_client.VariantOutputs): The outputs from the variant prediction.
        transcript_extractor (transcript_utils.TranscriptExtractor): Extractor for transcripts.
        outpath (str): Path to save the output figure.
        significant_types: set of significant types to plot, e.g. RNA_SEQ, SPLICE_JUNCTIONS, etc.
    """
    print(f'{variant!s} {significant_types}')
    plot_size = 2**10
    ref_output = vout.reference
    alt_output = vout.alternate
    ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}
    plot_elements = []
    plot_negative, plot_positive = False, False
    # Transcript annotation - tighten this window to limit the view to transcripts we'll actually plot
    if transcripts := transcript_extractor.extract(variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_100KB)):
        for transcript in transcripts:
            if transcript.is_negative_strand:
                plot_negative = True
            if transcript.is_positive_strand:
                plot_positive = True
        plot_elements.append(plot_components.TranscriptAnnotation(transcripts))

    if plot_negative:
        if 'SPLICE_JUNCTIONS' in significant_types:
            # Sashimi plots
            ref_sashimi = ref_output.splice_junctions.filter_to_strand('-')
            alt_sashimi = alt_output.splice_junctions.filter_to_strand('-')
            if ref_sashimi and len(ref_sashimi) > 0:
                plot_elements.append(
                    plot_components.Sashimi(
                        ref_sashimi, ylabel_template='Splice Junctions: Reference {biosample_name} ({strand})\n{name}'
                    )
                )
            if alt_sashimi and len(alt_sashimi) > 0:
                plot_elements.append(
                    plot_components.Sashimi(
                        alt_sashimi, ylabel_template='Splice Junctions: Alternate {biosample_name} ({strand})\n{name}'
                    )
                )

        if 'RNA_SEQ' in significant_types:
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

        if 'SPLICE_SITES' in significant_types:
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

        if 'SPLICE_SITE_USAGE' in significant_types:
            # Overlaid tracks — Splice Site Usage
            ref_usage = ref_output.splice_site_usage.filter_to_nonpositive_strand()
            alt_usage = alt_output.splice_site_usage.filter_to_nonpositive_strand()
            if ref_usage or alt_usage:
                plot_elements.append(
                    plot_components.OverlaidTracks(
                        tdata={'REF': ref_usage, 'ALT': alt_usage},
                        colors=ref_alt_colors,
                        ylabel_template='SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}',
                    )
                )

    if plot_positive:
        # if 'SPLICE_JUNCTIONS' in significant_types:
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

        if 'RNA_SEQ' in significant_types:
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

        if 'SPLICE_SITES' in significant_types:
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

        if 'SPLICE_SITE_USAGE' in significant_types:
            # Overlaid tracks — Splice Site Usage
            ref_usage = ref_output.splice_site_usage.filter_to_positive_strand()
            alt_usage = alt_output.splice_site_usage.filter_to_positive_strand()
            if ref_usage or alt_usage:
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
                fig_width=30,
                interval=vout.reference.splice_sites.interval.resize(plot_size),
                annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
                title=f'Predicted REF vs. ALT effects of variant in Kidney tissue\n'
                f'Significant types: {", ".join(sorted(significant_types))}',
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
    print(f'Loaded {len(variants)} variants from {input_variants!s}.')
    print(f'Using ontology terms: {ontology!s} with threshold {THRESHOLD}.')
    # Initialize counters for significant results
    # and a dictionary to track the number of significant variants per position
    sig_results_counter = 0
    sig_var_counter_dict: dict[tuple[str, int], int] = {}

    # Add a dictionary to store max_raw_score for each variant
    variant_max_scores = {}

    for var in variants:
        interval = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        key = (var.chromosome, var.position)
        variant_scores = model.score_variant(
            interval=interval,
            variant=var,
            variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[tool] for tool in SCORER_CHOICES],
        )

        tidied_scores = variant_scorers.tidy_scores(
            [variant_scores],
            match_gene_strand=True,
        )

        filtered_scores = tidied_scores[tidied_scores['quantile_score'] >= THRESHOLD]
        filtered_scores = filtered_scores[filtered_scores['ontology_curie'].isin(ontology)]
        max_raw_score = max(tidied_scores['raw_score'].values, default=0)

        # Store max_raw_score for later use
        variant_max_scores[var] = max_raw_score

        # If there are no scores above the threshold, skip to the next variant
        if filtered_scores.empty:
            print(f'No scores above threshold for variant {var}.')
            continue

        # either extend the significant results DataFrame or initialize it
        if significant_results is None:
            significant_results = filtered_scores
        else:
            # track the number of significant results for this variant
            sig_var_counter_dict[key] = sig_var_counter_dict.get(key, 0) + 1
            sig_results_counter += 1
            significant_results = pd.concat([significant_results, filtered_scores])

        # track the types of significant scores
        significant_types = set(filtered_scores.output_type.to_numpy().flatten())

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
            },
        )

        # track count
        n_tracks = variant_prediction.alternate.splice_sites.num_tracks
        if n_tracks == 0:
            warnings.warn(f'No tracks available for {ontology}; skipping variant {var}.', stacklevel=2)
            continue

        # indel alignment
        # If the variant is an indel, align the reference splice sites to the alternate
        length_alter = len(var.reference_bases) - len(var.alternate_bases)
        if length_alter != 0:
            align_reference_for_indel(var, interval, variant_prediction, length_alter)
        # Plot and save the variant tracks
        plot_variant_tracks(
            var,
            variant_prediction,
            transcript_extractor,
            f'{output_root!s}/{var!s}.png',
            significant_types=significant_types,
        )

    if significant_results is None:
        print('No significant results found.')
        return

    # Write the significant results to a TSV file via a StringIO buffer
    buffer = io.StringIO()
    significant_results.to_csv(buffer, sep='\t', index=False)
    buffer.seek(0)
    with to_anypath(f'{output_root}.tsv').open('w') as handle:
        handle.write(buffer.read())
    buffer.close()
    # Build top five summary string
    top_five_str = ''
    if sig_var_counter_dict:
        top_five = sorted(sig_var_counter_dict.items(), key=lambda x: x[1], reverse=True)[:5]
        for variant_key, count in top_five:
            chrom, pos = variant_key
            top_five_str += f'\nVariant {chrom}:{pos}: (count: {count})'

    summary_text = (
        f'Summary of Significant Variants: \n'
        f'{sig_results_counter} out of {len(variants)} variants\n '
        f'had significant scores above the threshold {THRESHOLD}.\n'
        f'Top 5 positions with significant variants:{top_five_str}'
    )
    pngs_to_pdf_streaming(output_root, summary_text, variant_max_scores)
    print(
        f'{sig_results_counter} out of {len(variants)} variants had significant scores above the threshold {THRESHOLD}.'
    )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the input MatrixTable', required=True)
    parser.add_argument('--out', help='prefix to write output files to)', required=True)
    parser.add_argument('--ontology', help='Ontology term to use for scoring', default=['UBERON:0002113'], nargs='+')
    parser.add_argument('--api_key', help='API key to authenticate', required=True)

    args = parser.parse_args()
    main(input_variants=args.i, output_root=args.out, ontology=args.ontology, api_key=args.api_key)
