import io
import warnings
from argparse import ArgumentParser
from csv import DictReader
from PIL import Image

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

'''
Script to read a TSV of variants, score them with AlphaGenome, and plot sashimi plots for significant ones.
This is going to be archived in favor of the more integrated approach using Bedgraphs in IGV.
This script generates separate sashimi plots for each strand if there are significant splice junction changes.
It then merges the PNGs vertically into a single image per variant, and combines all images into a single PDF.
Requires an API key for AlphaGenome, and a TSV file with columns: CHROM, POS, REF, ALT.

Usage:
python sashimi_sep_plots.py -i input_variants.tsv --out output_prefix --ontology UBERON:0001134 --api_key YOUR_API_KEY

To DO:
- Add error handling for file operations and API calls.
- Optimize image merging for large numbers of variants.
- Make sure gene window is appropriate for the variant type.
'''

# This is a configurable parameter that can be set in the config file.
# It determines the threshold for significance in variant scoring.
THRESHOLD = config.config_retrieve(['alphagenome_params', 'sig_threshold'], default=0.99)

# This is the path to the GTF file used for gene annotation.
gtf = config.config_retrieve(
    ['alphagenome_params', 'gtf'],
    default='/Users/johass/PycharmProjects/cpg-flow-seqr-loader/gencode.v46.annotation.gtf.gz.feather',
)

# This is a list of variant scorers that will be used to score the variants.
SCORER_CHOICES = [
    'RNA_SEQ',
    'SPLICE_JUNCTIONS',
    'SPLICE_SITES',
    'SPLICE_SITE_USAGE',
]


def pngs_to_pdf_streaming(directory, summary_text):
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
            c.showPage()  # finalize page

        c.save()


def save_figure(path_logdir, figs):
    """
    Save multiple figures as separate PNG files, merge them vertically, and save to output path.

    :param path_logdir: Output path for the merged image
    :param figs: List of matplotlib figure objects to merge
    """
    if not isinstance(figs, list):
        figs = [figs]  # Handle single figure case

    buffers = []
    images = []

    try:
        # Save each figure to a separate buffer
        for i, fig in enumerate(figs):
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            buf.seek(0)
            buffers.append(buf)

            # Load image from buffer using PIL
            images.append(Image.open(buf))
            plt.close(fig)

        if not images:
            return

        # Calculate dimensions for merged image
        max_width = max(img.width for img in images)
        total_height = sum(img.height for img in images)

        # Create new image for merging
        merged_image = Image.new('RGB', (max_width, total_height), 'white')

        # Paste images vertically
        y_offset = 0
        for img in images:
            # Center horizontally if image is smaller than max width
            x_offset = (max_width - img.width) // 2
            merged_image.paste(img, (x_offset, y_offset))
            y_offset += img.height

        # Save merged image
        with to_anypath(path_logdir).open('wb') as handle:
            output_buf = io.BytesIO()
            merged_image.save(output_buf, format='PNG', dpi=(300, 300))
            output_buf.seek(0)
            handle.write(output_buf.read())
            output_buf.close()

    finally:
        # Clean up buffers and images
        for buf in buffers:
            buf.close()
        for img in images:
            img.close()

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
    gtf_t = gene_annotation.filter_to_longest_transcript(gtf)
    return transcript_utils.TranscriptExtractor(gtf_t)


def plot_variant_tracks(variant, vout, transcript_extractor, outpath: str, significant_types: set[str]):
    ref_output = vout.reference
    alt_output = vout.alternate

    # Prepare plot elements for each strand
    negative_elements = []
    positive_elements = []
    transcripts = transcript_extractor.extract(variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_100KB))
    if transcripts:
        neg_transcripts = [t for t in transcripts if t.is_negative_strand]
        pos_transcripts = [t for t in transcripts if t.is_positive_strand]
        if neg_transcripts:
            negative_elements.append(plot_components.TranscriptAnnotation(neg_transcripts))
        if pos_transcripts:
            positive_elements.append(plot_components.TranscriptAnnotation(pos_transcripts))

    # Sashimi plots for negative strand
    ref_sashimi_neg = ref_output.splice_junctions.filter_to_strand('-')
    alt_sashimi_neg = alt_output.splice_junctions.filter_to_strand('-')
    if ref_sashimi_neg and len(ref_sashimi_neg) > 0:
        negative_elements.append(plot_components.Sashimi(ref_sashimi_neg, ylabel_template='Reference (-)'))
    if alt_sashimi_neg and len(alt_sashimi_neg) > 0:
        negative_elements.append(plot_components.Sashimi(alt_sashimi_neg, ylabel_template='Alternate (-)'))

    # Sashimi plots for positive strand
    ref_sashimi_pos = ref_output.splice_junctions.filter_to_strand('+')
    alt_sashimi_pos = alt_output.splice_junctions.filter_to_strand('+')
    if ref_sashimi_pos and len(ref_sashimi_pos) > 0:
        positive_elements.append(plot_components.Sashimi(ref_sashimi_pos, ylabel_template='Reference (+)'))
    if alt_sashimi_pos and len(alt_sashimi_pos) > 0:
        positive_elements.append(plot_components.Sashimi(alt_sashimi_pos, ylabel_template='Alternate (+)'))

    padding = 200

    # Determine which subplots to plot based on existence of alt_sashimi_pos/neg
    plot_elements = []
    plot_intervals = []
    strand_labels = []
    axes_count = 0

    if len(alt_sashimi_neg.junctions) > 0:
        plot_chrom_neg = alt_sashimi_neg.junctions[0].chromosome
        plot_start_neg = alt_sashimi_neg.junctions[0].start - padding
        plot_end_neg = alt_sashimi_neg.junctions[-1].end + padding
        plot_interval_neg = genome.Interval(plot_chrom_neg, plot_start_neg, plot_end_neg)
        plot_elements.append(negative_elements)
        plot_intervals.append(plot_interval_neg)
        strand_labels.append('Negative strand')
        axes_count += 1



    if len(alt_sashimi_pos.junctions) > 0:
        plot_chrom_pos = alt_sashimi_pos.junctions[0].chromosome
        plot_start_pos = alt_sashimi_pos.junctions[0].start - padding
        plot_end_pos = alt_sashimi_pos.junctions[-1].end + padding
        plot_interval_pos = genome.Interval(plot_chrom_pos, plot_start_pos, plot_end_pos)
        plot_elements.append(positive_elements)
        plot_intervals.append(plot_interval_pos)
        strand_labels.append('Positive strand')
        axes_count += 1

    if axes_count == 0:
        warnings.warn(f'No sashimi junctions to plot for variant {variant}.', stacklevel=2)
        return

    # Create figure
    plt.figure(figsize=(30, 8 * axes_count))
    # Plot for each available strand
    if axes_count == 2:
        plot1 = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(neg_transcripts),
                plot_components.Sashimi(ref_sashimi_neg, ylabel_template='Reference (-)'),
                plot_components.Sashimi(alt_sashimi_neg, ylabel_template='Alternate (-)')],
            interval=plot_interval_neg,
            annotations=[plot_components.VariantAnnotation([variant])])
        plot2 = plot_components.plot(
                [
                plot_components.TranscriptAnnotation(pos_transcripts),
                plot_components.Sashimi(ref_sashimi_pos, ylabel_template='Reference (+)'),
                plot_components.Sashimi(alt_sashimi_pos, ylabel_template='Alternate (+)')],
            interval=plot_interval_pos,
            annotations=[plot_components.VariantAnnotation([variant])])
        plots=[plot1, plot2]
    else:
        plot = plot_components.plot(
            plot_elements[0],
            interval=plot_intervals[0],
            annotations=[plot_components.VariantAnnotation([variant])])
        plots=[plot]
    save_figure(outpath, plots)
    plt.close()
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

    output_metadata = model.output_metadata(
        dna_client.Organism.HOMO_SAPIENS
    ).concatenate()

    print(f'Loaded {len(variants)} variants from {input_variants!s}.')
    print(f'Using ontology terms: {ontology!s} with threshold {THRESHOLD}.')
    # Initialize counters for significant results
    # and a dictionary to track the number of significant variants per position
    sig_results_counter = 0
    sig_var_counter_dict: dict[tuple[str, int], int] = {}
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
    pngs_to_pdf_streaming(output_root, summary_text)
    print(
        f'{sig_results_counter} out of {len(variants)} variants had significant scores above the threshold {THRESHOLD}.'
    )



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the input MatrixTable', default='/Users/johass/PycharmProjects/cpg-flow-seqr-loader/TEST_VARS.tsv')
    parser.add_argument('--out', help='prefix to write output files to)' ,default='/Users/johass/PycharmProjects/cpg-flow-seqr-loader/src/test_out_sashimi_blood')
    parser.add_argument('--ontology', help='Ontology term to use for scoring', default=['UBERON:0001134'], nargs='+')
    parser.add_argument('--api_key', help='API key to authenticate',  default='AIzaSyDYx7VMDPepU7qeJOm7i-AVm9GsrV-BbW8')

    args = parser.parse_args()
    main(input_variants=args.i, output_root=args.out, ontology=args.ontology, api_key=args.api_key)
