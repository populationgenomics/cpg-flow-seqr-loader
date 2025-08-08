#!/usr/bin/env python3
"""
alpha_variant_scan.py — Configurable & Stable AlphaGenome Variant Expression Scanner

need: alphagenome, pandas, numpy, matplotlib.

This script scans predicted expression effects of genetic variants across organs using the AlphaGenome model.
It loads variant data, applies a sliding window analysis to compare reference and alternate allele effects, identifies significant regions, generates summary tables, and optionally plots results for visualization.
The script supports configurable parameters for input files, organs, thresholds, window size, merging, and output formats.
"""

import argparse
import warnings
from pathlib import Path
from cloudpathlib.anypath import to_anypath

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

from cpg_utils import config


# Default values
DEFAULT_UBERON = [
    'UBERON:0000992',
    'UBERON:0002371',
    'UBERON:0000948',
    'UBERON:0000955',
    'UBERON:0001134',
    'UBERON:0001264',
]
# |ALT/REF−1| score threshold.
THRESHOLD = 0.5
# Min bp length of merged region.'
MIN_LENGTH = 1000
# Max bp gap to merge regions.
MERGE_DISTANCE = 300
# Sliding window size (bp).
WINDOW_SIZE = 100
# Bp each side of variant center to scan.
SCAN_SPAN = 50000


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    P = argparse.ArgumentParser(
        prog='alpha_variant_scan',
        description='Scan predicted expression effects of variants across organs using AlphaGenome.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    P.add_argument('--variants', required=True, help='TSV of variants.')
    P.add_argument('--output', required=True, help='Folder to write all results to.')
    P.add_argument('--api-key', required=True, help='AlphaGenome API key.')
    P.add_argument(
        '--gtf',
        default='https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather',
        help='GENCODE transcript annotation feather path/URL.',
    )
    return P


# ------------------------------------------------------------------
# Model / Annotation
# ------------------------------------------------------------------


def load_transcript_extractor(gtf_path: str):
    gtf = pd.read_feather(gtf_path)
    gtf_t = gene_annotation.filter_protein_coding(gtf)
    gtf_t = gene_annotation.filter_to_longest_transcript(gtf_t)
    return transcript_utils.TranscriptExtractor(gtf_t)


# ------------------------------------------------------------------
# Core computations
# ------------------------------------------------------------------


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


def compute_window_scores(alt_vals: np.ndarray, ref_vals: np.ndarray, start_idx: int, end_idx: int) -> np.ndarray:
    """Return (n_windows, n_tracks) of ALT/REF−1 window means."""

    window_size = config.config_retrieve(['alphagenome_params', 'window_size'], default=WINDOW_SIZE)

    n_bases, n_tracks = alt_vals.shape
    if start_idx < 0:
        start_idx = 0
    if end_idx > n_bases:
        end_idx = n_bases
    n_windows = end_idx - start_idx - window_size + 1
    if n_windows <= 0:
        return np.empty((0, n_tracks))
    # cumulative sums (prepend 0 row)
    alt_cs = np.cumsum(np.vstack([np.zeros((1, n_tracks)), alt_vals]), axis=0)
    ref_cs = np.cumsum(np.vstack([np.zeros((1, n_tracks)), ref_vals]), axis=0)
    # slice rolling windows
    a = (
        alt_cs[start_idx + window_size : start_idx + window_size + n_windows]
        - alt_cs[start_idx : start_idx + n_windows]
    )
    r = (
        ref_cs[start_idx + window_size : start_idx + window_size + n_windows]
        - ref_cs[start_idx : start_idx + n_windows]
    )
    a /= float(window_size)
    r /= float(window_size)
    return (a / (r + 1e-8)) - 1.0
    # Suggestions for alternative implementations that may be more representative of splicing events:
    # scores = np.empty((n_windows, n_tracks))
    # for i in range(n_windows):
    #    win_alt = alt_vals[start_idx + i:start_idx + i + window_size, :]
    #    win_ref = ref_vals[start_idx + i:start_idx + i + window_size, :]
    #    scores[i, :] = np.max(np.abs(win_alt - win_ref), axis=0)
    # return scores


def call_regions(scores: np.ndarray) -> list[tuple[int, int]]:
    """Identify high-score regions (inclusive window indices)."""

    # pull threshold, min_length, merge_distance from config
    threshold = config.config_retrieve(['alphagenome_params', 'threshold'], default=THRESHOLD)
    min_length = config.config_retrieve(['alphagenome_params', 'min_length'], default=MIN_LENGTH)
    merge_distance = config.config_retrieve(['alphagenome_params', 'merge_distance'], default=MERGE_DISTANCE)

    idx = np.where(np.abs(scores) > threshold)[0]
    if idx.size == 0:
        return []
    # compress consecutive indices
    regions = []
    start = prev = idx[0]
    for i in idx[1:]:
        if i == prev + 1:
            prev = i
        else:
            regions.append((start, prev))
            start = prev = i
    regions.append((start, prev))  # tail
    # merge
    merged = []
    cur_s, cur_e = regions[0]
    for s, e in regions[1:]:
        if s - cur_e <= merge_distance:
            cur_e = e
        else:
            if cur_e - cur_s + 1 >= min_length:
                merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    if cur_e - cur_s + 1 >= min_length:
        merged.append((cur_s, cur_e))
    return merged


# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------


def plot_variant_tracks(variant, interval, vout, transcript_extractor, plot_size: int, outpath: Path):
    transcripts = transcript_extractor.extract(interval)
    ref_output = vout.reference
    alt_output = vout.alternate
    ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}
    plot_elements = []

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
            plt.savefig(outpath, dpi=300, bbox_inches='tight')
            plt.close()
        except ValueError:  # raised if `y` is empty.
            print('No plot elements found. Skipping plot rendering.')
            pass


def plot_scores(scores: np.ndarray, regions: list[tuple[int, int]], threshold: float, outpath: Path, title: str):
    plt.figure(figsize=(8, 3))
    plt.plot(scores, label='score')
    plt.axhline(threshold, color='grey', lw=0.5, ls='--')
    plt.axhline(-threshold, color='grey', lw=0.5, ls='--')
    for s, e in regions:
        plt.axvspan(s, e, color='red', alpha=0.3)
    plt.xlabel('window index')
    plt.ylabel('ALT/REF−1')
    plt.title(title)
    plt.legend(loc='upper right', fontsize='x-small')
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# ------------------------------------------------------------------
# Table writer
# ------------------------------------------------------------------


def write_table(df: pd.DataFrame, path: str):
    ext = Path(path).suffix.lower()
    if ext in {'.tsv', '.txt'}:
        df.to_csv(path, sep='\t', index=False)
    elif ext in {'.xlsx', '.xls'}:
        df.to_excel(path, index=False)
    else:
        # default csv
        df.to_csv(path, index=False)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    output_as_path = to_anypath(args.output)

    # load inputs
    variants_df = df = pd.read_csv(to_anypath(args.variants), sep='\t')
    transcript_extractor = load_transcript_extractor(args.gtf)
    dna_model = dna_client.create(args.api_key)

    results_rows = []

    scan_span = config.config_retrieve(['alphagenome_params', 'scan_span'], default=SCAN_SPAN)
    window_size = config.config_retrieve(['alphagenome_params', 'window_size'], default=WINDOW_SIZE)

    for ontology in config.config_retrieve(['alphagenonme_params', 'organs'], DEFAULT_UBERON):
        # made a change here, using the DF index instead of a counter
        for number_rank, row in variants_df.iterrows():
            variant = genome.Variant(
                chromosome=row['CHROM'],
                position=row['POS'],
                reference_bases=row['REF'],
                alternate_bases=row['ALT'],
            )

            # 1Mb centered interval
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

            vout = dna_model.predict_variant(
                interval=interval,
                variant=variant,
                requested_outputs={
                    dna_client.OutputType.RNA_SEQ,
                    dna_client.OutputType.SPLICE_SITES,
                    dna_client.OutputType.SPLICE_SITE_USAGE,
                    dna_client.OutputType.SPLICE_JUNCTIONS,
                },
                ontology_terms=[ontology],
            )

            # track count
            alt_vals = vout.alternate.splice_sites.values
            ref_vals = vout.reference.splice_sites.values
            n_tracks = alt_vals.shape[1]
            if n_tracks == 0:
                warnings.warn(f'No tracks available for {ontology}; skipping variant {variant}.')
                continue

            # indel alignment
            length_alter = len(variant.reference_bases) - len(variant.alternate_bases)
            if length_alter != 0:
                align_reference_for_indel(variant, interval, vout, length_alter)
                ref_vals = vout.reference.splice_sites.values  # updated

            # scan window bounds (0-based indexes into interval arrays)
            center_idx = variant.position - interval.start
            start_idx = center_idx - scan_span
            end_idx = center_idx + scan_span  # exclusive after adjust in compute

            # compute window scores (vectorized)
            win_scores = compute_window_scores(alt_vals, ref_vals, start_idx, end_idx)
            if win_scores.size == 0:
                warnings.warn(f'Scan span/window too large/small near edges for {variant}.')
                continue

            # track names
            try:
                track_names = (
                    f'{vout.reference.splice_sites.metadata.name}:{vout.reference.splice_sites.metadata.strand}'
                )
            except Exception:
                track_names = [f'track_{i}' for i in range(n_tracks)]

            # scan all tracks always; early-stop if user did NOT set scan_all_tracks? We'll still scan for table,
            # but we skip plotting after first significant if flag not set (speed).
            for ti in range(n_tracks):
                scores = win_scores[:, ti]
                regions = call_regions(scores)
                is_sig = len(regions) > 0

                # gather rows (one per region; if none, n_regions=0 + region_index=-1)
                if is_sig:
                    for ri, (rs, re) in enumerate(regions):
                        # convert window index -> bp relative to variant
                        rel_start = (start_idx + rs) - center_idx
                        rel_end = (start_idx + re + window_size - 1) - center_idx
                        abs_start = variant.position + rel_start
                        abs_end = variant.position + rel_end
                        seg_scores = scores[rs : re + 1]
                        mean_s = float(np.nanmean(seg_scores))
                        max_s = float(np.nanmax(seg_scores))
                        min_s = float(np.nanmin(seg_scores))
                        if mean_s > 0 and min_s > 0:
                            direction = 'up'
                        elif mean_s < 0 and max_s < 0:
                            direction = 'down'
                        else:
                            direction = 'mixed'
                        results_rows.append(
                            {
                                'chrom': variant.chromosome,
                                'pos': variant.position,
                                'ref': variant.reference_bases,
                                'alt': variant.alternate_bases,
                                'ontology': ontology,
                                'track_name': track_names[ti] if ti < len(track_names) else f'track_{ti}',
                                'is_significant': True,
                                'n_regions': len(regions),
                                'region_index': ri,
                                'rel_start_bp': int(rel_start),
                                'rel_end_bp': int(rel_end),
                                'abs_start_bp': int(abs_start),
                                'abs_end_bp': int(abs_end),
                                'mean_score': mean_s,
                                'max_score': max_s,
                                'min_score': min_s,
                                'direction': direction,
                                'plot_file': None,  # filled after plotting variant-level below
                            }
                        )
                else:
                    results_rows.append(
                        {
                            'chrom': variant.chromosome,
                            'pos': variant.position,
                            'ref': variant.reference_bases,
                            'alt': variant.alternate_bases,
                            'ontology': ontology,
                            'track_name': track_names[ti] if ti < len(track_names) else f'track_{ti}',
                            'is_significant': False,
                            'n_regions': 0,
                            'region_index': -1,
                            'rel_start_bp': np.nan,
                            'rel_end_bp': np.nan,
                            'abs_start_bp': np.nan,
                            'abs_end_bp': np.nan,
                            'mean_score': float(np.nanmean(scores)),
                            'max_score': float(np.nanmax(scores)),
                            'min_score': float(np.nanmin(scores)),
                            'direction': 'none',
                            'plot_file': None,
                        }
                    )

            # decide plot size (reuse original heuristic)
            if abs(length_alter) >= 2**14:
                plot_size = abs(length_alter) * 4
            else:
                plot_size = 2**15

            # plot organ-level REF/ALT overlay once per variant-organ if any sig OR user asked plot_non_sig
            sig_any = any(
                (
                    r['is_significant']
                    and r['chrom'] == variant.chromosome
                    and r['pos'] == variant.position
                    and r['ontology'] == ontology
                )
                for r in results_rows[-n_tracks:]
            )
            if sig_any or args.plot_non_sig:
                plot_path = output_as_path / (f'{ontology}_{number_rank}_{variant.chromosome}_{variant.position}.png')
                plot_variant_tracks(variant, interval, vout, transcript_extractor, plot_size, plot_path)
                # back-fill plot_file for recent rows
                for r in results_rows[-n_tracks:]:
                    if r['chrom'] == variant.chromosome and r['pos'] == variant.position and r['ontology'] == ontology:
                        r['plot_file'] = str(plot_path)

    # --------------------------
    # build DataFrame & write
    # --------------------------
    if not results_rows:
        print('No results generated.')
        return 0

    res_df = pd.DataFrame(results_rows)
    detailed_results = str(output_as_path / 'alphagenome_scan_results.csv')
    write_table(res_df, detailed_results)

    # variant×organ summary
    agg = (
        res_df.groupby(['chrom', 'pos', 'ref', 'alt', 'ontology'])['is_significant']
        .any()
        .reset_index()
        .rename(columns={'is_significant': 'is_significant_any'})
    )
    # track list for each variant×organ
    sig_tracks = (
        res_df[res_df.is_significant]
        .groupby(['chrom', 'pos', 'ref', 'alt', 'ontology'])['track_name']
        .apply(lambda s: ','.join(sorted(set(s))))
        .reset_index()
    )
    agg = agg.merge(sig_tracks, how='left', on=['chrom', 'pos', 'ref', 'alt', 'ontology'])
    agg['track_name'] = agg['track_name'].fillna('')

    summary_output_path = str(output_as_path / 'alphagenome_scan_results_variant_organ_summary.csv')
    write_table(agg, summary_output_path)

    print('Done.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
