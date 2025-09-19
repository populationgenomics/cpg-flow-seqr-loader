"""
AlphaGenome BedGraph Plotter

This module generates browser-compatible genomic tracks from AlphaGenome API predictions
to visualize the functional impact of splice site variants. It creates BedGraph files
for splice sites, splice site usage, RNA-seq data, and BED files for junction visualization.

Example:
    python alphagenome_bedgraph_plotter.py \
        --var_file variants.tsv \
        --output_root ./tracks \
        --ontology UBERON:0001134 \
        --organization variant

Required TSV columns: CHROM, POS, REF, ALT
"""

# ruff: noqa: PD011, PLR0915, ARG002

import json
from csv import DictReader
from enum import Enum
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from cloudpathlib.anypath import to_anypath
import asyncio
import pandas as pd
from argparse import ArgumentParser, Namespace
from typing import List
import io


# Load from JSON file
with to_anypath('src/ontology_to_biosample_mapping.json').open() as file:
    ontology_to_biosample = json.load(file)


class OrganizationMode(Enum):
    """Output organization modes for generated tracks.

    BY_VARIANT: Group files by variant, then by track type
    BY_TRACK_TYPE: Group files by track type (reference/alternate/delta)
    """

    BY_VARIANT = 'variant'
    BY_TRACK_TYPE = 'track_type'


THRESHOLD = 0.01
api_key = 'AIzaSyDYx7VMDPepU7qeJOm7i-AVm9GsrV-BbW8'


def load_variants_table(path: str) -> list[genome.Variant]:
    """Load variants from a TSV file.

    Reads a tab-separated file containing variant information and converts
    each row into a genome.Variant object.

    Args:
        path: Path to the TSV file. Must contain columns: CHROM, POS, REF, ALT

    Returns:
        List of genome.Variant objects parsed from the input file

    Raises:
        FileNotFoundError: If the input file doesn't exist
        KeyError: If required columns are missing from the TSV
        ValueError: If POS column contains non-integer values

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


class BedGraphWriter:
    """Handles creation and writing of BedGraph files with track metadata.

    This class manages the output directory structure and generates properly
    formatted BedGraph files with browser-compatible headers and styling.

    Attributes:
        output_root: Root directory for all output files
        organization_mode: How to organize output files (by variant or track type)
        bedgraph_dir: Main directory for BedGraph files
        track_colors: Color mapping for different track types
    """

    def __init__(
        self,
        var_file: str,
        output_root: str,
        ontology: str,
        organization_mode: OrganizationMode = OrganizationMode.BY_VARIANT,
    ):
        self.output_root = to_anypath(output_root)
        self.organization_mode = organization_mode
        self.default_view_limits = '-1:1'
        self.positive_color = '0,0,255'  # blue
        self.negative_color = '220,20,60'  # red

        # Create main bedgraph directory
        var_file_path = to_anypath(var_file)
        self.bedgraph_dir = self.output_root / var_file_path.stem

        # Color mapping for track types
        self.track_colors = {
            'reference': ('128,128,128', '128,128,128'),  # grey
            'alternate': ('169,0,17', '169,0,17'),  # red
            'delta': ('255,192,0', '255,192,0'),  # yellow
        }

    """Initialize BedGraphWriter with output configuration.

    Args:
        var_file: Path to variant file (used for directory naming)
        output_root: Root directory for all output files
        ontology: Ontology term being processed
        organization_mode: How to organize output files
    """

    def get_output_path(self, variant, track_name: str, track_type: str, ontology: str, output_type: str):
        """Determine output path based on organization mode and output type."""
        # Define directory structure mappings
        dir_mappings = {
            'rna_seq': 'RNAseq',
            'splice_sites': 'splice_sites',
            'splice_site_usage': 'splice_site_usage',
            'junctions': 'junctions',
        }

        if self.organization_mode == OrganizationMode.BY_VARIANT:
            # Group by sample (variant)
            variant_dir = self.bedgraph_dir / f'{variant!s}'

            # Create ontology directory within variant directory
            ontology_dir_name = f'{ontology}:{ontology_to_biosample[ontology]}'
            ontology_dir = variant_dir / ontology_dir_name

            # Create specific directory within variant directory
            specific_dir = ontology_dir / dir_mappings.get(output_type, output_type)

            if output_type == 'junctions':
                filename = f'{track_type}_{track_name}.bed'
            else:
                filename = (
                    f'{output_type}_{track_type}_{track_name}.bedgraph'
                    if output_type == 'rna_seq'
                    else f'{track_type}_{track_name}.bedgraph'
                )
            return specific_dir / filename

        # BY_TRACK_TYPE
        # Group by track type (reference, alternate, delta)
        track_type_dir = self.bedgraph_dir / track_type

        if output_type == 'junctions':
            filename = f'{variant!s}_junctions_{track_name}.bed'
        else:
            filename = (
                f'{variant!s}_RNA_{track_name}.bedgraph'
                if output_type == 'rna_seq'
                else f'{variant!s}_{track_name}.bedgraph'
            )
        return track_type_dir / filename

    def write_bedgraph(
        self,
        variant: any,
        interval: any,
        values: np.ndarray,
        track_name: str,
        output_type: str,
        ontology: str,
        track_type: str = 'reference',
    ):
        """Write a complete BedGraph file with header and data."""
        output_path = self.get_output_path(variant, track_name, track_type, ontology, output_type)
        compressed_data = self._compress_values(values, interval.chromosome, interval.start)

        # For cloud storage, just write the file - "directories" are created automatically
        with to_anypath(output_path).open('w') as f:
            # Write header
            for line in self._create_header(interval):
                f.write(f'{line}\n')

            # Write track definition
            f.write(self._create_track_line(track_name, track_type))

            # Write data
            for chrom, start, end, value in compressed_data:
                f.write(f'{chrom}\t{start}\t{end}\t{value:.6f}\n')

        return output_path

    def _create_header(self, interval) -> list[str]:
        """Create BedGraph header with browser settings."""
        return [
            f'browser position {interval.chromosome}:{interval.start}-{interval.end}',
            'browser hide all',
            'browser pack refGene',
        ]

    def _create_track_line(self, track_name: str, track_type: str = 'reference') -> str:
        """Create track definition line with styling."""
        # Get colors based on track type
        color, alt_color = self.track_colors.get(track_type, (self.positive_color, self.negative_color))

        return (
            f'track type=bedGraph '
            f'name="{track_name}" '
            f'description="{track_name} ({track_type})" '
            f'visibility=full '
            f'viewLimits={self.default_view_limits} '
            f'autoScale=off '
            f'color={color} '
            f'altColor={alt_color}\n'
        )

    def _compress_values(self, values: np.ndarray, chrom: str, start: int) -> list[tuple[str, int, int, float]]:
        """Compress adjacent similar values into regions to reduce file size.

        Combines consecutive positions with similar values (within 2% of data range)
        into single regions, significantly reducing file sizes while preserving
        visual fidelity.

        Args:
            values: Array of values to compress
            chrom: Chromosome name
            start: Starting genomic position

        Returns:
            List of tuples: (chromosome, start, end, value) for each region
        """
        if len(values) == 0:
            return []

        # Handle all NaN case
        if np.all(np.isnan(values)):
            return []

        # Calculate compression threshold
        valid_values = values[~np.isnan(values)]
        if len(valid_values) == 0:
            return []

        data_range = np.nanmax(valid_values) - np.nanmin(valid_values)
        min_delta = data_range * 0.02 if data_range > 0 else 1e-10

        current_region_start = start
        current_region_value = values[0]
        rows_for_file = []
        final_pos = start

        for position, value in enumerate(values[1:], start=start + 1):
            if np.isnan(value):
                continue

            if abs(value - current_region_value) > min_delta:
                if not np.isnan(current_region_value):
                    rows_for_file.append((chrom, current_region_start, position, current_region_value))
                current_region_start = position
                current_region_value = value

            final_pos = position

        # Add final region
        if not np.isnan(current_region_value):
            rows_for_file.append((chrom, current_region_start, final_pos + 1, current_region_value))

        return rows_for_file


class BaseTrackGenerator:
    """Base class for track generators with shared functionality.

    Provides common methods for value processing, indel realignment,
    and array manipulation used by all track generator subclasses.
    """

    def __init__(self, bedgraph_writer: BedGraphWriter):
        self.writer = bedgraph_writer

    def _realign_indel_values(self, variant, interval, values: list[float]) -> list[float]:
        """Realign values for indel variants to maintain coordinate mapping.

        For insertions: removes values at insertion site and pads with NaN
        For deletions: inserts NaN values to maintain coordinate alignment

        Args:
            variant: Variant with potential indel
            interval: Genomic interval
            values: Original values array

        Returns:
            Realigned values array maintaining original length
        """
        ref_len = len(variant.reference_bases)
        alt_len = len(variant.alternate_bases)

        if ref_len == alt_len:
            return values

        alter_by = abs(alt_len - ref_len)
        original_length = len(values)
        midpoint = variant.position - interval.start

        if alt_len > ref_len:
            # Insertion
            if midpoint + alter_by <= len(values):
                mean_val = np.nanmean(values[midpoint : midpoint + alter_by])
                values = values[:midpoint] + values[midpoint + alter_by :] + [np.nan] * alter_by
                if midpoint < len(values):
                    values[midpoint] = mean_val
        else:
            # Deletion
            values = np.insert(values, midpoint, [np.nan] * alter_by).tolist()

        return values[:original_length]

    def _get_values_array(self, values_data, index: int) -> np.ndarray:
        """Extract values array from data structure."""
        return np.array([val[index] for val in values_data])

    def _process_track_values(self, variant, interval, values_raw: list[float]) -> np.ndarray:
        """Process track values with indel realignment if needed."""
        if len(variant.reference_bases) != len(variant.alternate_bases):
            values_processed = self._realign_indel_values(variant, interval, values_raw)
        else:
            values_processed = values_raw
        return np.array(values_processed)


class RNAseqTrackGenerator(BaseTrackGenerator):
    """Generates RNA-seq tracks for variant analysis.

    Creates reference, alternate, and delta (difference) tracks for all
    available RNA-seq datasets from the AlphaGenome predictions.
    """

    def generate_all_tracks(self, variant, interval, scores, ontology: str) -> list:
        """Generate all RNA-seq tracks for a variant.

        Creates reference, alternate, and delta tracks for each RNA-seq
        dataset present in the scores.

        Args:
            variant: Variant being analyzed
            interval: Genomic interval
            scores: AlphaGenome prediction scores
            ontology: Ontology term for file organization

        Returns:
            List of paths to generated BedGraph files
        """
        if not hasattr(scores.reference, 'rna_seq'):
            return []

        generated_files = []

        for index in range(scores.reference.rna_seq.num_tracks):
            track_name = scores.reference.rna_seq.names[index]
            track_strand = scores.reference.rna_seq.strands[index]
            full_track_name = f'{track_name} {track_strand}'

            # Generate all three track types
            ref_file = self._generate_reference_track(variant, interval, scores, index, full_track_name, ontology)
            alt_file = self._generate_alternate_track(variant, interval, scores, index, full_track_name, ontology)
            delta_file = self._generate_delta_track(variant, interval, scores, index, full_track_name, ontology)

            generated_files.extend([ref_file, alt_file, delta_file])

        return generated_files

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate reference sequence track."""
        ref_series = self._get_values_array(scores.reference.rna_seq.values, index)
        return self.writer.write_bedgraph(variant, interval, ref_series, track_name, 'rna_seq', ontology, 'reference')

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.rna_seq.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        return self.writer.write_bedgraph(variant, interval, alt_series, track_name, 'rna_seq', ontology, 'alternate')

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate delta (difference) track."""
        ref_series = self._get_values_array(scores.reference.rna_seq.values, index)
        alt_raw = [val[index] for val in scores.alternate.rna_seq.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(variant, interval, delta_series, track_name, 'rna_seq', ontology, 'delta')


class SpliceSiteUsageTrackGenerator(BaseTrackGenerator):
    """Generates splice site usage tracks for variants."""

    def generate_all_tracks(self, variant, interval, scores, ontology: str) -> list:
        """Generate reference, alternate, and delta tracks for all splice site usage."""
        if not hasattr(scores.reference, 'splice_sites'):
            return []

        generated_files = []

        for index in range(scores.reference.splice_sites.num_tracks):
            track_name = scores.reference.splice_sites.names[index]
            track_strand = scores.reference.splice_sites.strands[index]
            full_track_name = f'Splice Site Usage {track_name} {track_strand} {ontology}'

            # Generate all three track types
            ref_file = self._generate_reference_track(variant, interval, scores, index, full_track_name, ontology)
            alt_file = self._generate_alternate_track(variant, interval, scores, index, full_track_name, ontology)
            delta_file = self._generate_delta_track(variant, interval, scores, index, full_track_name, ontology)

            generated_files.extend([ref_file, alt_file, delta_file])

        return generated_files

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate reference sequence track."""
        ref_series = self._get_values_array(scores.reference.splice_sites.values, index)
        return self.writer.write_bedgraph(
            variant, interval, ref_series, track_name, 'splice_site_usage', ontology, 'reference'
        )

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        return self.writer.write_bedgraph(
            variant, interval, alt_series, track_name, 'splice_site_usage', ontology, 'alternate'
        )

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate delta (difference) track."""
        ref_series = self._get_values_array(scores.reference.splice_sites.values, index)
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(
            variant, interval, delta_series, track_name, 'splice_site_usage', ontology, 'delta'
        )


class SpliceSiteTrackGenerator(BaseTrackGenerator):
    """Generates splice site tracks for variants."""

    def generate_all_tracks(self, variant, interval, scores, ontology: str) -> list:
        """Generate reference, alternate, and delta tracks for all splice sites."""
        if not hasattr(scores.reference, 'splice_sites'):
            return []

        generated_files = []

        for index in range(scores.reference.splice_sites.num_tracks):
            track_name = scores.reference.splice_sites.names[index]
            track_strand = scores.reference.splice_sites.strands[index]
            full_track_name = f'Splice Sites {track_strand} {track_name}'

            # Generate all three track types
            ref_file = self._generate_reference_track(variant, interval, scores, index, full_track_name, ontology)
            alt_file = self._generate_alternate_track(variant, interval, scores, index, full_track_name, ontology)
            delta_file = self._generate_delta_track(variant, interval, scores, index, full_track_name, ontology)

            generated_files.extend([ref_file, alt_file, delta_file])

        return generated_files

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate reference sequence track."""
        ref_series = self._get_values_array(scores.reference.splice_sites.values, index)
        return self.writer.write_bedgraph(
            variant, interval, ref_series, track_name, 'splice_sites', ontology, 'reference'
        )

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        return self.writer.write_bedgraph(
            variant, interval, alt_series, track_name, 'splice_sites', ontology, 'alternate'
        )

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate delta (difference) track."""
        ref_series = self._get_values_array(scores.reference.splice_sites.values, index)
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(
            variant, interval, delta_series, track_name, 'splice_sites', ontology, 'delta'
        )


class JunctionTrackGenerator:
    """Generates BED tracks for splice junction visualization.

    Creates BED12 format files suitable for sashimi plot visualization
    in genome browsers, with junction arcs and expression values.
    """

    def __init__(self, bedgraph_writer: BedGraphWriter):
        self.writer = bedgraph_writer

    def generate_junction_bed_tracks(self, variant, interval, scores, ontology: str) -> list:
        """Generate BED files for junction visualization from JunctionData.
                Creates BED12 format files with junction data suitable for
        genome browser sashimi plot visualization.

        Args:
            variant: Variant being analyzed
            interval: Genomic interval
            scores: AlphaGenome prediction scores containing junction data
            ontology: Ontology term for file organization

        Returns:
            List of paths to generated BED files

        Note:
            Generated files use BED12 format with:
            - Block structure for junction visualization
            - Score values scaled for line thickness
            - Color coding by track type
        """
        generated_files: list = []

        if not hasattr(scores.reference, 'splice_junctions'):
            return generated_files

        junction_data = scores.reference.splice_junctions

        # Generate tracks for each tissue/track in metadata
        for track_index in range(junction_data.num_tracks):
            track_name = junction_data.metadata.iloc[track_index]['name']

            # Generate reference junction track
            ref_bed_file = self._generate_junction_bed_track(
                variant, interval, junction_data, track_index, track_name, 'reference', ontology
            )
            generated_files.append(ref_bed_file)

            # Generate alternate junction track if available
            if hasattr(scores.alternate, 'splice_junctions'):
                alt_junction_data = scores.alternate.splice_junctions
                alt_bed_file = self._generate_junction_bed_track(
                    variant, interval, alt_junction_data, track_index, track_name, 'alternate', ontology
                )
                generated_files.append(alt_bed_file)

        return generated_files

    def _generate_junction_bed_track(
            self, variant, interval, junction_data, track_index: int, track_name: str, track_type: str, ontology: str
    ):
        """Generate a single BED track file for junctions."""
        output_path = self.writer.get_output_path(variant, track_name, track_type, ontology, 'junctions')

        # Ensure parent directories exist for both local and cloud paths

        with output_path.open('w') as f:
            # Write BED track header for sashimi visualization with score labels
            f.write(f'track type=bed name="Sashimi {track_type} {ontology}" ')
            f.write(f'description="Junction {track_name} ({track_type})" ')
            f.write('graphType=junctions itemRgb=on visibility=full ')
            f.write('scoreFilter=0 useScore=1 scoreLabel=1 ')
            f.write('maxHeightPixels=128:60:11\n')

            # Process each junction
            for junction_idx, junction in enumerate(junction_data.junctions):
                score_value = junction_data.values[junction_idx, track_index]

                # Skip low-value junctions
                if score_value < 0.01 or np.isnan(score_value):
                    continue

                # Scale score significantly for visibility
                bed_score = min(1000, max(1, int(score_value * 10000)))

                # Junction coordinates
                chrom = junction.chromosome
                start = junction.start
                end = junction.end
                strand = junction.strand

                # Create junction ID that will be displayed as label
                junction_id = f'{score_value:.2f}'

                # Create blocks for junction visualization
                block_size = 1
                block_starts = f'0,{end - start - block_size}'
                block_sizes = f'{block_size},{block_size}'

                # Color coding
                color = '128,128,128' if track_type == 'reference' else '169,0,17'

                # Write BED12 format line
                bed_line = '\t'.join(
                    [
                        chrom,  # chromosome
                        str(start),  # chromStart
                        str(end),  # chromEnd
                        junction_id,  # name (score value for display)
                        str(bed_score),  # score (scaled up for thickness)
                        strand,  # strand
                        str(start),  # thickStart
                        str(end),  # thickEnd
                        color,  # itemRgb
                        '2',  # blockCount
                        block_sizes,  # blockSizes
                        block_starts,  # blockStarts
                    ]
                )

                f.write(bed_line + '\n')

        return output_path


def process_variant_with_globals(args):
    """Process a single variant with global generators"""
    variant_prediction, var, interval, ontology = args

    print(f'Processing variant: {var!s} in interval {interval.chromosome}:{interval.start}-{interval.end}')

    generated_files = []
    generated_files.extend(global_track_generator.generate_all_tracks(var, interval, variant_prediction, ontology))
    generated_files.extend(global_rna_track_generator.generate_all_tracks(var, interval, variant_prediction, ontology))
    generated_files.extend(
        global_usage_track_generator.generate_all_tracks(var, interval, variant_prediction, ontology)
    )
    generated_files.extend(
        global_junction_track_generator.generate_junction_bed_tracks(var, interval, variant_prediction, ontology)
    )

    return generated_files


def init_worker(track_gen, rna_gen, junction_gen, usage_gen):
    """Initialize worker process with shared objects"""
    global \
        global_track_generator, \
        global_rna_track_generator, \
        global_junction_track_generator, \
        global_usage_track_generator
    global_track_generator = track_gen
    global_rna_track_generator = rna_gen
    global_junction_track_generator = junction_gen
    global_usage_track_generator = usage_gen



async def process_variants_streaming(model, var_file: str, output_root: str):
    """Process variants one at a time to minimize memory usage."""
    variants = load_variants_table(var_file)
    print(f'Loaded {len(variants)} variants.')

    significant_count = 0
    all_significant_dfs = []  # Only store DataFrames, not full variant objects

    # Process variants one by one (or in small groups)
    semaphore = asyncio.Semaphore(30)  # Limit concurrent processing

    async def process_single_variant(var):
        async with semaphore:
            loop = asyncio.get_event_loop()
            interval = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

            # Run blocking operation in thread pool
            variant_scores = await loop.run_in_executor(
                None,
                lambda: model.score_variant(
                    interval=interval,
                    variant=var,
                    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITE_USAGE']],
                ),
            )

            tidied_scores = variant_scorers.tidy_scores([variant_scores], match_gene_strand=True)
            filtered_scores = tidied_scores[tidied_scores['raw_score'].values >= THRESHOLD]

            if not filtered_scores.empty:
                return filtered_scores, interval, var
            return None

    # Process variants in small batches to control memory
    batch_size = 10  # Smaller batches
    for i in range(0, len(variants), batch_size):
        batch = variants[i : i + batch_size]
        tasks = [process_single_variant(var) for var in batch]
        results = await asyncio.gather(*tasks)

        # Collect only the essential data
        for result in results:
            if result is not None:
                scores_df, interval, var = result
                all_significant_dfs.append(scores_df)
                significant_count += 1
                # Don't store the full variant object or interval unless needed

        # Optional: Force garbage collection
        del batch, tasks, results

    print(f'{significant_count} variants had significant scores.')

    # Only combine DataFrames at the end
    if all_significant_dfs:
        final_df = pd.concat(all_significant_dfs, ignore_index=True)

        # Save to file
        buffer = io.StringIO()
        final_df.to_csv(buffer, sep='\t', index=False)
        buffer.seek(0)

        output_path = to_anypath(f'{output_root}.tsv')
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('w') as f:
            f.write(buffer.read())

        buffer.close()
        del final_df  # Free memory

    return significant_count


async def generate_tracks_streaming(
    model, var_file: str, output_root: str, ontologies: List[str], significant_count: int
):
    """Generate tracks in a memory-efficient way."""
    if significant_count == 0:
        return []

    variants = load_variants_table(var_file)

    # Process track generation in small batches
    batch_size = 5  # Even smaller for track generation
    all_generated_files = []

    for i in range(0, len(variants), batch_size):
        batch = variants[i : i + batch_size]

        # Score this batch first to filter significant ones
        significant_in_batch = []
        for var in batch:
            interval = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            # Quick check if this variant was significant (you might want to optimize this)
            significant_in_batch.append((var, interval))

        if not significant_in_batch:
            continue

        # Process each ontology for this batch
        for ontology in ontologies:
            vars_batch = [item[0] for item in significant_in_batch]
            intervals_batch = [item[1] for item in significant_in_batch]

            # Generate predictions for this small batch
            loop = asyncio.get_event_loop()
            variant_predictions = await loop.run_in_executor(
                None,
                lambda: model.predict_variants(
                    intervals=intervals_batch,
                    variants=vars_batch,
                    ontology_terms=[ontology],
                    requested_outputs={dna_client.OutputType.SPLICE_SITES},
                ),
            )

            # Generate track files immediately and collect paths
            # (implement your track generation logic here)
            # Don't store the full prediction data, just generate files

    return all_generated_files


async def main_async(var_file: str, output_root: str, ontologies: List[str], organization: str = 'variant'):
    """Memory-efficient async main function."""
    model = dna_client.create(api_key)

    # Step 1: Process variants and find significant ones
    significant_count = await process_variants_streaming(model, var_file, output_root)

    # Step 2: Generate tracks for significant variants only
    generated_files = await generate_tracks_streaming(model, var_file, output_root, ontologies, significant_count)

    return generated_files


def run_main(args: Namespace):
    """Run the async main function."""
    return asyncio.run(
        main_async(
            var_file=args.var_file,
            output_root=args.output_root,
            ontologies=args.ontology,
            organization=args.organization,
        )
    )


if __name__ == '__main__':
    parser = ArgumentParser(description='Generate BedGraph tracks for splice site variants')
    parser.add_argument('--var_file', help='Path to variant TSV file', required=True)
    parser.add_argument('--output_root', help='Root output directory', required=True)
    parser.add_argument(
        '--ontology', nargs='+', default=['UBERON:0001134', 'UBERON:0002113', 'UBERON:0002369'], help='Ontology terms'
    )
    parser.add_argument(
        '--organization',
        choices=['variant', 'track_type'],
        default='variant',
        help="Organization mode: 'variant' groups by variant, 'track_type' groups by reference/alternate/delta",
    )
    args = parser.parse_args()
    result = run_main(args)
    print('Generated files:', result)