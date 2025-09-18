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
import asyncio
from argparse import Namespace
from typing import List, AsyncGenerator
import json
from argparse import ArgumentParser
from csv import DictReader
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from cloudpathlib.anypath import to_anypath

# Load from JSON file
with open('src/ontology_to_biosample_mapping.json') as file:
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
        self.output_root = output_root
        self.organization_mode = organization_mode
        self.default_view_limits = '-1:1'
        self.positive_color = '0,0,255'  # blue
        self.negative_color = '220,20,60'  # red

        # Create main bedgraph directory
        self.bedgraph_dir = Path(output_root) / Path(var_file).stem
        self.bedgraph_dir.mkdir(parents=True, exist_ok=True)

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

    def get_output_path(self, variant, track_name: str, track_type: str, ontology: str, output_type: str) -> Path:
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
            variant_dir.mkdir(exist_ok=True)

            # Create ontology directory within variant directory using global ontology_to_biosample
            ontology_dir_name = f'{ontology}:{ontology_to_biosample[ontology]}'
            ontology_dir = variant_dir / ontology_dir_name
            ontology_dir.mkdir(exist_ok=True)

            # Create specific directory within variant directory
            specific_dir = ontology_dir / dir_mappings.get(output_type, output_type)
            specific_dir.mkdir(exist_ok=True)

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
        track_type_dir.mkdir(exist_ok=True)

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
        variant: Any,
        interval: Any,
        values: np.ndarray,
        track_name: str,
        output_type: str,
        ontology: str,
        track_type: str = 'reference',
    ) -> Path:
        """Write a complete BedGraph file with header and data.

        Creates a properly formatted BedGraph file with browser directives,
        track definition line, and compressed genomic data.

        Args:
            variant: Variant being analyzed
            interval: Genomic interval for the data
            values: Array of numeric values for each position
            track_name: Display name for the track
            output_type: Type of output (rna_seq, splice_sites, etc.)
            ontology: Ontology term for organization
            track_type: Track category (reference, alternate, delta)

        Returns:
            Path to the created BedGraph file

        Raises:
            IOError: If file cannot be written
        """
        output_path = self.get_output_path(variant, track_name, track_type, ontology, output_type)
        compressed_data = self._compress_values(values, interval.chromosome, interval.start)

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

    def generate_all_tracks(self, variant, interval, scores, ontology: str) -> list[Path]:
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

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate reference sequence track."""
        ref_series = self._get_values_array(scores.reference.rna_seq.values, index)
        return self.writer.write_bedgraph(variant, interval, ref_series, track_name, 'rna_seq', ontology, 'reference')

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.rna_seq.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        return self.writer.write_bedgraph(variant, interval, alt_series, track_name, 'rna_seq', ontology, 'alternate')

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate delta (difference) track."""
        ref_series = self._get_values_array(scores.reference.rna_seq.values, index)
        alt_raw = [val[index] for val in scores.alternate.rna_seq.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(variant, interval, delta_series, track_name, 'rna_seq', ontology, 'delta')


class SpliceSiteUsageTrackGenerator(BaseTrackGenerator):
    """Generates splice site usage tracks for variants."""

    def generate_all_tracks(self, variant, interval, scores, ontology: str) -> list[Path]:
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

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate reference sequence track."""
        ref_series = self._get_values_array(scores.reference.splice_sites.values, index)
        return self.writer.write_bedgraph(
            variant, interval, ref_series, track_name, 'splice_site_usage', ontology, 'reference'
        )

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        return self.writer.write_bedgraph(
            variant, interval, alt_series, track_name, 'splice_site_usage', ontology, 'alternate'
        )

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
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

    def generate_all_tracks(self, variant, interval, scores, ontology: str) -> list[Path]:
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

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate reference sequence track."""
        ref_series = self._get_values_array(scores.reference.splice_sites.values, index)
        return self.writer.write_bedgraph(
            variant, interval, ref_series, track_name, 'splice_sites', ontology, 'reference'
        )

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]
        alt_series = self._process_track_values(variant, interval, alt_raw)
        return self.writer.write_bedgraph(
            variant, interval, alt_series, track_name, 'splice_sites', ontology, 'alternate'
        )

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str) -> Path:
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

    def generate_junction_bed_tracks(self, variant, interval, scores, ontology: str) -> list[Path]:
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
        generated_files: list[Path] = []

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
    ) -> Path:
        """Generate a single BED track file for junctions."""
        output_path = self.writer.get_output_path(variant, track_name, track_type, ontology, 'junctions')

        with to_anypath(output_path).open('w') as f:
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


async def load_variants_in_batches(var_file: str, batch_size: int = 30) -> AsyncGenerator[List, None]:
    """Load variants in batches to avoid memory issues."""
    variants = load_variants_table(var_file)
    for i in range(0, len(variants), batch_size):
        yield variants[i : i + batch_size]


async def process_variant_batch(model, variants_batch):
    """Process a batch of variants concurrently."""
    semaphore = asyncio.Semaphore(30)  # Limit to 30 concurrent operations

    async def score_variant(var):
        async with semaphore:  # Ensure only 30 run at once
            loop = asyncio.get_event_loop()
            interval = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

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

    # Process all variants in batch concurrently
    tasks = [score_variant(var) for var in variants_batch]
    results = await asyncio.gather(*tasks)
    return [r for r in results if r is not None]


async def main_async(var_file: str, output_root: str, ontologies: List[str], organization: str = 'variant'):
    """Memory-efficient async main function."""
    model = dna_client.create(api_key)
    variants = load_variants_table(var_file)
    print(f'Loaded {len(variants)} variants.')

    # Process one batch at a time to keep memory usage low
    all_significant_data = []  # Only store minimal data

    batch_count = 0
    async for batch in load_variants_in_batches(var_file):
        batch_count += 1
        print(f'Processing batch {batch_count}')

        batch_results = await process_variant_batch(model, batch)

        # Store only essential data
        for scores_df, interval, var in batch_results:
            all_significant_data.append({'scores': scores_df, 'interval': interval, 'variant': var})

        # Optional: Clear batch from memory
        del batch

    if not all_significant_data:
        print('No significant variants found.')
        return None

    # Combine only when needed for final output
    significant_results = pd.concat([item['scores'] for item in all_significant_data], ignore_index=True)
    significant_results.to_csv(f'{output_root}.tsv', sep='\t', index=False)

    print(f'Processed {len(all_significant_data)} significant variants.')

    # Generate tracks one batch at a time
    all_generated_files = []
    batch_size = 30

    for i in range(0, len(all_significant_data), batch_size):
        batch_data = all_significant_data[i : i + batch_size]
        intervals_batch = [item['interval'] for item in batch_data]
        vars_batch = [item['variant'] for item in batch_data]

        # Process each ontology for this batch
        for ontology in ontologies:
            variant_predictions = model.predict_variants(
                intervals=intervals_batch,
                variants=vars_batch,
                ontology_terms=[ontology],
                requested_outputs={dna_client.OutputType.SPLICE_SITES},
            )
            # Generate files...

    return all_generated_files


def run_main(args: Namespace):
    """Run the async main function."""
    asyncio.run(
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
    run_main(args)
