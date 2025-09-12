import io
from argparse import ArgumentParser
from csv import DictReader
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm.auto
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from cloudpathlib.anypath import to_anypath
from tqdm import tqdm


class OrganizationMode(Enum):
    BY_VARIANT = 'variant'
    BY_TRACK_TYPE = 'track_type'


API_KEY = 'get yer own'


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


class BedGraphWriter:
    """Handles creation and writing of BedGraph files with track metadata."""

    def __init__(
        self,
        var_file,
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

    def _get_splice_sites_output_path(self, variant, track_name: str, track_type: str, ontology: str) -> Path:
        """Determine output path based on organization mode."""
        if self.organization_mode == OrganizationMode.BY_VARIANT:
            # Group by sample (variant)
            variant_dir = self.bedgraph_dir / f'{variant!s}'
            variant_dir.mkdir(exist_ok=True)

            # Create track_type directory within variant directory
            ontology_dir = variant_dir / ontology
            ontology_dir.mkdir(exist_ok=True)

            # Create track_type directory within variant directory
            track_type_dir = ontology_dir / 'splice_sites'
            track_type_dir.mkdir(exist_ok=True)

            filename = f'splice_sites_{track_type}_{track_name}.bedgraph'
            return track_type_dir / filename
        if self.organization_mode == OrganizationMode.BY_TRACK_TYPE:
            # Group by track type (reference, alternate, delta)
            track_type_dir = self.bedgraph_dir / track_type
            track_type_dir.mkdir(exist_ok=True)
            filename = f'{variant!s}_splice_sites_{track_name}.bedgraph'
            return track_type_dir / filename

    def splice_site_usage_outpath(self, variant, track_name: str, track_type: str, ontology: str) -> Path:
        """Determine output path based on organization mode."""
        if self.organization_mode == OrganizationMode.BY_VARIANT:
            # Group by sample (variant)
            variant_dir = self.bedgraph_dir / f'{variant!s}'
            variant_dir.mkdir(exist_ok=True)

            # Create track_type directory within variant directory
            ontology_dir = variant_dir / ontology
            ontology_dir.mkdir(exist_ok=True)

            # Create track_type directory within variant directory
            track_type_dir = ontology_dir / 'splice_site_usage'
            track_type_dir.mkdir(exist_ok=True)

            filename = f'{track_type}_{track_name}.bedgraph'
            return track_type_dir / filename

        if self.organization_mode == OrganizationMode.BY_TRACK_TYPE:
            # Group by track type (reference, alternate, delta)
            track_type_dir = self.bedgraph_dir / track_type
            track_type_dir.mkdir(exist_ok=True)
            filename = f'{variant!s}_{track_name}.bedgraph'
            return track_type_dir / filename

    def _get_rnaseq_output_path(self, variant, track_name: str, track_type: str, ontology: str) -> Path:
        """Get output path for RNA-seq BED files."""
        if self.organization_mode == OrganizationMode.BY_VARIANT:
            variant_dir = self.bedgraph_dir / f'{variant!s}'
            variant_dir.mkdir(exist_ok=True)

            # Create track_type directory within variant directory
            ontology_dir = variant_dir / ontology
            ontology_dir.mkdir(exist_ok=True)

            # Create junctions directory within variant directory
            rna_dir = ontology_dir / 'RNAseq'
            rna_dir.mkdir(exist_ok=True)

            filename = f'RNA_{track_type}_{track_name}.bedgraph'
            return rna_dir / filename
        track_type_dir = self.bedgraph_dir / track_type
        track_type_dir.mkdir(exist_ok=True)
        filename = f'{variant!s}_RNA_{track_name}.bedgraph'
        return track_type_dir / filename

    def write_bedgraph(
        self,
        variant,
        interval,
        values: np.ndarray,
        track_name: str,
        output_type,
        ontology,
        track_type: str = 'reference',
    ):
        """Write a complete BedGraph file with header and data to organized directory."""
        if output_type == 'rna_seq':
            output_path = self._get_rnaseq_output_path(variant, track_name, track_type, ontology)
        if output_type == 'splice_sites':
            output_path = self._get_splice_sites_output_path(variant, track_name, track_type, ontology)
        if output_type == 'splice_site_usage':
            output_path = self.splice_site_usage_outpath(variant, track_name, track_type, ontology)
        compressed_data = self._compress_values(values, interval.chromosome, interval.start)
        with open(output_path, 'w') as f:
            # Write header
            for line in self._create_header(interval):
                f.write(f'{line}\n')

            # Write track definition
            f.write(self._create_track_line(track_name, track_type))

            # Write data
            for chrom, start, end, value in compressed_data:
                f.write(f'{chrom}\t{start}\t{end}\t{value}\n')

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
        # Define colors based on track type
        if track_type == 'reference':
            color = '128,128,128'  # grey
            alt_color = '128,128,128'  # grey
        elif track_type == 'alternate':
            color = '169,0,17'  # red
            alt_color = '169,0,17'  # red
        elif track_type == 'delta':
            color = '255,192,0'  # yellow
            alt_color = '255,192,0'  # yellow
        else:
            # fallback to default colors
            color = self.positive_color
            alt_color = self.negative_color

        return (
            f'track type=bedGraph '
            f'name=" {track_name}" '
            f'description=" {track_name} ({track_type})" '
            f'visibility=full '
            f'viewLimits={self.default_view_limits} '
            f'autoScale=off '
            f'color={color} '
            f'altColor={alt_color}\n'
        )

    def _compress_values(self, values, chrom: str, start: int) -> list[tuple[str, int, int, float]]:
        """Compress adjacent similar values into regions."""
        if len(values) == 0:
            return []

        data_range = np.nanmax(values) - np.nanmin(values)
        min_delta = data_range * 0.02

        current_region_start = start
        current_region_value = values[0]
        rows_for_file = []
        final_pos = start

        for position, value in enumerate(values[1:], start=start + 1):
            if np.isnan(value):
                continue

            if abs(value - current_region_value) > min_delta:
                rows_for_file.append((chrom, current_region_start, position, current_region_value))
                current_region_start = position
                current_region_value = value

            final_pos = position

        # Add final region
        rows_for_file.append((chrom, current_region_start, final_pos + 1, current_region_value))
        return rows_for_file


class RNAseqTrackGenerator:
    """Generates splice site tracks for variants."""

    def __init__(self, bedgraph_writer: BedGraphWriter):
        self.writer = bedgraph_writer

    def generate_all_tracks(self, variant, interval, scores, ontology):
        """Generate reference, alternate, and delta tracks for all splice sites."""
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
        ref_series = np.array([val[index] for val in scores.reference.rna_seq.values])
        return self.writer.write_bedgraph(variant, interval, ref_series, track_name, 'rna_seq', ontology, 'reference')

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.rna_seq.values]

        if len(variant.reference_bases) != len(variant.alternate_bases):
            alt_series = np.array(self._realign_indel_values(variant, interval, alt_raw))
        else:
            alt_series = np.array(alt_raw)

        return self.writer.write_bedgraph(variant, interval, alt_series, track_name, 'rna_seq', ontology, 'alternate')

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate delta (difference) track."""
        ref_series = np.array([val[index] for val in scores.reference.rna_seq.values])
        alt_raw = [val[index] for val in scores.alternate.rna_seq.values]

        if len(variant.reference_bases) != len(variant.alternate_bases):
            alt_series = np.array(self._realign_indel_values(variant, interval, alt_raw))
        else:
            alt_series = np.array(alt_raw)

        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(variant, interval, delta_series, track_name, 'rna_seq', ontology, 'delta')

    def _realign_indel_values(self, variant, interval, values) -> list[float]:
        """Realign values for indel variants to maintain coordinate mapping."""
        alter_by = abs(len(variant.alternate_bases) - len(variant.reference_bases))
        original_length = len(values)
        midpoint = variant.position - interval.start

        if len(variant.alternate_bases) > len(variant.reference_bases):
            mean_val = np.nanmean(values[midpoint : midpoint + alter_by])
            values = values[:midpoint] + values[midpoint + alter_by :] + [np.nan] * alter_by
            values[midpoint] = mean_val
        else:
            values = np.insert(values, midpoint, [np.nan] * alter_by)

        return values[:original_length]


class SpliceSiteUsageTrackGenerator:
    """Generates splice site tracks for variants."""

    def __init__(self, bedgraph_writer: BedGraphWriter):
        self.writer = bedgraph_writer

    def generate_all_tracks(self, variant, interval, scores, ontology):
        """Generate reference, alternate, and delta tracks for all splice sites."""
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
        ref_series = np.array([val[index] for val in scores.reference.splice_sites.values])
        return self.writer.write_bedgraph(
            variant, interval, ref_series, track_name, 'splice_site_usage', ontology, 'reference'
        )

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]

        if len(variant.reference_bases) != len(variant.alternate_bases):
            alt_series = np.array(self._realign_indel_values(variant, interval, alt_raw))
        else:
            alt_series = np.array(alt_raw)

        return self.writer.write_bedgraph(
            variant, interval, alt_series, track_name, 'splice_site_usage', ontology, 'alternate'
        )

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate delta (difference) track."""
        ref_series = np.array([val[index] for val in scores.reference.splice_sites.values])
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]

        if len(variant.reference_bases) != len(variant.alternate_bases):
            alt_series = np.array(self._realign_indel_values(variant, interval, alt_raw))
        else:
            alt_series = np.array(alt_raw)

        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(
            variant, interval, delta_series, track_name, 'splice_site_usage', ontology, 'delta'
        )

    def _realign_indel_values(self, variant, interval, values) -> list[float]:
        """Realign values for indel variants to maintain coordinate mapping."""
        alter_by = abs(len(variant.alternate_bases) - len(variant.reference_bases))
        original_length = len(values)
        midpoint = variant.position - interval.start

        if len(variant.alternate_bases) > len(variant.reference_bases):
            mean_val = np.nanmean(values[midpoint : midpoint + alter_by])
            values = values[:midpoint] + values[midpoint + alter_by :] + [np.nan] * alter_by
            values[midpoint] = mean_val
        else:
            values = np.insert(values, midpoint, [np.nan] * alter_by)

        return values[:original_length]


class JunctionTrackGenerator:
    """Generates BED tracks for splice junction visualization (sashimi plots)."""

    def __init__(self, bedgraph_writer: BedGraphWriter):
        self.writer = bedgraph_writer

    def generate_junction_bed_tracks(self, variant, interval, scores, ontology):
        """Generate BED files for junction visualization from JunctionData."""
        generated_files = []

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
        output_path = self._get_junction_output_path(variant, track_name, track_type, ontology)

        with open(output_path, 'w') as f:
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
                if score_value < 0.01:
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

                # Color coding: blue for positive strand, red for negative
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

    def _get_junction_output_path(self, variant, track_name: str, track_type: str, ontology: str) -> Path:
        """Get output path for junction BED files."""
        if self.writer.organization_mode == OrganizationMode.BY_VARIANT:
            variant_dir = self.writer.bedgraph_dir / f'{variant!s}'
            variant_dir.mkdir(exist_ok=True)

            # Create track_type directory within variant directory
            ontology_dir = variant_dir / ontology
            ontology_dir.mkdir(exist_ok=True)

            # Create junctions directory within variant directory
            junctions_dir = ontology_dir / 'junctions'
            junctions_dir.mkdir(exist_ok=True)

            filename = f'{track_type}_{track_name}.bed'
            return junctions_dir / filename
        track_type_dir = self.writer.bedgraph_dir / track_type
        track_type_dir.mkdir(exist_ok=True)
        filename = f'{variant!s}_junctions_{track_name}.bed'
        return track_type_dir / filename


class SpliceSiteTrackGenerator:
    """Generates splice site tracks for variants."""

    def __init__(self, bedgraph_writer: BedGraphWriter):
        self.writer = bedgraph_writer

    def generate_all_tracks(self, variant, interval, scores, ontology):
        """Generate reference, alternate, and delta tracks for all splice sites."""
        generated_files = []

        for index in range(scores.reference.splice_sites.num_tracks):
            track_name = scores.reference.splice_sites.names[index]
            track_strand = scores.reference.splice_sites.strands[index]
            full_track_name = f'Splice Sites  {track_strand} {track_name}'

            # Generate all three track types
            ref_file = self._generate_reference_track(variant, interval, scores, index, full_track_name, ontology)
            alt_file = self._generate_alternate_track(variant, interval, scores, index, full_track_name, ontology)
            delta_file = self._generate_delta_track(variant, interval, scores, index, full_track_name, ontology)

            generated_files.extend([ref_file, alt_file, delta_file])

        return generated_files

    def _generate_reference_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate reference sequence track."""
        ref_series = np.array([val[index] for val in scores.reference.splice_sites.values])
        return self.writer.write_bedgraph(
            variant, interval, ref_series, track_name, 'splice_sites', ontology, 'reference'
        )

    def _generate_alternate_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate alternate sequence track with indel realignment."""
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]

        if len(variant.reference_bases) != len(variant.alternate_bases):
            alt_series = np.array(self._realign_indel_values(variant, interval, alt_raw))
        else:
            alt_series = np.array(alt_raw)

        return self.writer.write_bedgraph(
            variant, interval, alt_series, track_name, 'splice_sites', ontology, 'alternate'
        )

    def _generate_delta_track(self, variant, interval, scores, index: int, track_name: str, ontology: str):
        """Generate delta (difference) track."""
        ref_series = np.array([val[index] for val in scores.reference.splice_sites.values])
        alt_raw = [val[index] for val in scores.alternate.splice_sites.values]

        if len(variant.reference_bases) != len(variant.alternate_bases):
            alt_series = np.array(self._realign_indel_values(variant, interval, alt_raw))
        else:
            alt_series = np.array(alt_raw)

        delta_series = np.subtract(alt_series, ref_series)
        return self.writer.write_bedgraph(
            variant, interval, delta_series, track_name, 'splice_sites', ontology, 'delta'
        )

    def _realign_indel_values(self, variant, interval, values) -> list[float]:
        """Realign values for indel variants to maintain coordinate mapping."""
        alter_by = abs(len(variant.alternate_bases) - len(variant.reference_bases))
        original_length = len(values)
        midpoint = variant.position - interval.start

        if len(variant.alternate_bases) > len(variant.reference_bases):
            mean_val = np.nanmean(values[midpoint : midpoint + alter_by])
            values = values[:midpoint] + values[midpoint + alter_by :] + [np.nan] * alter_by
            values[midpoint] = mean_val
        else:
            values = np.insert(values, midpoint, [np.nan] * alter_by)

        return values[:original_length]


import multiprocessing as mp


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


# If generators need to be passed as arguments:
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


def main(var_file: str, output_root: str, ontologies: list[str], organization: str = 'variant'):
    """Main function to process variants and generate BedGraph tracks."""
    model = dna_client.create(API_KEY)
    # as generated from the filter_mt_to_vars_of_interest.py script
    variants = load_variants_table(var_file)
    significant_results: pd.DataFrame | None = None
    THRESHOLD = 0.01
    print(f'Loaded {len(variants)} variants from {var_file!s}.')
    print(f'Using ontology terms: {ontologies!s} with threshold {THRESHOLD}.')
    # Initialize counters for significant results
    # and a dictionary to track the number of significant variants per position
    sig_results_counter = 0
    sig_var_counter_dict: dict[tuple[str, int], int] = {}
    intervals_to_score = []
    vars_to_score = []

    # Parse organization mode
    try:
        org_mode = OrganizationMode(organization)
    except ValueError:
        raise ValueError(f"Invalid organization mode: {organization}. Choose 'variant' or 'track_type'")
    for var in variants:
        interval = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        key = (var.chromosome, var.position)
        variant_scores = model.score_variant(
            interval=interval,
            variant=var,
            variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITE_USAGE']],
        )

        tidied_scores = variant_scorers.tidy_scores(
            [variant_scores],
            match_gene_strand=True,
        )

        filtered_scores = tidied_scores[tidied_scores['raw_score'].values >= THRESHOLD]

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
            intervals_to_score.append(interval)
            vars_to_score.append(var)

    if not intervals_to_score or not vars_to_score:
        print('No variants to score after filtering. Exiting.')
        return None
    for ontology in ontologies:
        # Initialize track generation components
        ontology = [ontology]
        bedgraph_writer = BedGraphWriter(var_file, output_root, ontology[0], org_mode)
        track_generator = SpliceSiteTrackGenerator(bedgraph_writer)
        rna_track_generator = RNAseqTrackGenerator(bedgraph_writer)
        junction_track_generator = JunctionTrackGenerator(bedgraph_writer)  # Add this line
        usage_track_generator = SpliceSiteUsageTrackGenerator(bedgraph_writer)

        variant_predictions = model.predict_variants(
            intervals=intervals_to_score,
            variants=vars_to_score,
            ontology_terms=ontology,
            requested_outputs={
                dna_client.OutputType.SPLICE_SITES,
                dna_client.OutputType.SPLICE_SITE_USAGE,
                dna_client.OutputType.RNA_SEQ,
                dna_client.OutputType.SPLICE_JUNCTIONS,
            },
            max_workers=5,  # Limit to 5 concurrent requests, may need to adjust later based on API quota
        )

        variant_args = list(
            zip(
                variant_predictions,
                vars_to_score,
                intervals_to_score,
                [ontology[0]] * len(variant_predictions),
                strict=False,
            )
        )
        num_processes = min(mp.cpu_count(), len(variant_args))

        with mp.Pool(
            processes=num_processes,
            initializer=init_worker,
            initargs=(track_generator, rna_track_generator, junction_track_generator, usage_track_generator),
        ) as pool:
            results = list(
                tqdm(
                    pool.imap(process_variant_with_globals, variant_args),
                    total=len(variant_args),
                    desc='Processing variants',
                    colour='green',
                )
            )

        all_generated_files = [file for sublist in results for file in sublist]

        print(f'Generated {len(all_generated_files)} track files in {bedgraph_writer.bedgraph_dir}')

    if significant_results is None:
        print('No significant results found.')
        return None

    # Write the significant results to a TSV file via a StringIO buffer
    buffer = io.StringIO()
    significant_results.to_csv(buffer, sep='\t', index=False)
    buffer.seek(0)
    with to_anypath(f'{output_root}.tsv').open('w') as handle:
        handle.write(buffer.read())
    buffer.close()
    print(
        f'{sig_results_counter} out of {len(variants)} variants had significant scores above the threshold {THRESHOLD}.'
    )
    return all_generated_files


if __name__ == '__main__':
    parser = ArgumentParser(description='Generate BedGraph tracks for splice site variants')
    parser.add_argument(
        '--var_file',
        help='Path to variant TSV file',
        default='',
    )
    parser.add_argument(
        '--output_root', help='Root output directory', default=''
    )
    parser.add_argument('--ontology', nargs='+', default=['UBERON:0001134'], help='Ontology terms')  # ,required=True)
    parser.add_argument(
        '--organization',
        choices=['variant', 'track_type'],
        default='variant',
        help="Organization mode: 'variant' groups by variant, 'track_type' groups by reference/alternate/delta",
    )

    args = parser.parse_args()
    main(args.var_file, args.output_root, args.ontology, args.organization)
