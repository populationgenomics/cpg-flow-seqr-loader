"""
Write a seqr-ready VCF for a single Dataset, from the multicohort dense MT and the VQSR sites-only VCF

No further annotation is carried out, seqr calculates everything it needs on import. This:
    - subsets the dense (already split) MT to the Dataset's SGs
    - removes rows without a non-ref call in this Dataset
    - takes the FILTER column from VQSR, and drops all other row-level fields (ID, QUAL, INFO are all '.')
    - keeps only GT, AD, and GQ entries
    - replaces SG IDs with participant external IDs, as required by seqr
"""

import json
from argparse import ArgumentParser

import loguru

import hail as hl

from cpg_utils import config, hail_batch, to_path


def main(mt_path: str, vqsr_vcf: str, id_map_file: str, output: str) -> None:
    """
    Args:
        mt_path (str): path to the multicohort dense MT, with multiallelics split
        vqsr_vcf (str): path to the sites-only VCF with SNP and Indel VQSR applied
        id_map_file (str): JSON file, a mapping of {SG ID: participant external ID} for all SGs to export
        output (str): path to write the VCF to, ending in .vcf.bgz
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['seqr_vcf', 'worker_memory']),
        worker_cores=config.config_retrieve(['seqr_vcf', 'worker_cores']),
        driver_memory=config.config_retrieve(['seqr_vcf', 'driver_memory']),
        driver_cores=config.config_retrieve(['seqr_vcf', 'driver_cores']),
    )

    sg_to_participant: dict[str, str] = json.loads(to_path(id_map_file).read_text())

    mt = hl.read_matrix_table(mt_path)

    if missing := set(sg_to_participant) - set(mt.s.collect()):
        raise ValueError(f'{len(missing)}/{len(sg_to_participant)} requested SG IDs are not in the MT: {missing}')

    loguru.logger.info(f'Subsetting MT to {len(sg_to_participant)} SGs')
    mt = mt.filter_cols(hl.literal(set(sg_to_participant)).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    # only the FILTER column is used, so the INFO fields (some of which Hail can't parse) are never read
    vqsr_ht = hl.import_vcf(vqsr_vcf, reference_genome=hail_batch.genome_build(), force_bgz=True).rows()

    # same treatment as AnnotateCohort: PASS is dropped from the set, and an empty set is written as PASS
    mt = mt.select_rows(filters=vqsr_ht[mt.row_key].filters.filter(lambda val: val != 'PASS'))
    mt = mt.select_entries('GT', 'AD', 'GQ')
    mt = mt.key_cols_by(s=hl.literal(sg_to_participant)[mt.s])

    loguru.logger.info(f'Writing VCF to {output}')
    # carry the VQSR tranche definitions into the header, so every FILTER value is defined
    filter_metadata = hl.get_vcf_metadata(vqsr_vcf)['filter']
    hl.export_vcf(mt, output, tabix=True, metadata={'filter': filter_metadata})


def cli_main():
    parser = ArgumentParser(description='Write a seqr-ready VCF for a single Dataset from the dense MT and VQSR VCF')
    parser.add_argument('--input', required=True, help='Path to the multicohort dense MT')
    parser.add_argument('--vqsr', required=True, help='Path to the sites-only VQSR VCF')
    parser.add_argument('--id_map', required=True, help='JSON mapping of SG ID to participant external ID')
    parser.add_argument('--output', required=True, help='Path to write the VCF to')
    args = parser.parse_args()
    main(mt_path=args.input, vqsr_vcf=args.vqsr, id_map_file=args.id_map, output=args.output)


if __name__ == '__main__':
    cli_main()
