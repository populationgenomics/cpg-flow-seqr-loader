"""
Script to parse JSON VEP results.

Decoding a JSON string is an expensive operation, so we provide a fixed schema.
"""

import argparse

import hail as hl

from cpg_utils import config, hail_batch


def vep_json_to_ht(vep_result_paths: list[str], out_path: str):
    """
    Parse results from VEP with annotations formatted in JSON, and write into a Hail Table.
    """

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),
    )

    json_schema = hl.dtype(
        """
        struct{
            minimised:int32,
            assembly_name:str,
            allele_string:str,
            ancestral:str,
            colocated_variants:array<struct{
                allele_string:str,
                clin_sig:array<str>,
                clin_sig_allele:str,
                end:int32,
                id:str,
                minimised:int32,
                minor_allele:str,
                minor_allele_freq:float64,
                phenotype_or_disease:int32,
                pubmed:array<int32>,
                seq_region_name:str,
                somatic:int32,
                start:int32,
                strand:int32
            }>,
            context:str,
            end:int32,
            id:str,
            input:str,
            intergenic_consequences:array<struct{
                allele_num:int32,
                consequence_terms:array<str>,
                impact:str,minimised:int32,
                variant_allele:str
            }>,
            most_severe_consequence:str,
            motif_feature_consequences:array<struct{
                allele_num:int32,
                consequence_terms:array<str>,
                high_inf_pos:str,
                impact:str,
                minimised:int32,
                motif_feature_id:str,
                motif_name:str,
                motif_pos:int32,
                motif_score_change:float64,
                strand:int32,
                transcription_factors:array<str>,
                variant_allele:str
            }>,
            regulatory_feature_consequences:array<struct{
                allele_num:int32,
                biotype:str,
                consequence_terms:array<str>,
                impact:str,
                minimised:int32,
                regulatory_feature_id:str,
                variant_allele:str
            }>,
            seq_region_name:str,
            start:int32,
            strand:int32,
            transcript_consequences:array<struct{
                allele_num:int32,
                amino_acids:str,
                appris:str,
                biotype:str,
                canonical:int32,
                mane_select:str,
                mane_plus_clinical:str,
                ccds:str,
                cdna_start:int32,
                cdna_end:int32,
                cds_end:int32,
                cds_start:int32,
                codons:str,
                consequence_terms:array<str>,
                distance:int32,
                domains:array<struct{
                    db:str,
                    name:str
                }>,
                exon:str,
                gene_id:str,
                gene_pheno:int32,
                gene_symbol:str,
                gene_symbol_source:str,
                hgnc_id:str,
                hgvsc:str,
                hgvsp:str,
                hgvs_offset:int32,
                impact:str,
                intron:str,
                lof:str,
                lof_flags:str,
                lof_filter:str,
                lof_info:str,
                existing_inframe_oorfs:int32,
                existing_outofframe_oorfs:int32,
                existing_uorfs:int32,
                5utr_consequence:str,
                5utr_annotation:dict<
                    str,
                    struct{
                        type:str,
                        KozakContext:str,
                        KozakStrength:str,
                        DistanceToCDS:str,
                        CapDistanceToStart:str,
                        DistanceToStop:str,
                        Evidence:str,
                        AltStop:str,
                        AltStopDistanceToCDS:str,
                        FrameWithCDS:str,
                        StartDistanceToCDS:str,
                        newSTOPDistanceToCDS:str,
                        alt_type:str,
                        alt_type_length:str,
                        ref_StartDistanceToCDS:str,
                        ref_type:str,
                        ref_type_length:str
                    }
                >,
                minimised:int32,
                mirna:array<str>,
                polyphen_prediction:str,
                polyphen_score:float64,
                protein_end:int32,
                protein_start:int32,
                protein_id:str,
                sift_prediction:str,
                sift_score:float64,
                strand:int32,
                swissprot:array<str>,
                transcript_id:str,
                trembl:array<str>,
                tsl:int32,
                uniparc:array<str>,
                uniprot_isoform:array<str>,
                variant_allele:str,
                am_class:str,
                am_pathogenicity:float64,
                source:str,
                flags:array<str>
            }>,
            variant_class:str
        }
    """,
    )

    ht = hl.import_table(paths=vep_result_paths, no_header=True, types={'f0': json_schema})
    ht = ht.transmute(vep=ht.f0)
    split_vep_input = ht.vep.input.split('\t')
    ht = ht.annotate(
        locus=hl.locus(ht.vep.seq_region_name, hl.parse_int(split_vep_input[1])),
        alleles=[split_vep_input[3], split_vep_input[4]],
    )
    ht = ht.key_by(ht.locus, ht.alleles)
    ht.write(out_path, overwrite=True)


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='VEP results JSON', required=True, nargs='+')
    parser.add_argument('--output', help='Output Hail table', required=True)
    args = parser.parse_args()

    vep_json_to_ht(vep_result_paths=args.input, out_path=args.output)


if __name__ == '__main__':
    cli_main()
