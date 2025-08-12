"""
Standalone script to run the AnnotateCohort process in the seqr_loader workflow.

Annotation with gnomAD 4 is a bit tricky. There is a Hail Table containing all the data, so we're using that as a source
Within that HT the data is compressed in a cool way

(I haven't found a README to explain this, but I'm recalling a conversation with a member of the Broad Seqr team)

Each row (keyed on [locus, alleles], i.e. 'row_key') contains a schema, e.g.
`ht[row.key].joint.freq` is the frequency fields in the joint (exome & genome combined) dataset

instead of our current tables where `ref_ht[mt.row_key].gnomad_exomes` is a struct of
{
    AC=int32,
    AN=int32,
    Af=float64,
    ...
}

the gnomad4 data is stored as an array of structs, one for each dataset, e.g.
+------------------------------------------------------------------------------+
| joint.freq                                                                   |
+------------------------------------------------------------------------------+
| array<struct{AC: int32, AF: float64, AN: int32, homozygote_count: int32}>    |
+------------------------------------------------------------------------------+
| [(0,0.00e+00,56642,0),(2,1.65e-05,121094,0),(0,0.00e+00,25546,0),(0,0.00e... |
+------------------------------------------------------------------------------+

Each element of this array is a struct containing the AC, AF, AN, and homozygote_count for a particular gnomAD sub-pop.
The globals of the HT contain a mapping of the dataset name to the index in this array, which is a massive space saving
compared to keeping the full struct for each dataset in the row.

For our purposes here, I'm going to assume we want the joint (exomes and genomes) frequency data, linked to the adjusted
population 'adj' -  genotype calls have been adjusted to account for potential technical artifacts or biases:

1. identify the name(s) of the populations we're interested in:
    - target_pop = "adj"

2. identify the part of the schema we're interested in:
    - "joint.freq"

3. use the globals to find the index of the 'adj' population for joint.freq:
    - "target_index = ht.globals.joint_globals.freq_index_dict[target_pop]"

4. pull out the relevant struct for this population:
    - "ht[row_key].joint.freq[target_index]"

This same process should be repeated for the "adj_XY" population for hemi counts.

The popmax AF/fafmax should be taken from the separate 'max' fields as they are useful aggregated max values:
- joint.grpmax.AF: the max AF seen in any sub-population
- joint.fafmax.faf95_ma: the max filtering allele frequency seen in any sub-population
"""

from argparse import ArgumentParser
from os.path import join

import loguru
from cpg_flow import utils
from cpg_utils import Path, config, hail_batch, to_path

import hail as hl

from cpg_seqr_loader.hail_scripts import variant_id, vep

# adj is the adjusted population, which is the default for gnomAD v4
# all samples across all groups, adjusted for technical artifacts
GNOMAD_TARGET_POP = 'adj'
GNOMAD_XY_TARGET_POP = 'XY_adj'


def load_vqsr(vcf_path: str, ht_path: Path) -> hl.Table:
    """
    Convert VQSR VCF to HT, and checkpoints

    Args:
        vcf_path ():
        ht_path ():
    """
    if utils.can_reuse(ht_path):
        loguru.logger.info(f'Reading VQSR checkpoint from {ht_path}')
        return hl.read_table(str(ht_path))
    loguru.logger.info(f'AS-VQSR: importing annotations from a site-only VCF {vcf_path}')
    vqsr_ht = hl.import_vcf(
        vcf_path,
        reference_genome=hail_batch.genome_build(),
        force_bgz=True,
    ).rows()

    # two comments copied in from the previous implementation, unsure if these are still valid

    # VCF has SB fields as float in header:
    # > ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
    # Even though they are lists of ints, e.g. SB=6,11,2,0
    # Hail would fail to parse it, throwing:
    # > java.lang.NumberFormatException: For input string: "6,11,2,0"
    # To mitigate this, we can drop the SB field before the HT is (lazily) parsed.
    # In order words, dropping it before calling ht.write() makes sure that Hail would
    # never attempt to actually parse it.

    # Dropping also all INFO/AS* annotations as well as InbreedingCoeff, as they are
    # causing problems splitting multiallelics after parsing by Hail, when Hail attempts
    # to subset them by allele index. For example, for these 2 variants:
    # chr1    10145   .       AAC     A,TAC   .       VQSRTrancheINDEL99.50to99.90    AC=0,0;AC_raw=1,1;AS_FS=0,0;AS_FilterStatus=VQSRTrancheINDEL99.50to99.90;AS_MQ=45.5636,46.0067;AS_MQRankSum=0.092,1.34;AS_QD=3.64286,1;AS_QUALapprox=51,20;AS_ReadPosRankSum=0.657,1.128;AS_SB_TABLE=15,15|1,1|1,1;AS_SOR=0.693147,0.693147;AS_VQSLOD=-1.9389;AS_VarDP=14,20;AS_culprit=AS_MQRankSum;DP=1908;FS=0;MQ=45.7857;MQRankSum=1.34;NEGATIVE_TRAIN_SITE;QD=2.08824;QUALapprox=71;ReadPosRankSum=1.128;SB=15,15,2,2;SOR=0.693147;VarDP=34
    # chr1    10146   .       AC      A       .       VQSRTrancheINDEL99.50to99.90    AC=4;AC_raw=5;AS_FS=5.75068;AS_FilterStatus=VQSRTrancheINDEL99.50to99.90;AS_MQ=41.3793;AS_MQRankSum=1.033;AS_QD=12.1209;AS_QUALapprox=1103;AS_ReadPosRankSum=-0.875;AS_SB_TABLE=18,12|28,33;AS_SOR=0.611231;AS_VQSLOD=-2.0660;AS_VarDP=91;AS_culprit=AS_MQRankSum;DP=1727;FS=5.75068;MQ=41.3793;MQRankSum=1.033;NEGATIVE_TRAIN_SITE;QD=12.1209;QUALapprox=1103;ReadPosRankSum=-0.875;SB=18,12,28,33;SOR=0.611231;VarDP=91
    # The first one has 2 alleles, and one AS_FilterStatus value, same as the second
    # one, with one AS_FilterStatus value. So for the first one indexing by allele
    # index would work, but for the second one it would throw an index out of bounds:
    # `HailException: array index out of bounds: index=1, length=1`
    vqsr_ht = vqsr_ht.annotate(info=vqsr_ht.info.drop(*[f for f in vqsr_ht.info if (f.startswith('AS_') or f == 'SB')]))
    return vqsr_ht.checkpoint(str(ht_path), overwrite=True)


def annotate_ourdna(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Annotate the MatrixTable with OurDNA data (exomes & genomes)."""
    if not (
        (ourdna_exome_ht_path := config.config_retrieve(['ourdna', 'exome'], None))
        and (ourdna_genome_ht_path := config.config_retrieve(['ourdna', 'genome'], None))
    ):
        loguru.logger.info('One or both of the OurDNA HTs are not configured, skipping annotation')
        return mt

    ourdna_exomes = hl.read_table(ourdna_exome_ht_path)
    ourdna_genomes = hl.read_table(ourdna_genome_ht_path)

    # the ht.freq block contains a list of Structs, each AC/AF/AN/homozygote_count
    # within this list of structs, the first element is the 'adj' adjusted (qc pass) population, which is what we want
    return mt.annotate_rows(
        ourdna_genomes=ourdna_genomes[mt.row_key].freq[0],
        ourdna_exomes=ourdna_exomes[mt.row_key].freq[0],
    )


def annotate_gnomad4(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    All steps relating to the annotation of gnomAD v4(.1) data

    Args:
        mt (hl.MatrixTable): the input MT

    Returns:
        same MT, with gnomAD 4 annotations placed into the INFO struct as a nested Struct
    """

    gnomad4_ht = hl.read_table(config.reference_path('gnomad_4.1_joint_ht'))

    # the index of the target populations in the joint.freq array
    target_index = hl.eval(gnomad4_ht.globals.joint_globals.freq_index_dict[GNOMAD_TARGET_POP])
    target_xy_index = hl.eval(gnomad4_ht.globals.joint_globals.freq_index_dict[GNOMAD_XY_TARGET_POP])

    return mt.annotate_rows(
        gnomad_genomes=hl.struct(
            # these are taken explicitly from the adj population (across all of gnomAD)
            AC=gnomad4_ht[mt.row_key].joint.freq[target_index].AC,
            AN=gnomad4_ht[mt.row_key].joint.freq[target_index].AN,
            AF=gnomad4_ht[mt.row_key].joint.freq[target_index].AF,
            Hom=gnomad4_ht[mt.row_key].joint.freq[target_index].homozygote_count,
            # a couple of max-value entries
            FAF_AF=gnomad4_ht[mt.row_key].joint.fafmax.faf95_max,
            AF_POPMAX_OR_GLOBAL=gnomad4_ht[mt.row_key].joint.grpmax.AF,
            # not 100% sure about this one... target the `adj_XY` population
            Hemi=gnomad4_ht[mt.row_key].joint.freq[target_xy_index].AC,
        ),
    )


def annotate_with_external_sources(
    mt: hl.MatrixTable,
    vep_ht_path: str,
) -> hl.MatrixTable:
    """
    Annotate MatrixTable with external Tables
    """

    ref_ht_path = config.config_retrieve(['references', 'seqr_combined_reference_data'])
    clinvar_ht_path = config.config_retrieve(['references', 'seqr_clinvar'])

    loguru.logger.info(
        f"""
        Annotating Variant Matrix Table with additional datasets:
        \tVEP annotations from {vep_ht_path}
        \tReference data from {ref_ht_path}
        \tClinVar data from {clinvar_ht_path}
        """,
    )

    clinvar_ht = hl.read_table(clinvar_ht_path)
    ref_ht = hl.read_table(ref_ht_path)
    vep_ht = hl.read_table(vep_ht_path)

    # join all annotation sources into the MatrixTable
    return mt.annotate_rows(
        clinvar_data=clinvar_ht[mt.row_key],
        ref_data=ref_ht[mt.row_key],
        vep=vep_ht[mt.row_key].vep,
    )


def annotate_cohort(
    mt_path: str,
    out_mt_path: str,
    vep_ht_path: str,
    checkpoint_prefix: str,
    vqsr_vcf_path: str | None = None,
) -> None:
    """Convert VCF to matrix table, annotate for Seqr Loader, add VEP, OurDNA, and VQSR annotations."""

    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),
    )

    mt = hl.read_matrix_table(mt_path)
    loguru.logger.info(f'Imported MT from {mt_path} as {mt.n_partitions()} partitions')

    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
    # Map sequencing type to Seqr-style string in global variables
    # https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
    mt = mt.annotate_globals(
        sourceFilePath=mt_path,
        genomeVersion=hail_batch.genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
        sampleType={
            'genome': 'WGS',
            'exome': 'WES',
            'single_cell': 'RNA',
        }.get(sequencing_type, ''),
    )

    if vqsr_vcf_path:
        loguru.logger.info('Adding VQSR annotations into the Matrix Table')
        vqsr_checkpoint = to_path(checkpoint_prefix) / 'vqsr.ht'
        vqsr_ht = load_vqsr(vqsr_vcf_path, vqsr_checkpoint)
        mt = mt.annotate_globals(**vqsr_ht.index_globals())
        mt = mt.annotate_rows(
            info=vqsr_ht[mt.row_key].info,
            filters=vqsr_ht[mt.row_key].filters.filter(lambda val: val != 'PASS'),
        )
        mt = mt.checkpoint(output=join(checkpoint_prefix, 'mt_vep_vqsr.mt'), overwrite=True)

    # join all annotation sources into the MatrixTable
    mt = annotate_with_external_sources(mt, vep_ht_path)
    mt = annotate_ourdna(mt)

    # # annotate all the gnomAD v4 fields in a separate function
    # mt = annotate_gnomad4(mt)

    mt = mt.checkpoint(output=join(checkpoint_prefix, 'mt_plus_ext_tables.mt'), overwrite=True)

    # update AF attributes from observed callset frequencies
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(
        # annotate the info struct with the variant QC data - Hail's variant_qc includes the REF allele, so skip it
        info=mt.info.annotate(
            AC=[mt.variant_qc.AC[1]],
            AF=[mt.variant_qc.AF[1]],
            AN=mt.variant_qc.AN,
        ),
        # taking a single value here for downstream compatibility in Seqr
        AC=mt.info.AC[1],
        AF=mt.info.AF[1],
        AN=mt.info.AN,
    )
    mt = mt.drop('variant_qc')

    loguru.logger.info('Reformatting annotated fields')
    mt = mt.annotate_rows(
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        sortedTranscriptConsequences=vep.get_expr_for_vep_sorted_transcript_consequences_array(mt.vep),
        variantId=variant_id.get_expr_for_variant_id(mt),
        contig=variant_id.get_expr_for_contig(mt.locus),
        pos=mt.locus.position,
        start=mt.locus.position,
        end=mt.locus.position + hl.len(mt.alleles[0]) - 1,
        ref=mt.alleles[0],
        alt=mt.alleles[1],
        xpos=variant_id.get_expr_for_xpos(mt.locus),
        xstart=variant_id.get_expr_for_xpos(mt.locus),
        xstop=variant_id.get_expr_for_xpos(mt.locus) + hl.len(mt.alleles[0]) - 1,
    )

    # this was previously executed in the MtToEs job, as it wasn't possible on QoB
    loguru.logger.info('Adding GRCh37 coords')
    liftover_path = config.reference_path('liftover_38_to_37')
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(liftover_path, rg37)
    mt = mt.annotate_rows(rg37_locus=hl.liftover(mt.locus, 'GRCh37'))

    # only remove InbreedingCoeff if present (post-VQSR)
    if 'InbreedingCoeff' in mt.info:
        mt = mt.annotate_rows(info=mt.info.drop('InbreedingCoeff'))

    loguru.logger.info(
        'Annotating with seqr-loader fields: round 2 (expanding sortedTranscriptConsequences, ref_data, clinvar_data)',
    )
    mt = mt.annotate_rows(
        domains=vep.get_expr_for_vep_protein_domains_set_from_sorted(mt.sortedTranscriptConsequences),
        transcriptConsequenceTerms=vep.get_expr_for_vep_consequence_terms_set(mt.sortedTranscriptConsequences),
        transcriptIds=vep.get_expr_for_vep_transcript_ids_set(mt.sortedTranscriptConsequences),
        mainTranscript=vep.get_expr_for_worst_transcript_consequence_annotations_struct(
            mt.sortedTranscriptConsequences,
        ),
        geneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences),
        codingGeneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences, only_coding_genes=True),
        cadd=mt.ref_data.cadd,
        dbnsfp=mt.ref_data.dbnsfp,
        geno2mp=mt.ref_data.geno2mp,
        gnomad_exomes=mt.ref_data.gnomad_exomes,
        gnomad_exome_coverage=mt.ref_data.gnomad_exome_coverage,
        gnomad_genomes=mt.ref_data.gnomad_genomes,
        gnomad_genome_coverage=mt.ref_data.gnomad_genome_coverage,
        eigen=mt.ref_data.eigen,
        exac=mt.ref_data.exac,
        g1k=mt.ref_data.g1k,
        mpc=mt.ref_data.mpc,
        primate_ai=mt.ref_data.primate_ai,
        splice_ai=mt.ref_data.splice_ai,
        clinvar=hl.struct(
            allele_id=mt.clinvar_data.info.ALLELEID,
            clinical_significance=hl.delimit(mt.clinvar_data.info.CLNSIG),
            gold_stars=mt.clinvar_data.gold_stars,
        ),
    )

    loguru.logger.info('Final Structure:')
    mt.write(out_mt_path, overwrite=True)
    loguru.logger.info(f'Written final matrix table into {out_mt_path}')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input MatrixTable to annotate')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    parser.add_argument('--vep', required=True, help='HT with VEP annotations')
    parser.add_argument('--checkpoint', required=True, help='Checkpoint prefix')
    parser.add_argument('--vqsr', required=False, help='Site-only VQSR VCF')
    args = parser.parse_args()
    annotate_cohort(
        mt_path=args.input,
        out_mt_path=args.output,
        vep_ht_path=args.vep,
        checkpoint_prefix=args.checkpoint,
        vqsr_vcf_path=args.vqsr,
    )


if __name__ == '__main__':
    cli_main()
