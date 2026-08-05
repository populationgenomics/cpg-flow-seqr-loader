#!/usr/bin/env python3
"""
Merge input gVCFs into a single multisample VCF that seqr can load, using bcftools in a Hail Batch job.

seqr expects a *joint-called* VCF (spec:
https://storage.googleapis.com/seqr-reference-data/seqr-vcf-info.pdf). This produces a pragmatic
bcftools-merged equivalent - it is NOT GATK/DRAGEN joint genotyping: calls are stitched together
per-sample rather than re-genotyped across the cohort, QUAL/PL are not recomputed and no VQSR is
applied. seqr's load-time validations still pass; use this for small cohorts (e.g. the synthesized
synthetic-proband trios) where full joint calling is unnecessary.

Each input gVCF is pulled in together with its adjacent .tbi index (assumed to always exist beside
the gVCF), so bcftools merge has the index it needs without an in-job indexing step.

The batch job streams a single pipe, in the bcftools image (no on-disk intermediate):
    1. bcftools merge --gvcf <ref>   - uses INFO/END reference blocks so samples that are hom-ref at
                                       a site read 0/0 (confident) rather than ./. (missing)
    2. bcftools norm -m -any -f <ref> --check-ref x   - split multiallelics + left-align so every row
                                       is one normalised variant with AD length == REF+ALT, dropping
                                       rows whose REF disagrees with the FASTA (non-ACGTN alleles)
    3. bcftools view   - drop the gVCF <NON_REF>/<*> symbolic rows
    4. bcftools norm -d exact   - drop duplicate rows (same POS/REF/ALT) that merge leaves behind for
                                  sites the inputs represent inconsistently (especially indels in
                                  tandem repeats / STRs); seqr's validate_no_duplicate_variants rejects
                                  the whole callset if any variant key repeats. Keeps the first row, so
                                  genotypes carried only on dropped rows are lost - acceptable for a
                                  pragmatic synthetic-trio callset, NOT a substitute for joint genotyping.
    5. tabix index the output

Normalising each input gVCF *before* merge was tried and made this worse (independent left-alignment
fragments tandem-repeat indels into more inconsistent representations, so merge emits more
unreconcilable rows) - dedup after the merge is the reliable approach.

This assumes each input gVCF's header already declares GT/AD/GQ (every GATK/DRAGEN gVCF and the
synthetic gVCF does) and relies on bcftools always emitting ##FILTER=<ID=PASS>; together
those satisfy seqr's header requirements without a separate reheader step.
"""

import argparse
import logging
import sys

from cpg_utils import to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch

GIGABYTE = 1024**3


def main():
    parser = argparse.ArgumentParser(description='Merge gVCFs into a seqr-loadable multisample VCF')
    parser.add_argument('--input', required=True, nargs='+', help='Input gVCF paths (.g.vcf.gz)')
    parser.add_argument('--output', required=True, help='Output multisample VCF path (.vcf.gz)')
    parser.add_argument('--reference', required=True, help='Reference FASTA (.fa/.fasta); expects a sibling .fai')
    args = parser.parse_args()

    if len(args.input) < 2:
        logging.error('At least two input gVCFs are required to merge')
        sys.exit(1)

    batch = get_batch()
    job = batch.new_job(name='Merge gVCFs to VCF')
    job.image(image_path('bcftools'))

    # size storage to hold the localised inputs + reference plus headroom for the streamed output.
    # merge -> norm -> view -> dedup is piped (see below), so there is no multi-GB on-disk
    # intermediate to budget for - only the final output.vcf.gz lands on disk alongside the inputs.
    input_bytes = sum(to_path(p).stat().st_size for p in args.input)
    reference_bytes = to_path(args.reference).stat().st_size
    job.storage(input_bytes + reference_bytes + 10 * GIGABYTE)

    # pull each gVCF in together with its adjacent .tbi (shared root) so merge finds the index
    local_inputs = [
        batch.read_input_group(**{'g.vcf.gz': gvcf, 'g.vcf.gz.tbi': f'{gvcf}.tbi'})['g.vcf.gz']
        for gvcf in args.input
    ]
    inputs_arg = ' '.join(str(gvcf) for gvcf in local_inputs)

    # colocate the .fai beside the .fa (read_input_group localises to a shared root) so `-f` finds it
    reference_ext = args.reference.rsplit('.', 1)[-1]
    reference_group = batch.read_input_group(
        **{reference_ext: args.reference, f'{reference_ext}.fai': f'{args.reference}.fai'},
    )
    reference = reference_group[reference_ext]

    job.declare_resource_group(output={'vcf': '{root}.vcf.gz', 'tbi': '{root}.vcf.gz.tbi'})

    job.command(
        f"""
        set -euo pipefail

        # Stream merge -> norm -> view -> dedup in a single pipe with -Ou (uncompressed BCF) between
        # stages, so no multi-GB intermediate is written to disk. The previous version materialised
        # merged.vcf.gz (~7GB for a WGS trio); on the job disk that filled alongside the localised
        # inputs + reference, silently truncating the intermediate, which norm then failed to read
        # (BGZF read error, exit 9).
        #
        # merge --gvcf uses INFO/END reference blocks so hom-ref samples read 0/0 (confident), not ./.
        # norm -m -any splits multiallelics + left-aligns; --check-ref x drops rows whose REF doesn't
        # match the FASTA (e.g. IUPAC 'Y' where the reference is 'N'), which seqr's non-ACGTN
        # reference-allele validation would otherwise reject. view then drops the gVCF <NON_REF>/<*>
        # symbolic rows so every remaining row is one normalised variant with AD length == REF + ALT.
        #
        # norm -d exact removes duplicate rows the merge leaves behind: bcftools merge stitches samples
        # per-record rather than re-genotyping, so a site that the input gVCFs represent inconsistently
        # (common for indels in tandem repeats / STRs) can land on multiple rows sharing the same
        # POS/REF/ALT once the multiallelics are split. seqr's validate_no_duplicate_variants rejects
        # the whole callset if any variant key repeats. -d exact keeps the FIRST such row and drops the
        # rest, so genotypes carried only on dropped rows are lost - acceptable for a pragmatic
        # synthetic-trio callset (the affected sites are low-confidence repeat-region indels), NOT a
        # substitute for real joint genotyping. It runs last, after the split, so it sees biallelic
        # rows. Normalising each input *before* merge was tried and made this worse (see module docstring).
        bcftools merge --gvcf {reference} -Ou {inputs_arg} \\
            | bcftools norm -m -any -f {reference} --check-ref x -Ou \\
            | bcftools view -e 'ALT="<NON_REF>" || ALT="<*>"' -Ou \\
            | bcftools norm -d exact -Oz -o {job.output.vcf}

        tabix -p vcf {job.output.vcf}
        # {job.output.tbi} is produced by the tabix call above (tabix has no -o); named for hail
        """,
    )

    batch.write_output(job.output.vcf, args.output)
    batch.write_output(job.output.tbi, args.output + '.tbi')

    batch.run(wait=False)


if __name__ == '__main__':
    main()
