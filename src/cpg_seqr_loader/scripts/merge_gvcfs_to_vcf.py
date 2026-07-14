#!/usr/bin/env python3
"""
Merge input gVCFs into a single multisample VCF that seqr can load, using bcftools in a Hail Batch job.

seqr expects a *joint-called* VCF (spec:
https://storage.googleapis.com/seqr-reference-data/seqr-vcf-info.pdf). This produces a pragmatic
bcftools-merged equivalent - it is NOT GATK/DRAGEN joint genotyping: calls are stitched together
per-sample rather than re-genotyped across the cohort, QUAL/PL are not recomputed and no VQSR is
applied. seqr's load-time validations still pass; use this for small cohorts (e.g. the synthesized
dummy-proband trios) where full joint calling is unnecessary.

Each input gVCF is pulled in together with its adjacent .tbi index (assumed to always exist beside
the gVCF), so bcftools merge has the index it needs without an in-job indexing step.

The batch job runs, in the bcftools image:
    1. bcftools merge --gvcf <ref>   - uses INFO/END reference blocks so samples that are hom-ref at
                                       a site read 0/0 (confident) rather than ./. (missing)
    2. bcftools norm -m -any -f <ref> | drop <NON_REF>/<*>   - split multiallelics + left-align so
                                       every row is one normalised variant with AD length == REF+ALT
    3. tabix index the output

This assumes each input gVCF's header already declares GT/AD/GQ (every GATK/DRAGEN gVCF and the
dummy-synthesis output does) and relies on bcftools always emitting ##FILTER=<ID=PASS>; together
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

    # size storage to hold the localised inputs + reference plus room for intermediates/output
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

        # --gvcf uses INFO/END reference blocks so hom-ref samples read 0/0 (confident), not ./.
        bcftools merge --gvcf {reference} -Oz -o $BATCH_TMPDIR/merged.vcf.gz {inputs_arg}

        # split multiallelics + left-align, then drop the gVCF <NON_REF>/<*> symbolic rows so every
        # remaining row is one normalised variant with AD length == REF + ALT
        bcftools norm -m -any -f {reference} -Ou $BATCH_TMPDIR/merged.vcf.gz \\
            | bcftools view -e 'ALT="<NON_REF>" || ALT="<*>"' -Oz -o {job.output.vcf}

        tabix -p vcf {job.output.vcf}
        # {job.output.tbi} is produced by the tabix call above (tabix has no -o); named for hail
        """,
    )

    batch.write_output(job.output.vcf, args.output)
    batch.write_output(job.output.tbi, args.output + '.tbi')

    batch.run()


if __name__ == '__main__':
    main()
