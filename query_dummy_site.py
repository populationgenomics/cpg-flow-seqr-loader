"""
Ad-hoc query: at a given locus, extract GT/GQ/DP/AD from the maternal, paternal, and dummy
proband gVCFs and check the dummy against the synthesise_dummy_probands.py expectations:

  script_expected_GQ  = min(mother.GQ, father.GQ) at overlapping records  (line 392)
  matt_expected_GQ    = GQ of the parent(s) whose record ACTUALLY carries the ALT

The gap between these two matters when one parent is on a ref block over the variant site:
the ref block's GQ is the band-wide MINIMUM, not the site-specific GQ, so min()ing it against
the contributing parent's variant-record GQ pulls the dummy's GQ below the contributor's.
Same argument for DP -> MIN_DP.

The MATT_* flags in the output identify sites where the dummy's GQ/DP equal the non-
contributing parent's ref-block value and are strictly lower than the contributor's - direct
evidence for Matt's theory.
"""

import click
from cpg_utils.config import image_path, output_path
from cpg_utils.hail_batch import get_batch


AWK_ANALYSER = r'''
BEGIN { FS = "\t"; OFS = "\t" }
function to_num(x) { return (x == "." || x == "") ? 0 : x + 0 }
function ab_from_ad(ad,    n, a, i, sum, alt_sum) {
    if (ad == "." || ad == "") return -1
    n = split(ad, a, ",")
    sum = 0; alt_sum = 0
    for (i = 1; i <= n; i++) { sum += a[i]; if (i > 1) alt_sum += a[i] }
    return (sum > 0) ? alt_sum / sum : -1
}
function expected_ab(gt) {
    if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") return 0.5
    if (gt == "1/1" || gt == "1|1" || gt == "1")                  return 1.0
    if (gt == "0/0" || gt == "0|0" || gt == "0")                  return 0.0
    return -1
}
function worst_case_gt(m_carries, f_carries, have_m, have_f) {
    if (!have_m || !have_f) return "?"
    if (m_carries && f_carries) return "1/1"
    if (m_carries || f_carries) return "0/1"
    return "0/0"
}
# Find the index of `sample`'s record overlapping target_pos. Prefer a variant record over a
# ref block when both overlap (shouldn't happen in a well-formed gVCF, but defensive).
function find_overlap(sample, target_pos,    i, n, best) {
    best = 0; n = count[sample] + 0
    for (i = 1; i <= n; i++) {
        if (pos_[sample, i] <= target_pos && end_[sample, i] >= target_pos) {
            if (best == 0 || (is_ref_[sample, best] == "yes" && is_ref_[sample, i] == "no")) best = i
        }
    }
    return best
}

{
    sample = $1; chrom = $2; pos = $3 + 0; info_end = $4; ref = $5; alt = $6
    gt = $7; gq = $8; dp = $9; min_dp = $10; ad = $11
    is_ref = (alt == "<NON_REF>" || alt == "." || alt == "") ? "yes" : "no"
    end = (info_end != "." && info_end != "") ? info_end + 0 : pos + length(ref) - 1
    dp_eff = (dp == "." || dp == "") ? min_dp : dp
    ab = ab_from_ad(ad)
    ab_s = (ab < 0) ? "." : sprintf("%.3f", ab)

    n = ++count[sample]
    chrom_[sample, n] = chrom; pos_[sample, n] = pos; end_[sample, n] = end
    ref_[sample, n] = ref; alt_[sample, n] = alt; gt_[sample, n] = gt
    gq_[sample, n] = gq; dp_[sample, n] = dp_eff; ad_[sample, n] = ad
    is_ref_[sample, n] = is_ref; ab_[sample, n] = ab; ab_s_[sample, n] = ab_s

    raw_n++
    raw_rows[raw_n] = sample "\t" chrom "\t" pos "\t" end "\t" ref "\t" alt "\t" gt "\t" gq "\t" dp_eff "\t" ad "\t" ab_s "\t" is_ref
}

END {
    print "=== raw records (per sample, sorted by (sample, POS)) ==="
    print "sample\tCHROM\tPOS\tEND\tREF\tALT\tGT\tGQ\tDP_or_MINDP\tAD\tobs_AB\tref_block"
    for (i = 1; i <= raw_n; i++) print raw_rows[i]

    print ""
    print "=== per-dummy-variant comparison ==="
    print "CHROM\tPOS\tREF\tALT\tdummy(GT/GQ/DP/AB)\tmother_overlap(kind:GT/GQ/DP)\tfather_overlap(kind:GT/GQ/DP)\tscript_exp(GT/GQ/AB)\tmatt_exp(GQ/DP)\tflags"

    n_dummy_variants = 0; n_gt_wrong = 0; n_gq_ne_script = 0; n_ab_ne_script = 0
    n_matt_gq_confirmed = 0; n_matt_dp_confirmed = 0
    n_check = count["dummy"] + 0
    for (i = 1; i <= n_check; i++) {
        if (is_ref_["dummy", i] == "yes") continue
        n_dummy_variants++

        d_chrom = chrom_["dummy", i]; d_pos = pos_["dummy", i]; d_ref = ref_["dummy", i]
        d_alt = alt_["dummy", i]; d_gt = gt_["dummy", i]
        d_gq = to_num(gq_["dummy", i]); d_dp = to_num(dp_["dummy", i]); d_ab = ab_["dummy", i]

        m_idx = find_overlap("mother", d_pos); f_idx = find_overlap("father", d_pos)

        have_m = (m_idx > 0); have_f = (f_idx > 0)
        m_is_var = have_m && (is_ref_["mother", m_idx] == "no")
        f_is_var = have_f && (is_ref_["father", f_idx] == "no")
        # a parent "carries" the ALT if its overlapping record is a variant record with a non-zero GT allele
        m_carries = m_is_var && (gt_["mother", m_idx] ~ /1/)
        f_carries = f_is_var && (gt_["father", f_idx] ~ /1/)

        m_gq = have_m ? to_num(gq_["mother", m_idx]) : -1
        f_gq = have_f ? to_num(gq_["father", f_idx]) : -1
        m_dp = have_m ? to_num(dp_["mother", m_idx]) : -1
        f_dp = have_f ? to_num(dp_["father", f_idx]) : -1

        exp_gt = worst_case_gt(m_carries, f_carries, have_m, have_f)

        # script-expected: min across both overlaps (mimics emit_variant on lines 392-393)
        if (m_gq >= 0 && f_gq >= 0) script_gq = (m_gq < f_gq ? m_gq : f_gq)
        else if (m_gq >= 0)         script_gq = m_gq
        else if (f_gq >= 0)         script_gq = f_gq
        else                        script_gq = -1

        # matt-expected: min over just the contributing (variant-record) parents
        matt_gq = -1; matt_dp = -1
        if (m_carries) { matt_gq = m_gq; matt_dp = m_dp }
        if (f_carries) {
            if (matt_gq < 0 || f_gq < matt_gq) matt_gq = f_gq
            if (matt_dp < 0 || f_dp < matt_dp) matt_dp = f_dp
        }

        script_ab = expected_ab(exp_gt)

        flags = ""
        if (exp_gt != "?" && d_gt != exp_gt) { flags = flags " GT_ne_expected(" d_gt "!=" exp_gt ")"; n_gt_wrong++ }
        if (script_gq >= 0 && d_gq != script_gq) { flags = flags " GQ_ne_script(" d_gq "!=" script_gq ")"; n_gq_ne_script++ }
        if (script_ab >= 0 && d_ab >= 0 && (d_ab - script_ab > 0.001 || script_ab - d_ab > 0.001)) {
            flags = flags " AB_ne_script(" sprintf("%.3f", d_ab) "!=" sprintf("%.3f", script_ab) ")"
            n_ab_ne_script++
        }

        # ---- Matt's theory test: one contributor + one ref block, dummy's GQ/DP match the ref block's ----
        n_contrib = m_carries + f_carries
        if (n_contrib == 1 && have_m && have_f && (m_is_var != f_is_var)) {
            if (m_carries) { c_gq = m_gq; c_dp = m_dp; o_gq = f_gq; o_dp = f_dp; other = "father_refblk" }
            else           { c_gq = f_gq; c_dp = f_dp; o_gq = m_gq; o_dp = m_dp; other = "mother_refblk" }
            if (c_gq > o_gq && d_gq == o_gq) {
                flags = flags " MATT_GQ_CONFIRMED(" other "_GQ=" o_gq "<contrib_GQ=" c_gq ")"; n_matt_gq_confirmed++
            }
            if (c_dp > o_dp && d_dp == o_dp) {
                flags = flags " MATT_DP_CONFIRMED(" other "_DP=" o_dp "<contrib_DP=" c_dp ")"; n_matt_dp_confirmed++
            }
        }
        if (flags == "") flags = " ok"

        m_kind = have_m ? (m_is_var ? "VAR" : "REFBLK") : "-"
        f_kind = have_f ? (f_is_var ? "VAR" : "REFBLK") : "-"
        m_desc = have_m ? m_kind ":" gt_["mother", m_idx] "/" m_gq "/" m_dp : "-"
        f_desc = have_f ? f_kind ":" gt_["father", f_idx] "/" f_gq "/" f_dp : "-"
        d_desc = d_gt "/" d_gq "/" d_dp "/" (d_ab < 0 ? "." : sprintf("%.3f", d_ab))
        script_exp = exp_gt "/" (script_gq < 0 ? "." : script_gq) "/" (script_ab < 0 ? "." : sprintf("%.3f", script_ab))
        matt_exp = (matt_gq < 0 ? "." : matt_gq) "/" (matt_dp < 0 ? "." : matt_dp)

        print d_chrom, d_pos, d_ref, d_alt, d_desc, m_desc, f_desc, script_exp, matt_exp, substr(flags, 2)
    }

    print ""
    print "=== summary ==="
    print "dummy_variant_records_checked\t" n_dummy_variants
    print "  GT_mismatches\t\t\t" n_gt_wrong
    print "  GQ_mismatches_vs_script\t" n_gq_ne_script
    print "  AB_mismatches_vs_script\t" n_ab_ne_script
    print "  MATT_GQ_theory_confirmed\t" n_matt_gq_confirmed "  (dummy_GQ = refblk_GQ < contrib_GQ)"
    print "  MATT_DP_theory_confirmed\t" n_matt_dp_confirmed "  (dummy_DP = refblk_MIN_DP < contrib_DP)"
}
'''


def _expand_variant_to_region(variant: str, window: int) -> str:
    """Convert a CPG-style variant key 'chr7:72955124:T:TTTC' into a bcftools region string.

    The window pads either side of the anchor POS to catch nearby records (e.g. a spanning ref
    block whose POS is earlier, or a differently-anchored indel from the other parent).
    """
    parts = variant.split(':')
    if len(parts) < 2:
        raise click.BadParameter(f'--variant must be at least chrom:pos, got {variant!r}')
    chrom, pos_s = parts[0], parts[1]
    pos = int(pos_s)
    return f'{chrom}:{max(1, pos - window)}-{pos + window}'


@click.command()
@click.option('--mother', required=True, help='gs:// path to maternal gVCF (.g.vcf.gz)')
@click.option('--father', required=True, help='gs:// path to paternal gVCF (.g.vcf.gz)')
@click.option('--dummy', required=True, help='gs:// path to synthesised dummy proband gVCF')
@click.option(
    '--variant',
    'variants',
    multiple=True,
    help='Variant key chr:pos[:ref[:alt]] - auto-expanded to a window around the anchor. '
    'Repeat the flag for multiple sites, e.g. --variant chr7:72955124:T:TTTC --variant chr1:12345:A:G',
)
@click.option(
    '--region',
    'regions',
    multiple=True,
    help='Explicit region chr:start-end. Repeatable. Use this or --variant (or both).',
)
@click.option(
    '--window',
    default=15,
    show_default=True,
    help='Padding (bp) added either side of each --variant anchor when expanding to a region.',
)
@click.option('--out-name', default='site-check', help='basename for the output TSV')
@click.option(
    '--out-category',
    default='analysis',
    show_default=True,
    help='cpg-utils bucket category, e.g. "analysis" or "tmp" (tmp bucket auto-expires after ~8 days).',
)
def main(
    mother: str,
    father: str,
    dummy: str,
    variants: tuple[str, ...],
    regions: tuple[str, ...],
    window: int,
    out_name: str,
    out_category: str,
) -> None:
    all_regions = list(regions) + [_expand_variant_to_region(v, window) for v in variants]
    if not all_regions:
        raise click.UsageError('provide at least one --variant or --region')
    # bcftools -r accepts a comma-separated list; one call per file is much cheaper than one per site
    region_arg = ','.join(all_regions)

    b = get_batch()
    j = b.new_job(f'inspect {len(all_regions)} site(s)')
    j.image(image_path('bcftools'))

    inputs = {
        'mother': b.read_input_group(vcf=mother, tbi=mother + '.tbi'),
        'father': b.read_input_group(vcf=father, tbi=father + '.tbi'),
        'dummy': b.read_input_group(vcf=dummy, tbi=dummy + '.tbi'),
    }

    query_fmt = r'%CHROM\t%POS\t%INFO/END\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%MIN_DP\t%AD]\n'
    extract_cmds = [
        f'bcftools query -r {region_arg} -f \'{query_fmt}\' {group.vcf} '
        f'| awk -v s={label} \'BEGIN{{OFS="\\t"}}{{print s,$0}}\''
        for label, group in inputs.items()
    ]

    j.command(
        f'''set -euo pipefail
cat > /tmp/analyse.awk << 'ANALYSER_EOF'
{AWK_ANALYSER}
ANALYSER_EOF

{{ {" ; ".join(extract_cmds)} ; }} | sort -k1,1 -k3,3n | awk -f /tmp/analyse.awk > {j.out}
'''
    )

    b.write_output(j.out, output_path(f'{out_name}.tsv', out_category))
    b.run(wait=False)


if __name__ == '__main__':
    main()
