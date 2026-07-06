"""
Ad-hoc query: at a given locus, extract GT/GQ/DP/AD from the maternal, paternal, and dummy
proband gVCFs and check the dummy against the synthesise_dummy_probands.py expectations:

  expected_GT   from worst-case rule (both het same ALT -> 1/1, one het -> 0/1, ...)
  expected_GQ   min(mother_GQ, father_GQ)                          (see line 392 of the script)
  expected_AB   0.5 for a called 0/1, 1.0 for 1/1                  (fabricated AD, lines 242-253)

Runs entirely in Hail Batch on the analysis-runner image; the only thing that returns to your
laptop is a small TSV in gs://<dataset>-analysis/... (no patient sequence).
"""

import click
from cpg_utils.config import image_path, output_path
from cpg_utils.hail_batch import get_batch


# awk analyser: reads a merged TSV with a leading sample-label column, produces three sections
# (raw records, per-site comparison, mismatches). Kept self-contained so this script has no
# extra image requirements beyond bcftools + awk (both present in the CPG bcftools image).
AWK_ANALYSER = r'''
BEGIN {
    FS = "\t"; OFS = "\t"
    print "=== raw records (sorted by POS) ==="
    print "sample\tCHROM\tPOS\tREF\tALT\tGT\tGQ\tDP_or_MINDP\tAD\tobserved_AB\tis_ref_block"
}
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
function worst_case_gt(m_gt, f_gt) {
    # crude classifier that only handles biallelic ALTs at the same anchor - enough to flag
    # "both parents 0/1 should give 1/1" style mismatches. anything more nuanced falls through
    # to "?" and is skipped in the mismatch section.
    m_has = (m_gt ~ /1/); f_has = (f_gt ~ /1/)
    if (m_has && f_has) return "1/1"
    if (m_has || f_has) return "0/1"
    if (m_gt != "" && f_gt != "") return "0/0"
    return "?"
}
{
    sample = $1; chrom = $2; pos = $3; ref = $4; alt = $5; gt = $6; gq = $7; dp = $8; min_dp = $9; ad = $10
    dp_or_min = (dp == "." || dp == "") ? min_dp : dp
    is_ref = (alt == "<NON_REF>" || alt == "." || alt == "") ? "yes" : "no"
    ab = ab_from_ad(ad)
    ab_s = (ab < 0) ? "." : sprintf("%.3f", ab)
    print sample, chrom, pos, ref, alt, gt, gq, dp_or_min, ad, ab_s, is_ref

    key = chrom "@" pos "@" ref
    if (!(key in seen)) { seen[key] = 1; order[++n_keys] = key }
    got[sample, key] = 1
    gt_[sample, key] = gt
    gq_[sample, key] = gq
    dp_[sample, key] = dp_or_min
    ad_[sample, key] = ad
    ab_[sample, key] = ab
    alt_[sample, key] = alt
    is_ref_[sample, key] = is_ref
}
END {
    print ""
    print "=== per-site comparison (only sites where dummy has a variant record) ==="
    print "CHROM\tPOS\tREF\tmother_ALT/GT/GQ/AB\tfather_ALT/GT/GQ/AB\tdummy_ALT/GT/GQ/AB\texpected_GT\texpected_GQ\texpected_AB\tstatus"

    n_mismatch = 0
    for (i = 1; i <= n_keys; i++) {
        key = order[i]
        if (!((("dummy", key) in got)) || is_ref_["dummy", key] == "yes") continue
        split(key, kp, "@"); chrom = kp[1]; pos = kp[2]; ref = kp[3]

        m_present = ("mother", key) in got; f_present = ("father", key) in got
        m_gt = m_present ? gt_["mother", key] : ""
        f_gt = f_present ? gt_["father", key] : ""
        d_gt = gt_["dummy", key]
        m_gq = m_present ? to_num(gq_["mother", key]) : -1
        f_gq = f_present ? to_num(gq_["father", key]) : -1
        d_gq = to_num(gq_["dummy", key])
        d_ab = ab_["dummy", key]

        exp_gt = worst_case_gt(m_gt, f_gt)
        if (m_gq >= 0 && f_gq >= 0)      exp_gq = (m_gq < f_gq ? m_gq : f_gq)
        else if (m_gq >= 0)              exp_gq = m_gq
        else if (f_gq >= 0)              exp_gq = f_gq
        else                             exp_gq = -1
        exp_ab = expected_ab(exp_gt)

        status = "ok"
        reasons = ""
        if (exp_gt != "?" && d_gt != exp_gt) { status = "MISMATCH"; reasons = reasons "GT(" d_gt "!=" exp_gt ") " }
        if (exp_gq >= 0 && d_gq != exp_gq)   { status = "MISMATCH"; reasons = reasons "GQ(" d_gq "!=" exp_gq ") " }
        if (exp_ab >= 0 && d_ab >= 0 && (d_ab - exp_ab > 0.001 || exp_ab - d_ab > 0.001)) {
            status = "MISMATCH"; reasons = reasons "AB(" sprintf("%.3f", d_ab) "!=" sprintf("%.3f", exp_ab) ") "
        }
        if (!m_present || !f_present) {
            note = ""; if (!m_present) note = note "no_mother_record "; if (!f_present) note = note "no_father_record "
            reasons = reasons note; if (status == "ok") status = "PARENT_MISSING"
        }
        if (status != "ok") n_mismatch++

        m_field = m_present ? gt_["mother", key] "/" (m_gq < 0 ? "." : m_gq) "/" (ab_["mother", key] < 0 ? "." : sprintf("%.3f", ab_["mother", key])) : "-"
        f_field = f_present ? gt_["father", key] "/" (f_gq < 0 ? "." : f_gq) "/" (ab_["father", key] < 0 ? "." : sprintf("%.3f", ab_["father", key])) : "-"
        d_field = d_gt "/" d_gq "/" (d_ab < 0 ? "." : sprintf("%.3f", d_ab))
        m_col = m_present ? alt_["mother", key] "  " m_field : m_field
        f_col = f_present ? alt_["father", key] "  " f_field : f_field
        d_col = alt_["dummy", key] "  " d_field
        exp_ab_s = (exp_ab < 0) ? "." : sprintf("%.3f", exp_ab)
        exp_gq_s = (exp_gq < 0) ? "." : exp_gq

        print chrom, pos, ref, m_col, f_col, d_col, exp_gt, exp_gq_s, exp_ab_s, status "  " reasons
    }

    print ""
    print "=== summary ==="
    print "sites_with_dummy_variant\t" n_keys "  (raw seen keys, including ref blocks)"
    print "mismatches\t" n_mismatch
}
'''


@click.command()
@click.option('--mother', required=True, help='gs:// path to maternal gVCF (.g.vcf.gz)')
@click.option('--father', required=True, help='gs:// path to paternal gVCF (.g.vcf.gz)')
@click.option('--dummy', required=True, help='gs:// path to synthesised dummy proband gVCF')
@click.option('--region', required=True, help='e.g. chr7:72955120-72955135 (window around the anchor)')
@click.option('--out-name', default='site-check', help='basename for the output TSV')
@click.option(
    '--out-category',
    default='analysis',
    show_default=True,
    help='cpg-utils bucket category to write into, e.g. "analysis" or "tmp". '
    '"tmp" lands in gs://cpg-<dataset>-<access>-tmp/... which auto-expires.',
)
def main(mother: str, father: str, dummy: str, region: str, out_name: str, out_category: str) -> None:
    b = get_batch()
    j = b.new_job(f'inspect {region}')
    j.image(image_path('bcftools'))

    inputs = {
        'mother': b.read_input_group(vcf=mother, tbi=mother + '.tbi'),
        'father': b.read_input_group(vcf=father, tbi=father + '.tbi'),
        'dummy': b.read_input_group(vcf=dummy, tbi=dummy + '.tbi'),
    }

    query_fmt = r'%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%MIN_DP\t%AD]\n'
    extract_cmds = [
        f'bcftools query -r {region} -f \'{query_fmt}\' {group.vcf} '
        f'| awk -v s={label} \'BEGIN{{OFS="\\t"}}{{print s,$0}}\''
        for label, group in inputs.items()
    ]

    # write the awk analyser to a file inside the container to avoid shell-escape gymnastics
    j.command(
        f'''set -euo pipefail
cat > /tmp/analyse.awk << 'ANALYSER_EOF'
{AWK_ANALYSER}
ANALYSER_EOF

{{ {" ; ".join(extract_cmds)} ; }} | sort -k2,2 -k3,3n | awk -f /tmp/analyse.awk > {j.out}
'''
    )

    b.write_output(j.out, output_path(f'{out_name}.tsv', out_category))
    b.run(wait=False)


if __name__ == '__main__':
    main()
