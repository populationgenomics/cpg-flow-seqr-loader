# noqa: D100

import logging
from typing import Any, Optional, Union

import hail as hl

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def to_phred(linear_expr: hl.expr.NumericExpression) -> hl.expr.Float64Expression:
    """
    Compute the phred-scaled value of the linear-scale input.

    :param linear_expr: input
    :return: Phred-scaled value
    """
    return -10 * hl.log10(linear_expr)


def get_lowqual_expr(
    alleles: hl.expr.ArrayExpression,
    qual_approx_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
    snv_phred_threshold: int = 30,
    snv_phred_het_prior: int = 30,  # 1/1000
    indel_phred_threshold: int = 30,
    indel_phred_het_prior: int = 39,  # 1/8,000
) -> Union[hl.expr.BooleanExpression, hl.expr.ArrayExpression]:
    """
    Compute lowqual threshold expression for either split or unsplit alleles based on QUALapprox or AS_QUALapprox.

    .. note::

        When running This lowqual annotation using QUALapprox, it differs from the GATK LowQual filter.
        This is because GATK computes this annotation at the site level, which uses the least stringent prior for mixed sites.
        When run using AS_QUALapprox, this implementation can thus be more stringent for certain alleles at mixed sites.

    :param alleles: Array of alleles
    :param qual_approx_expr: QUALapprox or AS_QUALapprox
    :param snv_phred_threshold: Phred-scaled SNV "emission" threshold (similar to GATK emission threshold)
    :param snv_phred_het_prior: Phred-scaled SNV heterozygosity prior (30 = 1/1000 bases, GATK default)
    :param indel_phred_threshold: Phred-scaled indel "emission" threshold (similar to GATK emission threshold)
    :param indel_phred_het_prior: Phred-scaled indel heterozygosity prior (30 = 1/1000 bases, GATK default)
    :return: lowqual expression (BooleanExpression if `qual_approx_expr`is Numeric, Array[BooleanExpression] if `qual_approx_expr` is ArrayNumeric)
    """
    min_snv_qual = snv_phred_threshold + snv_phred_het_prior
    min_indel_qual = indel_phred_threshold + indel_phred_het_prior
    min_mixed_qual = max(min_snv_qual, min_indel_qual)

    if isinstance(qual_approx_expr, hl.expr.ArrayNumericExpression):
        return hl.range(1, hl.len(alleles)).map(
            lambda ai: hl.if_else(
                hl.is_snp(alleles[0], alleles[ai]),
                qual_approx_expr[ai - 1] < min_snv_qual,
                qual_approx_expr[ai - 1] < min_indel_qual,
            )
        )
    else:
        return (
            hl.case()
            .when(
                hl.range(1, hl.len(alleles)).all(lambda ai: hl.is_snp(alleles[0], alleles[ai])),
                qual_approx_expr < min_snv_qual,
            )
            .when(
                hl.range(1, hl.len(alleles)).all(lambda ai: hl.is_indel(alleles[0], alleles[ai])),
                qual_approx_expr < min_indel_qual,
            )
            .default(qual_approx_expr < min_mixed_qual)
        )


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.

    Defaults correspond to gnomAD values.
    """
    return (
        (gq_expr >= adj_gq)
        & hl.if_else(gt_expr.is_haploid(), dp_expr >= haploid_adj_dp, dp_expr >= adj_dp)
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            .default((ad_expr[gt_expr[0]] / dp_expr >= adj_ab) & (ad_expr[gt_expr[1]] / dp_expr >= adj_ab))
        )
    )


def annotation_type_is_numeric(t: Any) -> bool:
    """
    Given an annotation type, return whether it is a numerical type or not.

    :param t: Type to test
    :return: If the input type is numeric
    """
    return t in (hl.tint32, hl.tint64, hl.tfloat32, hl.tfloat64)


def annotation_type_in_vcf_info(t: Any) -> bool:
    """
    Given an annotation type, returns whether that type can be natively exported to a VCF INFO field.

    .. note::

        Types that aren't natively exportable to VCF will be converted to String on export.

    :param t: Type to test
    :return: If the input type can be exported to VCF
    """
    return (
        annotation_type_is_numeric(t)
        or t in (hl.tstr, hl.tbool)
        or (isinstance(t, (hl.tarray, hl.tset)) and annotation_type_in_vcf_info(t.element_type))
    )


def fs_from_sb(
    sb: Union[hl.expr.ArrayNumericExpression, hl.expr.ArrayExpression],
    normalize: bool = True,
    min_cell_count: int = 200,
    min_count: int = 4,
    min_p_value: float = 1e-320,
) -> hl.expr.Int64Expression:
    """
    Compute `FS` (Fisher strand balance) annotation from  the `SB` (strand balance table) field.

    `FS` is the phred-scaled value of the double-sided Fisher exact test on strand balance.

    Using default values will have the same behavior as the GATK implementation, that is:
    - If sum(counts) > 2*`min_cell_count` (default to GATK value of 200), they are normalized
    - If sum(counts) < `min_count` (default to GATK value of 4), returns missing
    - Any p-value < `min_p_value` (default to GATK value of 1e-320) is truncated to that value

    In addition to the default GATK behavior, setting `normalize` to `False` will perform a chi-squared test
    for large counts (> `min_cell_count`) instead of normalizing the cell values.

    .. note::

        This function can either take
        - an array of length four containing the forward and reverse strands' counts of ref and alt alleles: [ref fwd, ref rev, alt fwd, alt rev]
        - a two dimensional array with arrays of length two, containing the counts: [[ref fwd, ref rev], [alt fwd, alt rev]]

    GATK code here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/FisherStrand.java

    :param sb: Count of ref/alt reads on each strand
    :param normalize: Whether to normalize counts is sum(counts) > min_cell_count (normalize=True), or use a chi sq instead of FET (normalize=False)
    :param min_cell_count: Maximum count for performing a FET
    :param min_count: Minimum total count to output FS (otherwise null it output)
    :return: FS value
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(lambda x: hl.flatten(x), sb)

    sb_sum = hl.bind(lambda x: hl.sum(x), sb)

    # Normalize table if counts get too large
    if normalize:
        fs_expr = hl.bind(
            lambda sb, sb_sum: hl.if_else(
                sb_sum <= 2 * min_cell_count,
                sb,
                sb.map(lambda x: hl.int(x / (sb_sum / min_cell_count))),
            ),
            sb,
            sb_sum,
        )

        # FET
        fs_expr = to_phred(
            hl.max(
                hl.fisher_exact_test(fs_expr[0], fs_expr[1], fs_expr[2], fs_expr[3]).p_value,
                min_p_value,
            )
        )
    else:
        fs_expr = to_phred(
            hl.max(
                hl.contingency_table_test(sb[0], sb[1], sb[2], sb[3], min_cell_count=min_cell_count).p_value,
                min_p_value,
            )
        )

    # Return null if counts <= `min_count`
    return hl.or_missing(
        sb_sum > min_count,
        hl.max(0, fs_expr),  # Needed to avoid -0.0 values
    )


def sor_from_sb(sb: Union[hl.expr.ArrayNumericExpression, hl.expr.ArrayExpression]) -> hl.expr.Float64Expression:
    """
    Compute `SOR` (Symmetric Odds Ratio test) annotation from  the `SB` (strand balance table) field.

    .. note::

        This function can either take
        - an array of length four containing the forward and reverse strands' counts of ref and alt alleles: [ref fwd, ref rev, alt fwd, alt rev]
        - a two dimensional array with arrays of length two, containing the counts: [[ref fwd, ref rev], [alt fwd, alt rev]]

    GATK code here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandOddsRatio.java

    :param sb: Count of ref/alt reads on each strand
    :return: SOR value
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(lambda x: hl.flatten(x), sb)

    sb = sb.map(lambda x: hl.float64(x) + 1)

    ref_fw = sb[0]
    ref_rv = sb[1]
    alt_fw = sb[2]
    alt_rv = sb[3]
    symmetrical_ratio = ((ref_fw * alt_rv) / (alt_fw * ref_rv)) + ((alt_fw * ref_rv) / (ref_fw * alt_rv))
    ref_ratio = hl.min(ref_rv, ref_fw) / hl.max(ref_rv, ref_fw)
    alt_ratio = hl.min(alt_fw, alt_rv) / hl.max(alt_fw, alt_rv)
    sor = hl.log(symmetrical_ratio) + hl.log(ref_ratio) - hl.log(alt_ratio)

    return sor


def pab_max_expr(
    gt_expr: hl.expr.CallExpression,
    ad_expr: hl.expr.ArrayExpression,
    la_expr: Optional[hl.expr.ArrayExpression] = None,
    n_alleles_expr: Optional[hl.expr.Int32Expression] = None,
) -> hl.expr.ArrayExpression:
    """
    Compute the maximum p-value of the binomial test for the alternate allele balance (PAB) for each allele.

    .. note::

        This function can take a `gt_expr` and `ad_expr` that use local or global
        alleles. If they use local alleles, `la_expr` and `n_alleles_expr` should be
        provided to transform `gt_expr` and `ad_expr` to global alleles.

    :param gt_expr: Genotype call expression.
    :param ad_expr: Allele depth expression.
    :param la_expr: Allele local index expression. When provided `gt_expr` and
        `ad_expr` are transformed from using local alleles to global alleles using
        `la_expr`.
    :param n_alleles_expr: Number of alleles expression. Required when 'la_expr' is
        provided.
    :return: Array expression of maximum p-values.
    """
    if la_expr is not None:
        if n_alleles_expr is None:
            raise ValueError('Must provide `n_alleles_expr` if `la_expr` is provided!')

        ad_expr = hl.vds.local_to_global(ad_expr, la_expr, n_alleles_expr, fill_value=0, number='R')
        gt_expr = hl.vds.lgt_to_gt(gt_expr, la_expr)

    expr = hl.agg.array_agg(
        lambda x: hl.agg.filter(
            gt_expr.is_het(),
            hl.agg.max(hl.binom_test(x, hl.sum(ad_expr), 0.5, 'two-sided')),
        ),
        ad_expr[1:],  # Skip ref allele
    )

    return expr
