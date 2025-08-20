import logging

import hail as hl

POP_NAMES = {
    "afr": "African/African-American",
    "ami": "Amish",
    "amr": "Latino",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "eur": "European",
    "fin": "Finnish",
    # NOTE: mde is kept for historical purposes, in gnomAD v3.1 mid was used instead
    "mde": "Middle Eastern",
    "mid": "Middle Eastern",
    "nfe": "Non-Finnish European",
    "oth": "Other",
    "sas": "South Asian",
    "uniform": "Uniform",
    "sas_non_consang": "South Asian (F < 0.05)",
    "consanguineous": "South Asian (F > 0.05)",
    "exac": "ExAC",
    "bgr": "Bulgarian (Eastern European)",
    "est": "Estonian",
    "gbr": "British",
    "nwe": "North-Western European",
    "seu": "Southern European",
    "swe": "Swedish",
    "kor": "Korean",
    "sgp": "Singaporean",
    "jpn": "Japanese",
    "oea": "Other East Asian",
    "oeu": "Other European",
    "onf": "Other Non-Finnish European",
    "unk": "Unknown",
}
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SORT_ORDER = [
    "subset",
    "downsampling",
    "popmax",
    "grpmax",
    "pop",
    "gen_anc",
    "subpop",
    "sex",
    "group",
]
"""
Order to sort subgroupings during VCF export.
Ensures that INFO labels in VCF are in desired order (e.g., raw_AC_afr_female).
"""

GROUPS = ["adj", "raw"]
"""
Group names used to generate labels for high quality genotypes and all raw genotypes. Used in VCF export.
"""

HISTS = ["gq_hist_alt", "gq_hist_all", "dp_hist_alt", "dp_hist_all", "ab_hist_alt"]
"""
Quality histograms used in VCF export.
"""

FAF_POPS = ["afr", "amr", "eas", "nfe", "sas"]
"""
Global populations that are included in filtering allele frequency (faf) calculations. Used in VCF export.
"""

SEXES = ["XX", "XY"]
"""
Sample sexes used in VCF export.

Used to stratify frequency annotations (AC, AN, AF) for each sex.
Note that sample sexes in gnomAD v3 and earlier were 'male' and 'female'.
"""


INFO_VCF_AS_PIPE_DELIMITED_FIELDS = [
    "AS_QUALapprox",
    "AS_VarDP",
    "AS_MQ_DP",
    "AS_RAW_MQ",
    "AS_SB_TABLE",
]

"""
Dictionary used during VCF export to export MatrixTable entries.
"""


def adjust_vcf_incompatible_types(
    ht: hl.Table,
    pipe_delimited_annotations: list[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.Table:
    """
    Create a Table ready for vcf export.

    In particular, the following conversions are done:
        - All int64 are coerced to int32
        - Fields specified by `pipe_delimited_annotations` are converted from arrays to pipe-delimited strings

    :param ht: Input Table.
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct).
    :return: Table ready for VCF export.
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    ht = ht.key_by("locus", "alleles")

    info_type_convert_expr = {}
    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            logger.warning(
                "Coercing field info.%s from int64 to int32 for VCF output. Value will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update({f: hl.int32(hl.min(2**31 - 1, ht.info[f]))})
        elif ft == hl.dtype("array<int64>"):
            logger.warning(
                "Coercing field info.%s from array<int64> to array<int32> for VCF"
                " output. Array values will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update({f: ht.info[f].map(lambda x: hl.int32(hl.min(2**31 - 1, x)))})

    ht = ht.annotate(info=ht.info.annotate(**info_type_convert_expr))

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, it
    # is handled below
    for f in pipe_delimited_annotations:
        if f in ht.info and f != "AS_SB_TABLE":
            info_expr[f] = "|" + get_pipe_expr(ht.info[f])

    # Flatten SB if it is an array of arrays
    if "SB" in ht.info and not isinstance(ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr["SB"] = ht.info.SB[0].extend(ht.info.SB[1])

    if "AS_SB_TABLE" in ht.info:
        info_expr["AS_SB_TABLE"] = get_pipe_expr(ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ",")))

    # Annotate with new expression
    return ht.annotate(info=ht.info.annotate(**info_expr))
