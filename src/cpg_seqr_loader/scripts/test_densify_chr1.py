"""
Test the post-densification pipeline on chr1 of a real checkpoint MT.

Reads an existing densification checkpoint (pre-info, pre-split),
filters to chr1, then runs two approaches:

  1. CURRENT: info → adjust_vcf_types → sparse_split_multi (no info adjustment)
  2. PROPOSED: info → sparse_split_multi(filter_changed_loci) → adjust_info → adjust_vcf_types

Reports diagnostics for both: schema, row counts, drop accounting,
and sample values so you can compare without running the full pipeline.

Usage:
    python -m cpg_seqr_loader.scripts.test_densify_chr1 \
        --checkpoint gs://path/to/densified_checkpoint.mt \
        --output gs://path/to/test_output_chr1.mt
"""

import argparse

import loguru
from cpg_utils import config, hail_batch

import hail as hl

from cpg_seqr_loader.hail_scripts.sparse_mt import default_compute_info
from cpg_seqr_loader.hail_scripts.vcf import adjust_vcf_incompatible_types


def _triangle_index(j: int, k: int) -> int:
    if k < j:
        j, k = k, j
    return k * (k + 1) // 2 + j


def _downcode_pl(lpl: list[int], la: list[int], a_index: int) -> list[int] | None:
    """Manually downcode LPL to biallelic PL for split allele a_index.

    For each biallelic genotype g in {0=hom-ref, 1=het, 2=hom-alt},
    PL[g] = min over all local genotypes (j,k) that downcode to g of LPL[tri(j,k)].
    Downcoding maps: the local allele at position la.index(a_index) -> 1 (alt),
    everything else -> 0 (ref).
    """
    if a_index not in la:
        return None
    local_a = la.index(a_index)
    n_local = len(la)
    pl: list[int | None] = [None, None, None]
    for j in range(n_local):
        for k in range(j, n_local):
            idx = _triangle_index(j, k)
            if idx >= len(lpl):
                continue
            dj = 1 if j == local_a else 0
            dk = 1 if k == local_a else 0
            biallelic_g = dj + dk
            if pl[biallelic_g] is None or lpl[idx] < pl[biallelic_g]:
                pl[biallelic_g] = lpl[idx]
    return pl


def verify_lpl_to_pl(
    mt_pre: hl.MatrixTable,
    mt_post: hl.MatrixTable,
    label: str,
    n_samples: int = 5,
):
    """Verify LPL->PL conversion at a multiallelic site by comparing
    pre-split local allele entries against post-split biallelic entries."""

    multi_rows = mt_pre.filter_rows(hl.len(mt_pre.alleles) > 2).rows()
    site = multi_rows.take(1)
    if not site:
        loguru.logger.warning(f'{label}: no multiallelic sites for LPL->PL verification')
        return

    locus = site[0].locus
    alleles = list(site[0].alleles)
    loguru.logger.info(f'\n=== {label}: LPL->PL verification at {locus} {alleles} ===')

    col_key_field = list(mt_pre.col_key.dtype.keys())[0]

    # Pre-split entries at this locus (prefer non-ref samples for interesting data)
    pre_ht = mt_pre.filter_rows(mt_pre.locus == locus).entries()
    pre_nonref = pre_ht.filter(hl.is_defined(pre_ht.LGT) & pre_ht.LGT.is_non_ref())
    pre_data = pre_nonref.take(n_samples)
    if len(pre_data) < 2:
        pre_data = pre_ht.take(n_samples)
    if not pre_data:
        loguru.logger.warning(f'{label}: no entries at {locus}')
        return

    sample_ids = [getattr(e, col_key_field) for e in pre_data]

    loguru.logger.info(f'  Pre-split ({len(sample_ids)} samples):')
    for e in pre_data:
        sid = getattr(e, col_key_field)
        la = list(e.LA) if e.LA is not None else None
        lpl = list(e.LPL) if e.LPL is not None else None
        lgt = e.LGT
        lad = list(e.LAD) if e.LAD is not None else None
        loguru.logger.info(f'    {sid}: LA={la} LGT={lgt} LPL={lpl} LAD={lad} DP={e.DP}')

    # Post-split entries at same locus for same samples
    post_ht = mt_post.filter_rows(mt_post.locus == locus).entries()
    sample_set = hl.literal(set(sample_ids))
    post_ht = post_ht.filter(sample_set.contains(post_ht[col_key_field]))
    post_data = post_ht.collect()

    if not post_data:
        loguru.logger.warning(
            f'{label}: no post-split entries at {locus} (locus may have shifted via min_rep)'
        )
        return

    by_a_index: dict[int, list] = {}
    for e in post_data:
        by_a_index.setdefault(e.a_index, []).append(e)

    n_checks = 0
    n_pass = 0
    n_fail = 0

    for a_idx in sorted(by_a_index.keys()):
        loguru.logger.info(f'\n  Post-split a_index={a_idx} (allele: {alleles[a_idx]}):')
        for e in by_a_index[a_idx]:
            sid = getattr(e, col_key_field)
            pl = list(e.PL) if e.PL is not None else None
            gt = e.GT
            ad = list(e.AD) if e.AD is not None else None
            loguru.logger.info(f'    {sid}: GT={gt} PL={pl} AD={ad}')

            pre_e = next((p for p in pre_data if getattr(p, col_key_field) == sid), None)
            if pre_e is None or pre_e.LPL is None or pre_e.LA is None:
                continue

            la = list(pre_e.LA)
            lpl = list(pre_e.LPL)
            expected_pl = _downcode_pl(lpl, la, a_idx)
            n_checks += 1

            if expected_pl is None:
                if pl is None:
                    loguru.logger.info(f'      [PASS] allele {a_idx} not in LA={la}, PL correctly missing')
                    n_pass += 1
                else:
                    loguru.logger.error(f'      [FAIL] allele {a_idx} not in LA={la}, expected missing PL, got {pl}')
                    n_fail += 1
            else:
                if pl is not None and list(pl) == expected_pl:
                    loguru.logger.info(f'      [PASS] PL={pl} matches manual downcode from LPL={lpl} LA={la}')
                    n_pass += 1
                else:
                    loguru.logger.error(
                        f'      [FAIL] expected PL={expected_pl} from LPL={lpl} LA={la}, got {pl}'
                    )
                    n_fail += 1

            if pl is not None:
                if len(pl) != 3:
                    loguru.logger.error(f'      [FAIL] PL has {len(pl)} elements, expected 3')
                elif min(pl) != 0:
                    loguru.logger.warning(f'      [WARN] min(PL)={min(pl)}, expected 0 (normalization)')

                if gt is not None:
                    gt_idx = gt[0] + gt[1]
                    if pl[gt_idx] != min(pl):
                        loguru.logger.warning(
                            f'      [WARN] GT={gt} -> PL[{gt_idx}]={pl[gt_idx]} '
                            f'but min(PL)={min(pl)} at index {pl.index(min(pl))}'
                        )

            if ad is not None and pre_e.DP is not None:
                ad_sum = sum(ad)
                if ad_sum > pre_e.DP:
                    loguru.logger.warning(f'      [WARN] sum(AD)={ad_sum} > DP={pre_e.DP}')

    loguru.logger.info(f'\n  {label} LPL->PL summary: {n_pass} passed, {n_fail} failed out of {n_checks} checks')


def report_schema(mt: hl.MatrixTable, label: str):
    loguru.logger.info(f'\n=== {label}: info field schema ===')
    for field, dtype in mt.info.dtype.items():
        loguru.logger.info(f'  info.{field}: {dtype}')

    if 'AS_lowqual' in mt.row:
        loguru.logger.info(f'  AS_lowqual: {mt.AS_lowqual.dtype}')
    if 'lowqual' in mt.row:
        loguru.logger.info(f'  lowqual: {mt.lowqual.dtype}')

    loguru.logger.info(f'\n=== {label}: entry field schema ===')
    for field, dtype in mt.entry.dtype.items():
        loguru.logger.info(f'  {field}: {dtype}')


def report_row_counts(mt: hl.MatrixTable, label: str):
    n_rows = mt.count_rows()
    n_multiallelic = mt.aggregate_rows(hl.agg.count_where(mt.was_split))
    n_biallelic = n_rows - n_multiallelic
    loguru.logger.info(
        f'{label}: {n_rows} total rows '
        f'({n_biallelic} originally biallelic, {n_multiallelic} from multiallelic splits)'
    )
    return n_rows


def adjust_info_after_split(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Post-split: convert allele-specific info arrays to scalars using a_index."""

    SCALAR_INFO_FIELDS = {
        'QUALapprox', 'VarDP', 'ReadPosRankSum', 'MQRankSum',
        'MQ', 'QD', 'FS', 'SOR', 'SB', 'DP',
    }

    info_updates = {}
    for field, dtype in mt.info.dtype.items():
        if field in SCALAR_INFO_FIELDS:
            continue
        elif field == 'AS_SB_TABLE':
            info_updates[field] = [mt.info.AS_SB_TABLE[0], mt.info.AS_SB_TABLE[mt.a_index]]
        elif isinstance(dtype, hl.tarray):
            info_updates[field] = mt.info[field][mt.a_index - 1]

    mt = mt.annotate_rows(info=mt.info.annotate(**info_updates))

    if 'AS_lowqual' in mt.row:
        mt = mt.annotate_rows(AS_lowqual=mt.AS_lowqual[mt.a_index - 1])

    return mt


def compute_info(mt: hl.MatrixTable) -> tuple[hl.MatrixTable, hl.Table]:
    """Compute info fields and annotate onto MT. Returns (annotated MT, info_ht)."""
    info_ht = default_compute_info(mt, site_annotations=True, n_partitions=mt.n_partitions())
    info_ht = info_ht.annotate(info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp))
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
    mt = mt.drop('gvcf_info')
    return mt, info_ht


def run_current_approach(mt: hl.MatrixTable):
    """Mirrors what's on main right now: adjust_vcf_types BEFORE split, no info adjustment after."""
    loguru.logger.info('\n====== CURRENT APPROACH ======')

    mt, info_ht = compute_info(mt)

    info_ht = adjust_vcf_incompatible_types(info_ht, pipe_delimited_annotations=[])
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    pre_split_rows = mt.count_rows()
    expected_rows = mt.aggregate_rows(hl.agg.sum(hl.len(mt.alleles) - 1))
    loguru.logger.info(f'Pre-split: {pre_split_rows} rows, expecting {expected_rows} after split')

    mt_pre_split = mt
    mt = hl.experimental.sparse_split_multi(mt)

    verify_lpl_to_pl(mt_pre_split, mt, 'CURRENT')
    report_schema(mt, 'CURRENT (post-split)')
    actual_rows = report_row_counts(mt, 'CURRENT')

    dropped = expected_rows - actual_rows
    loguru.logger.info(f'CURRENT: {dropped} variants dropped by min_rep locus shift')

    return mt


def run_proposed_approach(mt: hl.MatrixTable):
    """Proposed fix: split FIRST, adjust info AFTER, then VCF types last."""
    loguru.logger.info('\n====== PROPOSED APPROACH ======')

    mt, _ = compute_info(mt)

    pre_split_rows = mt.count_rows()
    expected_rows = mt.aggregate_rows(hl.agg.sum(hl.len(mt.alleles) - 1))
    loguru.logger.info(f'Pre-split: {pre_split_rows} rows, expecting {expected_rows} after split')

    mt_pre_split = mt
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    verify_lpl_to_pl(mt_pre_split, mt, 'PROPOSED')

    actual_rows = mt.count_rows()
    dropped = expected_rows - actual_rows
    if dropped > 0:
        loguru.logger.warning(f'{dropped} variants dropped by min_rep locus shift during split')
    else:
        loguru.logger.info('0 variants dropped during min_rep locus shift')

    mt = adjust_info_after_split(mt)

    report_schema(mt, 'PROPOSED (post-split, post-info-adjust)')

    info_ht = mt.rows().select('info')
    info_ht = adjust_vcf_incompatible_types(info_ht, pipe_delimited_annotations=[])
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    report_schema(mt, 'PROPOSED (post-VCF-adjust)')
    report_row_counts(mt, 'PROPOSED')

    return mt


def main(checkpoint_path: str, output_path: str | None = None):
    hail_batch.init_batch(
        worker_memory=config.config_retrieve(['combiner', 'worker_memory']),
        worker_cores=config.config_retrieve(['combiner', 'worker_cores']),
        driver_memory=config.config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config.config_retrieve(['combiner', 'driver_cores']),
    )

    loguru.logger.info(f'Reading checkpoint: {checkpoint_path}')
    mt = hl.read_matrix_table(checkpoint_path)

    n_samples = mt.count_cols()
    n_rows_full = mt.count_rows()
    loguru.logger.info(f'Full checkpoint: {n_rows_full} rows, {n_samples} samples')
    loguru.logger.info(f'Entry schema: {dict(mt.entry.dtype.items())}')

    loguru.logger.info('Filtering to chr1...')
    mt_chr1 = mt.filter_rows(mt.locus.contig == 'chr1')
    n_chr1 = mt_chr1.count_rows()
    loguru.logger.info(f'chr1: {n_chr1} rows ({n_chr1/n_rows_full*100:.1f}% of total)')

    n_multiallelic = mt_chr1.aggregate_rows(hl.agg.count_where(hl.len(mt_chr1.alleles) > 2))
    loguru.logger.info(f'chr1 multiallelic sites: {n_multiallelic}')

    max_alleles = mt_chr1.aggregate_rows(hl.agg.max(hl.len(mt_chr1.alleles)))
    loguru.logger.info(f'chr1 max allele count at a single site: {max_alleles}')

    run_current_approach(mt_chr1)

    mt_proposed = run_proposed_approach(mt_chr1)

    if output_path:
        loguru.logger.info(f'Writing proposed chr1 result to {output_path}')
        mt_proposed.write(output_path, overwrite=True)
        loguru.logger.info('Done.')


def cli_main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--checkpoint', required=True, help='Path to the densified checkpoint MT')
    parser.add_argument('--output', help='Optional: write the proposed-fix chr1 MT here')
    args = parser.parse_args()
    main(checkpoint_path=args.checkpoint, output_path=args.output)


if __name__ == '__main__':
    cli_main()