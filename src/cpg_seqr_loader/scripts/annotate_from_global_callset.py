"""Annotate a densified cohort MT with row-level annotations from a global annotate_cohort.mt.

Used by workflows that skip the standard annotation stack (VEP / VQSR / gnomAD / clinvar joins)
because their cohort includes synthetic samples whose AC/AN/AF would pollute the global stats.
Instead of re-annotating from scratch, we join row-wise against the standard seqr-loader's
annotate_cohort.mt (produced from real samples only) and copy every row annotation across.

Design (after multiple attempts at re-key + union + computed-key lookup all hit
LowerDistributedSort failures):
  - The bulk join uses input_mt.row_key directly, matching global_rows' key structure. Where
    possible Hail streams this as a partition-aligned zip - no shuffle.
  - For the ~124 input rows whose exact keys don't appear in global (sparse_split_multi padding
    on the global side), we collect their global annotations driver-side (during recovery), then
    broadcast a small hl.literal dict of overrides.
  - hl.coalesce(direct_join_result, broadcast_override) picks the direct-join annotation where
    it landed, else the broadcast override for the 124 recovered rows.

Consequence: output MT keeps min-rep keys throughout. Those 124 variants appear in seqr under
min-rep form rather than the global's non-minimal padding - arguably more canonical (matches
ClinVar / gnomAD). No variants are lost.

Invariant: every input row must exist in the global (directly or via trim-recovery). If any row
can't be recovered, the recovery step raises loudly - either create_synthetic_proband_gvcf.py has
produced a variant no real sample carries, or temporal drift has moved the parental gVCFs beyond
what the global captured.
"""

import argparse

import loguru
from cpg_utils import config, hail_batch

import hail as hl


def annotate_from_global_callset(
    input_mt_path: str,
    global_mt_path: str,
    output_mt_path: str,
) -> None:
    """Row-join the input densified MT against a global annotate_cohort.mt and write the result."""
    loguru.logger.info(f'Reading input MT: {input_mt_path}')
    input_mt = hl.read_matrix_table(input_mt_path)

    loguru.logger.info(f'Reading global annotate_cohort MT: {global_mt_path}')
    global_mt = hl.read_matrix_table(global_mt_path)
    global_rows = global_mt.rows()

    # Drop row annotations produced by densify - their AC/AN/AF etc. are computed over the
    # synthetic cohort and would be misleading. Global values replace them.
    input_mt = input_mt.drop('info', 'site_dp', 'ANS')

    # Find rewrites for rows whose keys don't exact-match global. Also captures the full global
    # annotation struct for each recovered variant so we can broadcast it later. Raises loudly if
    # any row is unrecoverable (invented variant or temporal drift beyond the global's build).
    rewrites = _build_rewrite_list(input_mt, global_rows)

    loguru.logger.info('Joining global row annotations onto input MT')
    annotated_mt = _annotate_via_direct_join_with_broadcast(input_mt, global_rows, rewrites)

    global_metadata = global_mt.index_globals()
    annotated_mt = annotated_mt.annotate_globals(
        sourceFilePath=input_mt_path,
        genomeVersion=global_metadata.genomeVersion,
        sampleType=global_metadata.sampleType,
        hail_version=hl.__version__,
    )

    loguru.logger.info(f'Writing annotated MT: {output_mt_path}')
    annotated_mt.write(output_mt_path, overwrite=True)


def _trim_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Trim common right bases then common left (shifting position). Mirrors hl.min_rep semantics."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    shift = 0
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        shift += 1
    return pos + shift, ref, alt


def _build_rewrite_list(
    input_mt: hl.MatrixTable,
    global_rows: hl.Table,
    window_bp: int = 50,
) -> list[dict] | None:
    """Find rewrites for input_mt rows that fail exact-match against global.

    For each input row not in global, search nearby global rows (locus +/- window_bp), trim their
    alleles common-bases-style, and find one whose trimmed form matches the input row's key. We
    also capture the matching global row's FULL annotation struct so it can be broadcast later
    (no need for a second Hail lookup on the join path).

    Returns a list of dicts, one per recovered row:
        {
            'orig_contig': str, 'orig_pos': int, 'orig_ref': str, 'orig_alt': str,
            'global_ann': hl.Struct(...)  # global's row-value payload, without locus / alleles
        }
    Or None if no rows need rewriting.

    Raises ValueError if any row is unrecoverable - the workflow must not silently drop rows the
    user expected to appear in seqr.

    Cost: O(mismatches) locus-window queries on global_rows, each using the row-key partition
    index. No dependency on global size.
    """
    missing_rows = input_mt.rows().anti_join(global_rows).select().collect()
    if not missing_rows:
        loguru.logger.info('All input MT rows exact-match global - no key recovery needed')
        return None

    loguru.logger.info(
        f'Attempting key recovery via locus-window search for {len(missing_rows)} unmatched rows',
    )

    rewrites: list[dict] = []
    unrecovered: list[tuple[str, int, str, str]] = []
    for row in missing_rows:
        contig = row.locus.contig
        pos = row.locus.position
        ref, alt = row.alleles[0], row.alleles[1]
        window = hl.locus_interval(
            contig,
            max(1, pos - window_bp),
            pos + window_bp,
            reference_genome='GRCh38',
        )
        # Collect full rows (with annotations) - the matching one gets broadcast back to the MT.
        candidates = global_rows.filter(window.contains(global_rows.locus)).collect()
        matched = False
        for cand in candidates:
            if cand.locus.contig != contig:
                continue
            cand_pos, cand_ref, cand_alt = _trim_alleles(
                cand.locus.position,
                cand.alleles[0],
                cand.alleles[1],
            )
            if cand_pos == pos and cand_ref == ref and cand_alt == alt:
                ann_fields = {k: v for k, v in cand.items() if k not in ('locus', 'alleles')}
                rewrites.append(
                    {
                        'orig_contig': contig,
                        'orig_pos': pos,
                        'orig_ref': ref,
                        'orig_alt': alt,
                        'global_ann': hl.Struct(**ann_fields),
                    },
                )
                matched = True
                break
        if not matched:
            unrecovered.append((contig, pos, ref, alt))

    loguru.logger.info(
        f'Recovered keys for {len(rewrites)}/{len(missing_rows)} unmatched rows via trim-normalisation',
    )

    if unrecovered:
        example_strs = [f'{c}:{p} {r}>{a}' for c, p, r, a in unrecovered[:10]]
        raise ValueError(
            f'AnnotateFromGlobalCallset invariant violation: {len(unrecovered)} variant rows in the '
            f'input MT do not exist in the global annotate_cohort.mt, even after key-recovery via '
            f'trim-normalisation. Either create_synthetic_proband_gvcf.py has produced variants no '
            f'real sample carries, or temporal drift has moved the parental gVCFs beyond what the '
            f'global captured. Investigate one of the unrecovered variants in the parental gVCFs.\n'
            f'First {len(example_strs)} unrecovered variants:\n  ' + '\n  '.join(example_strs),
        )

    return rewrites if rewrites else None


def _annotate_via_direct_join_with_broadcast(
    input_mt: hl.MatrixTable,
    global_rows: hl.Table,
    rewrites: list[dict] | None,
) -> hl.MatrixTable:
    """Direct-key join for 99.999% of rows; broadcast dict fallback for the ~124 recovered ones.

    - `global_rows[input_mt.row_key]` uses the input MT's actual row key. When Hail can prove
      the two tables' keys align, it streams this without shuffling. If it does need to
      co-partition, it's still the best-optimised join case.
    - For the 124 rows whose exact keys don't appear in global, hl.coalesce falls through to a
      broadcast dict of overrides. The dict is a few MB - broadcast to each executor once, then
      per-row lookups are local.
    """
    # Direct join. Most rows get a real annotation struct here; the 124 mismatched rows get null.
    annotated_mt = input_mt.annotate_rows(direct_ann_tmp=global_rows[input_mt.row_key])

    if not rewrites:
        # No recovered rows - just unpack the direct-join struct into the row schema.
        annotated_mt = annotated_mt.annotate_rows(**annotated_mt.direct_ann_tmp)
        return annotated_mt.drop('direct_ann_tmp')

    # Broadcast dict keyed by (contig, pos, ref, alt) - a struct of primitives, safe as a
    # Hail dict key. Value is the full global annotation struct for that recovered variant.
    key_dtype = hl.tstruct(
        contig=hl.tstr,
        pos=hl.tint32,
        ref=hl.tstr,
        alt=hl.tstr,
    )
    override_dict = {
        hl.Struct(
            contig=r['orig_contig'],
            pos=r['orig_pos'],
            ref=r['orig_ref'],
            alt=r['orig_alt'],
        ): r['global_ann']
        for r in rewrites
    }
    override_lit = hl.literal(
        override_dict,
        dtype=hl.tdict(key_dtype, global_rows.row_value.dtype),
    )

    lookup_key = hl.struct(
        contig=annotated_mt.locus.contig,
        pos=annotated_mt.locus.position,
        ref=annotated_mt.alleles[0],
        alt=annotated_mt.alleles[1],
    )

    # Prefer the direct-join result; fall back to the broadcast override for the 124 recovered rows.
    annotated_mt = annotated_mt.annotate_rows(
        combined_ann_tmp=hl.coalesce(
            annotated_mt.direct_ann_tmp,
            override_lit.get(lookup_key),
        ),
    )
    annotated_mt = annotated_mt.annotate_rows(**annotated_mt.combined_ann_tmp)
    return annotated_mt.drop('direct_ann_tmp', 'combined_ann_tmp')


def cli_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to the densified input MT')
    parser.add_argument('--global_mt', required=True, help='Path to the global annotate_cohort.mt to join against')
    parser.add_argument('--output', required=True, help='Path to write the annotated MT')
    args = parser.parse_args()

    # driver_cores + worker_cores tune QoB resourcing. Defaults per Ed's guidance for large-MT
    # joins; override via config keys for retry ladder (2/1 -> 2/2 -> 4/2).
    hail_batch.init_batch(
        driver_cores=config.config_retrieve(['annotate_from_global_callset', 'driver_cores'], 2),
        worker_cores=config.config_retrieve(['annotate_from_global_callset', 'worker_cores'], 1),
    )
    annotate_from_global_callset(
        input_mt_path=args.input,
        global_mt_path=args.global_mt,
        output_mt_path=args.output,
    )


if __name__ == '__main__':
    cli_main()
