"""
suggested location for any utility methods or constants used across multiple stages
"""

import datetime
import functools
import hashlib
from dataclasses import dataclass
from typing import TYPE_CHECKING

import loguru
from cpg_flow import targets
from cpg_flow.metamist import get_metamist
from cpg_utils import Path, config, hail_batch, to_path
from metamist.graphql import gql, query

import hail as hl

if TYPE_CHECKING:
    from hailtop.batch.resource import ResourceGroup

DATE_STRING: str = datetime.datetime.now().strftime('%y-%m')  # noqa: DTZ005


TRAINING_PER_JOB: int = config.config_retrieve(['vqsr', 'vqsr_training_fragments_per_job'])
RECALIBRATION_PER_JOB: int = config.config_retrieve(['vqsr', 'vqsr_apply_fragments_per_job'])
INDEL_RECAL_DISC_SIZE: int = config.config_retrieve(['vqsr', 'indel_recal_disc_size'])
SNPS_RECAL_DISC_SIZE: int = config.config_retrieve(['vqsr', 'snps_recal_disc_size'])
SNPS_GATHER_DISC_SIZE: int = config.config_retrieve(['vqsr', 'snps_gather_disc_size'])

# some file extension constants
VCF_BGZ = 'vcf.bgz'
VCF_BGZ_TBI = 'vcf.bgz.tbi'
VCF_GZ = 'vcf.gz'
VCF_GZ_TBI = 'vcf.gz.tbi'

STANDARD_FEATURES = [
    'ReadPosRankSum',
    'MQRankSum',
    'QD',
    'FS',
    'SOR',
]
SNP_STANDARD_FEATURES = [*STANDARD_FEATURES, 'MQ']
INDEL_STANDARD_FEATURES = STANDARD_FEATURES

ALLELE_SPECIFIC_FEATURES = [
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
    'AS_FS',
    'AS_SOR',
    # Not using depth for the following reasons:
    # 1. The Broad pipelines don't use it;
    # 2. -G AS_StandardAnnotation flag to GenotypeGVCFs doesn't include it;
    # 3. For exomes, depth is an irrelevant feature and should be skipped:
    # 'AS_VarDP'
    # Note that for consistency, we also skip it for WGS.
]
SNP_ALLELE_SPECIFIC_FEATURES = [*ALLELE_SPECIFIC_FEATURES, 'AS_MQ']
INDEL_ALLELE_SPECIFIC_FEATURES = ALLELE_SPECIFIC_FEATURES

SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]


LATEST_ANALYSIS_QUERY = gql(
    """
    query LatestAnalysisEntry($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                meta
                output
                sequencingGroups {
                    id
                }
                timestampCompleted
            }
        }
    }
""",
)

SPECIFIC_VDS_QUERY = gql(
    """
    query getVDSByAnalysisId($vds_id: Int!) {
        analyses(id: {eq: $vds_id}) {
            output
            sequencingGroups {
                id
            }
        }
    }
""",
)

FAMILY_EXTERNAL_IDS_QUERY = gql(
    """
    query FamilyExternalIds($project: String!) {
        project(name: $project) {
            families {
                id
                externalId
            }
        }
    }
""",
)


def read_bed_file_as_intervals(bed_path: str) -> list[hl.Interval]:
    """Manually interpret an input BED file as a series of Intervals."""
    # read intervals BED file manually
    intervals: list[hl.Interval] = []
    with to_path(bed_path).open() as bed_handle:
        for line in bed_handle:
            stripped = line.strip()
            if not stripped:
                continue

            chrom, start, end = stripped.split()[:3]

            start_locus = hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38')
            end_locus = hl.Locus(chrom, int(end), reference_genome='GRCh38')
            intervals.append(hl.Interval(start_locus, end_locus, includes_start=True, includes_end=True))
    return intervals


@functools.cache
def run_annotate_dataset(dataset: str) -> bool:
    """Use all 3 config entries and make a single decision on whether to run the annotate_dataset stage."""
    write_vcf_datasets = config.config_retrieve(['workflow', 'write_vcf'])
    write_es_datasets = config.config_retrieve(['workflow', 'create_es_index_for_datasets'])
    write_mt_datasets = config.config_retrieve(['workflow', 'write_mt_for_datasets'])

    all_datasets = set.union(set(write_vcf_datasets), set(write_es_datasets), set(write_mt_datasets))

    return dataset in all_datasets


def query_for_specific_vds(vds_id: int) -> tuple[str, set[str]] | None:
    """
    query for a specific analysis of type entry_type for a dataset
    if found, return the set of SG IDs in the VDS (using the metadata)

    - stolen from the cpg_workflows.large_cohort.combiner Stage, but duplicated here so we can split pipelines without
      further code changes

    Args:
        vds_id (int): analysis id to query for

    Returns:
        either None if the analysis wasn't found, or a set of SG IDs in the VDS
    """

    # query for the exact, single analysis entry
    query_results: dict[str, dict] = query(SPECIFIC_VDS_QUERY, variables={'vds_id': vds_id})

    if not query_results['analyses']:
        return None
    vds_path: str = query_results['analyses'][0]['output']
    sg_ids = {sg['id'] for sg in query_results['analyses'][0]['sequencingGroups']}
    return vds_path, sg_ids


def query_for_latest_vds(dataset: str, entry_type: str = 'combiner') -> dict | None:
    """
    query for the latest analysis of type entry_type for a dataset
    Args:
        dataset (str): project to query for
        entry_type (str): type of analysis entry to query for
    Returns:
        str, the path to the latest analysis
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(LATEST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': entry_type})

    analyses_by_date = {}

    for analysis in result['project']['analyses']:
        if analysis['output'] and (
            analysis['meta']['sequencing_type'] == config.config_retrieve(['workflow', 'sequencing_type'])
        ):
            analyses_by_date[analysis['timestampCompleted']] = analysis

    if not analyses_by_date:
        loguru.logger.warning(f'No analysis of type {entry_type} found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort as strings
    return analyses_by_date[sorted(analyses_by_date)[-1]]


ANNOTATE_COHORT_STAGE_NAME = 'AnnotateCohort'


def query_for_latest_annotate_cohort_mt(dataset: str) -> str:
    """Find the most recent global AnnotateCohort matrixtable path in metamist.

    Reuses LATEST_ANALYSIS_QUERY (type=matrixtable) and filters Python-side to entries
    produced by this repo's AnnotateCohort stage (via meta.stage), for the current
    workflow's sequencing_type. The meta.stage filter is what excludes legacy runs from
    production-pipelines' AnnotateCohortSmallVariantsWithHailQuery, whose row schema may
    not match what this repo's downstream stages expect.

    Raises ValueError if no matching analysis exists — we don't want the downstream join
    to run against a stale or wrong callset, so this must fail loud at DAG-planning time.
    """
    query_dataset = dataset
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(LATEST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': 'matrixtable'})
    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])

    candidates = {
        analysis['timestampCompleted']: analysis
        for analysis in result['project']['analyses']
        if analysis['output']
        and analysis['meta'].get('stage') == ANNOTATE_COHORT_STAGE_NAME
        and analysis['meta'].get('sequencing_type') == sequencing_type
    }

    if not candidates:
        raise ValueError(
            f'No completed {ANNOTATE_COHORT_STAGE_NAME} matrixtable analysis found in metamist project '
            f'{query_dataset!r} for sequencing_type={sequencing_type!r}. '
            f'The AnnotateFromGlobalCallset stage cannot proceed without a source annotate_cohort.mt.',
        )

    latest = candidates[sorted(candidates)[-1]]
    loguru.logger.info(
        f'Latest global annotate_cohort.mt: {latest["output"]} (completed {latest["timestampCompleted"]})',
    )
    return latest['output']


@functools.lru_cache(1)
def get_localised_resources_for_vqsr() -> dict[str, 'ResourceGroup']:
    """Get the resources required for VQSR, once per run."""

    return {
        key: hail_batch.get_batch().read_input_group(
            base=config.config_retrieve(['references', f'{key}_vcf']),
            index=config.config_retrieve(['references', f'{key}_vcf_index']),
        )
        for key in [
            'axiom_poly',
            'dbsnp',
            'hapmap',
            'mills',
            'omni',
            'one_thousand_genomes',
        ]
    }


@functools.lru_cache(2)
def get_all_fragments_from_manifest(manifest_file: Path) -> list['ResourceGroup']:
    """
    read the manifest file, and return all the fragment resources as an ordered list
    this is a cached method as we don't want to localise every fragment once per task

    Args:
        manifest_file ():

    Returns:
        an ordered list of all the fragment VCFs and corresponding indices
    """

    resource_objects = []
    manifest_folder: Path = manifest_file.parent
    with manifest_file.open() as f:
        for line in f:
            vcf_path = manifest_folder / line.strip()
            resource_objects.append(
                hail_batch.get_batch().read_input_group(
                    **{
                        VCF_GZ: vcf_path,
                        VCF_GZ_TBI: f'{vcf_path}.tbi',
                    }
                ),
            )
    return resource_objects


@functools.cache
def get_family_sequencing_groups(dataset: targets.Dataset) -> dict | None:
    """
    Get the subset of sequencing groups that are in the specified families for a dataset
    Returns a dict containing the sequencing groups and a name suffix for the outputs
    """
    if not config.config_retrieve(['workflow', dataset.name, 'only_families'], []):
        return None
    only_family_ids = set(config.config_retrieve(['workflow', dataset.name, 'only_families'], []))
    # keep only the SG IDs for the families in the only_families list
    loguru.logger.info(f'Finding sequencing groups for families {only_family_ids} in dataset {dataset.name}')
    family_sg_ids = [sg.id for sg in dataset.get_sequencing_groups() if sg.pedigree.fam_id in only_family_ids]
    if not family_sg_ids:
        raise ValueError(f'No sequencing groups found for families {only_family_ids} in dataset {dataset.name}.')
    loguru.logger.info(f'Keeping only {len(family_sg_ids)} SGs from families {len(only_family_ids)} in {dataset}:')
    loguru.logger.info(only_family_ids)
    loguru.logger.info(family_sg_ids)

    h = hashlib.sha256(''.join(sorted(family_sg_ids)).encode()).hexdigest()[:4]
    name_suffix = f'{len(family_sg_ids)}_sgs-{len(only_family_ids)}_families-{h}'

    return {'family_sg_ids': family_sg_ids, 'name_suffix': name_suffix}


# ---------------------------------------------------------------------------
# Synthetic proband workflow helpers (see synthetic_proband_gvcf_workflow.py)
# ---------------------------------------------------------------------------

SYNTHETIC_GVCF_ANALYSIS_TYPE = 'synthetic_gvcf'


@dataclass(frozen=True)
class SyntheticProbandFamily:
    """One duo family in the multicohort, plus the invented name for its synthetic proband.

    Both metamist family IDs are kept:
      - family_id is the internal metamist ID (e.g. "18958") - stable, used for dict keys.
      - external_family_id is the collaborator-facing ID (e.g. "F000012345") - used in file
        names, sample names, PED rows, and Batch job labels so the pre-workflow's artefacts read
        as the identifiers our collaborators already know.

    The synthetic_sample_name string is the single source of truth for the synthetic proband's
    identity across the pre-workflow: it is embedded in the gVCF header (as the sample column
    name), written into the PED file as the proband row's individual ID, and later appears in the
    seqr MatrixTable as a sample ID. The three MUST match or seqr's trio inheritance filters will
    silently fail (pedigree references a proband that isn't present in the callset).
    """

    family_id: str
    external_family_id: str
    mother_sg: targets.SequencingGroup
    father_sg: targets.SequencingGroup
    synthetic_sample_name: str


@functools.cache
def get_family_external_id_map(dataset_name: str) -> dict[str, str]:
    """Return `{internal_family_id_str: external_family_id}` for every family in `dataset_name`.

    cpg-flow's pedigree query passes `replaceWithFamilyExternalIds: false`, so anything the
    framework surfaces (e.g. sg.pedigree.fam_id) is the internal metamist ID. This helper does
    the extra lookup so callers can convert to the external ID collaborators know.

    Uses metamist's GraphQL API rather than the auto-generated FamilyApi REST client. The REST
    client silently drops external IDs recorded under an empty-string source key (which is what
    create_test_subset.py produces), while the GraphQL `externalId` scalar returns them cleanly.

    Families with no external ID at all are omitted from the map (upstream callers should
    log-and-skip). Access-level suffixing (`-test`) is applied via cpg-flow's
    `get_metamist_proj` so the same call works at test and standard access levels.
    """
    metamist_proj = get_metamist().get_metamist_proj(dataset_name)
    result = query(FAMILY_EXTERNAL_IDS_QUERY, variables={'project': metamist_proj})
    return {
        str(family['id']): family['externalId'] for family in result['project']['families'] if family.get('externalId')
    }


def _get_synthetic_proband_couples_config(
    multicohort: targets.MultiCohort,
) -> dict[str, dict]:
    """Read explicit synthetic-proband couples from config, keyed by external family ID.

    Config shape (per dataset):
        [workflow.<dataset>]
        synthetic_proband_couples = [
          { family_id = "F000012345", mother_sg = "CPGxxx", father_sg = "CPGyyy" },
          { family_id = "F000067890", skip = true },
        ]

    Two entry shapes:
      - couple:  {family_id, mother_sg, father_sg} - forces synthesis using those SGs (used for
                 families whose composition isn't a clean 1M+1F duo, e.g. 4-person families or
                 same-sex pairs)
      - skip:    {family_id, skip = true}          - opts a family out of synthesis even when it
                 would otherwise auto-qualify (e.g. a 1M+1F duo we don't want a proband for)
    """
    entries_by_external_fid: dict[str, dict] = {}
    for dataset in multicohort.get_datasets():
        entries = config.config_retrieve(
            ['workflow', dataset.name, 'synthetic_proband_couples'],
            [],
        )
        for entry in entries:
            entries_by_external_fid[entry['family_id']] = entry
    return entries_by_external_fid


def _couple_from_auto_duo(
    members: list[targets.SequencingGroup],
) -> tuple[targets.SequencingGroup, targets.SequencingGroup] | None:
    """Return the (mother, father) pair if `members` is a clean 1M+1F duo, else None."""
    if len(members) != 2:
        return None
    sex_to_sg = {m.pedigree.sex.name: m for m in members}
    if set(sex_to_sg) != {'MALE', 'FEMALE'}:
        return None
    return sex_to_sg['FEMALE'], sex_to_sg['MALE']


def _couple_from_config(
    family_id: str,
    external_family_id: str,
    members: list[targets.SequencingGroup],
    config_entry: dict,
) -> tuple[targets.SequencingGroup, targets.SequencingGroup] | None:
    """Resolve the (mother, father) pair from an explicit config entry.

    Returns None (with a WARNING logged) if the entry is malformed, names SGs that aren't in the
    family, or names SGs with the wrong pedigree sex.
    """
    mother_id = config_entry.get('mother_sg')
    father_id = config_entry.get('father_sg')
    if not mother_id or not father_id:
        loguru.logger.warning(
            f'Skipping family {family_id} ({external_family_id}): '
            'config entry has no mother_sg/father_sg and skip is not set',
        )
        return None

    by_id = {m.id: m for m in members}
    mother_sg = by_id.get(mother_id)
    father_sg = by_id.get(father_id)
    if mother_sg is None or father_sg is None:
        loguru.logger.warning(
            f'Skipping family {family_id} ({external_family_id}): '
            f'config names mother_sg={mother_id!r} father_sg={father_id!r} '
            f'but the family has SGs {sorted(by_id)}',
        )
        return None
    if mother_sg.pedigree.sex.name != 'FEMALE' or father_sg.pedigree.sex.name != 'MALE':
        loguru.logger.warning(
            f'Skipping family {family_id} ({external_family_id}): '
            f'config-selected mother_sg has sex {mother_sg.pedigree.sex.name}, '
            f'father_sg has sex {father_sg.pedigree.sex.name} (need FEMALE / MALE)',
        )
        return None

    return mother_sg, father_sg


def _select_couple_for_family(
    family_id: str,
    external_family_id: str,
    members: list[targets.SequencingGroup],
    config_entry: dict | None,
) -> tuple[targets.SequencingGroup, targets.SequencingGroup] | None:
    """Pick the (mother_sg, father_sg) pair for one family, or None if it should be skipped.

    Skip decisions are logged inside this function so the caller can just filter Nones. The
    caller separately tracks the "multi-member family with no config entry" case for its
    end-of-run summary warning.
    """
    if config_entry and config_entry.get('skip'):
        loguru.logger.info(
            f'Skipping family {family_id} ({external_family_id}): '
            'opted out via synthetic_proband_couples config (skip = true)',
        )
        return None

    if len(members) == 1:
        loguru.logger.warning(
            f'Skipping family {family_id} ({external_family_id}): only 1 SG in the family, '
            'nothing to build a synthetic proband from',
        )
        return None

    if config_entry is None:
        # No override -> auto-qualify a clean 1M+1F duo, or warn+skip anything else.
        couple = _couple_from_auto_duo(members)
        if couple is not None:
            return couple
        loguru.logger.warning(
            f'Skipping family {family_id} ({external_family_id}): '
            f'{len(members)} SGs and no synthetic_proband_couples entry. '
            'Add one under workflow.<dataset>.synthetic_proband_couples to '
            'either include this family or explicitly skip it.',
        )
        return None

    return _couple_from_config(family_id, external_family_id, members, config_entry)


@functools.cache
def get_families_for_synthetic_probands(
    multicohort: targets.MultiCohort,
) -> list[SyntheticProbandFamily]:
    """Enumerate the families in the multicohort that should have a synthetic proband built.

    Selection is composition-driven - affected/phenotype status is deliberately NOT consulted,
    because collaborators may re-label affected status in metamist after the fact.

    Tiers:
      1. If the family appears in `synthetic_proband_couples` with `skip = true`, skip it.
      2. Auto-qualify a family with exactly two SGs when one is MALE and one is FEMALE.
      3. Skip (and log) a family with a single SG - nothing to build from.
      4. For any other composition (>2 SGs, or 2 SGs that are same-sex), look the family up in
         `workflow.<dataset>.synthetic_proband_couples`. If a matching entry names a valid
         (mother_sg, father_sg) pair from within the family, qualify it. Otherwise log a WARNING
         and skip - the family's real SGs still flow through to the combiner and PED.

    In every tier we still require: both chosen parents have a gVCF registered, and the family has
    an external ID in metamist (needed for the collaborator-facing sample name that ends up in
    gVCF headers, PED rows, and seqr).
    """
    # Union family-external-id maps across every dataset represented in the multicohort. Usually
    # there is only one dataset, but the code handles a mixed multicohort.
    external_id_by_internal: dict[str, str] = {}
    for dataset in multicohort.get_datasets():
        external_id_by_internal.update(get_family_external_id_map(dataset.name))

    config_by_external_fid = _get_synthetic_proband_couples_config(multicohort)

    grouped: dict[str, list[targets.SequencingGroup]] = {}
    for sg in multicohort.get_sequencing_groups():
        family_id = sg.pedigree.fam_id
        if not family_id:
            loguru.logger.warning(f'Skipping SG {sg.id}: no family_id set in pedigree')
            continue
        grouped.setdefault(family_id, []).append(sg)

    families: list[SyntheticProbandFamily] = []
    skipped_needs_config: list[str] = []

    for family_id, members in grouped.items():
        external_family_id = external_id_by_internal.get(family_id)
        if not external_family_id:
            loguru.logger.warning(
                f'Skipping family {family_id}: no external family ID recorded in metamist',
            )
            continue

        config_entry = config_by_external_fid.get(external_family_id)
        couple = _select_couple_for_family(family_id, external_family_id, members, config_entry)
        if couple is None:
            if len(members) > 1 and config_entry is None:
                skipped_needs_config.append(external_family_id)
            continue

        mother_sg, father_sg = couple
        if not mother_sg.gvcf or not father_sg.gvcf:
            missing = [m.id for m in (mother_sg, father_sg) if not m.gvcf]
            loguru.logger.warning(
                f'Skipping family {family_id} ({external_family_id}): no gVCF registered for {missing}',
            )
            continue

        families.append(
            SyntheticProbandFamily(
                family_id=family_id,
                external_family_id=external_family_id,
                mother_sg=mother_sg,
                father_sg=father_sg,
                synthetic_sample_name=f'{external_family_id}_synthetic_proband',
            ),
        )

    loguru.logger.info(
        f'Synthetic proband selection: {len(families)} families qualified ({[f.external_family_id for f in families]})',
    )
    if skipped_needs_config:
        loguru.logger.warning(
            f'Synthetic proband selection: {len(skipped_needs_config)} multi-member families '
            f'skipped for lack of a synthetic_proband_couples config entry: '
            f'{skipped_needs_config}',
        )
    return families


def build_gvcf_manifest_content(paths: list[str]) -> str:
    """Newline-separated gVCF paths, in the format Hail's combiner expects via --gvcf_add_file.

    Trailing newline is intentional - the combiner tolerates it and it makes downstream line-count
    / diff tooling less surprising.
    """
    return '\n'.join(paths) + '\n'


_PED_FIELDS = ('Family.ID', 'Individual.ID', 'Father.ID', 'Mother.ID', 'Sex', 'Phenotype')


def build_synthetic_pedigree_content(
    multicohort: targets.MultiCohort,
    families: list[SyntheticProbandFamily],
) -> str:
    """6-column TSV PED content, no header.

    One row per real sequencing group in the multicohort (via cpg-flow's
    sg.pedigree.get_ped_dict, mirroring what multicohort.write_ped_file would produce), plus one
    row per qualifying duo family's synthetic proband (male, affected, referencing its parents).

    Family.ID column is remapped from the internal metamist ID (what cpg-flow surfaces) to the
    external family ID our collaborators know. All rows for the same family MUST share the same
    Family.ID string or seqr won't group them into a trio - the synthetic proband row uses the
    external ID, so the real parent rows must too.

    Row-inclusion is deliberately broader than the qualifying-family set so seqr sees every real
    SG whose gVCF appears in the combiner manifest, whether or not their family qualified for
    synthetic proband synthesis. Non-qualifying real families still get whatever trio linkage
    their metamist pedigree happens to encode.

    The synthetic proband's individual ID (`<external_family_id>_synthetic_proband`) MUST match
    the sample name embedded in the gVCF header and the sample ID that appears in the loaded seqr
    MatrixTable, or seqr will silently disable trio inheritance filtering for the family.
    """
    external_id_by_internal: dict[str, str] = {}
    for dataset in multicohort.get_datasets():
        external_id_by_internal.update(get_family_external_id_map(dataset.name))

    lines: list[str] = []
    for sg in multicohort.get_sequencing_groups():
        ped_dict = sg.pedigree.get_ped_dict()
        # Remap Family.ID from internal → external where we have a mapping; SGs whose family
        # has no external ID fall through with the internal value (they won't be grouped with
        # any synthetic-proband row, which is fine — they're unqualifying families anyway).
        internal_fam = str(ped_dict['Family.ID'])
        family_id_for_ped = external_id_by_internal.get(internal_fam, internal_fam)
        row_values = [family_id_for_ped] + [str(ped_dict[field]) for field in _PED_FIELDS[1:]]
        lines.append('\t'.join(row_values))
    for family in families:
        lines.append(
            f'{family.external_family_id}\t{family.synthetic_sample_name}'
            f'\t{family.father_sg.id}\t{family.mother_sg.id}\t1\t2',
        )
    return '\n'.join(lines) + '\n'


def write_gvcf_manifest(paths: list[str], out_path: str) -> None:
    """Driver-side writer: dump one gVCF path per line to `out_path` (accepts local or gs:// paths)."""
    with to_path(out_path).open('w') as write_handle:
        write_handle.write(build_gvcf_manifest_content(paths))


def manually_find_ids_from_vds(vds_path: str) -> set[str]:
    """
    during development and the transition to input_cohorts over input_datasets, there are some instances
    where we have VDS entries in Metamist, but the analysis entry contains SG IDs which weren't combined into the VDS

    This check bypasses the quick "get all SG IDs in the VDS analysis entry" check,
    and instead checks the exact contents of the VDS

    Args:
        vds_path (str): path to the VDS. Assuming it exists, this will be checked before calling this method

    Returns:
        set[str]: the set of sample IDs in the VDS
    """
    hail_batch.init_batch()
    vds = hl.vds.read_vds(vds_path)

    # find the samples in the Variant Data MT
    return set(vds.variant_data.s.collect())
