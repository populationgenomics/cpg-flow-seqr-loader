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
from metamist.apis import FamilyApi
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
      - external_family_id is the collaborator-facing ID (e.g. "F000045871") - used in file
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

    Families can carry multiple external IDs (one per source registry, e.g. seqr, molpath). We
    return the first value in that dict; for ravenscroft-rpl each family has only one. Families
    with no external IDs at all are omitted from the map (upstream callers should log-and-skip).

    Uses metamist's FamilyApi REST client rather than a bespoke GraphQL query so we don't
    duplicate what the framework already exposes. Access-level suffixing (`-test`) is applied via
    cpg-flow's `get_metamist_proj` so the same call works at test and standard access levels.
    """
    metamist_proj = get_metamist().get_metamist_proj(dataset_name)
    families = FamilyApi().get_families(project=metamist_proj)
    result: dict[str, str] = {}
    for family in families:
        external_ids = getattr(family, 'external_ids', None) or {}
        if not external_ids:
            continue
        result[str(family.id)] = next(iter(external_ids.values()))
    return result


@functools.cache
def get_families_for_synthetic_probands(
    multicohort: targets.MultiCohort,
) -> list[SyntheticProbandFamily]:
    """Enumerate the duo-only families in the multicohort that qualify for synthetic-proband generation.

    A family qualifies iff all of:
      - exactly two sequencing groups share its family_id
      - one has pedigree sex MALE, one has pedigree sex FEMALE
      - neither member's pedigree references a parent (i.e. both members are themselves parents;
        no real proband is present in the family)
      - both members have a gVCF registered in metamist
      - the family has at least one external ID recorded in metamist (needed for the
        collaborator-facing sample name that ends up in gVCF headers, PED rows, and seqr)

    Families that don't match are logged and skipped, not raised, so a mixed multicohort still
    produces synthetic probands for the families that do qualify.
    """
    # Union family-external-id maps across every dataset represented in the multicohort. Usually
    # there is only one dataset (e.g. ravenscroft-rpl) but the code handles a mixed multicohort.
    external_id_by_internal: dict[str, str] = {}
    for dataset in multicohort.get_datasets():
        external_id_by_internal.update(get_family_external_id_map(dataset.name))

    grouped: dict[str, list[targets.SequencingGroup]] = {}
    for sg in multicohort.get_sequencing_groups():
        family_id = sg.pedigree.fam_id
        if not family_id:
            loguru.logger.warning(f'Skipping SG {sg.id}: no family_id set in pedigree')
            continue
        grouped.setdefault(family_id, []).append(sg)

    families: list[SyntheticProbandFamily] = []
    for family_id, members in grouped.items():
        if len(members) != 2:
            loguru.logger.warning(
                f'Skipping family {family_id}: expected 2 members, got {len(members)}',
            )
            continue

        # Any member with dad or mom set is a child within this family (indicates a real proband),
        # so this isn't a parental duo we should synthesize for.
        if any(m.pedigree.dad or m.pedigree.mom for m in members):
            loguru.logger.warning(
                f'Skipping family {family_id}: at least one member has parents set in the pedigree '
                '(suggests a real proband is present, not a parental duo)',
            )
            continue

        # Compare sex by enum name (Sex.__str__ returns .name) to avoid importing the enum, which
        # is defined in cpg_workflows.targets and not necessarily re-exported by cpg_flow.
        sex_to_sg = {m.pedigree.sex.name: m for m in members}
        if set(sex_to_sg) != {'MALE', 'FEMALE'}:
            observed = [m.pedigree.sex.name for m in members]
            loguru.logger.warning(
                f'Skipping family {family_id}: need one male and one female parent, got {observed}',
            )
            continue

        mother_sg = sex_to_sg['FEMALE']
        father_sg = sex_to_sg['MALE']

        if not mother_sg.gvcf or not father_sg.gvcf:
            missing = [m.id for m in (mother_sg, father_sg) if not m.gvcf]
            loguru.logger.warning(
                f'Skipping family {family_id}: no gVCF registered for {missing}',
            )
            continue

        external_family_id = external_id_by_internal.get(family_id)
        if not external_family_id:
            loguru.logger.warning(
                f'Skipping family {family_id}: no external family ID recorded in metamist',
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
        f'Found {len(families)} duo families eligible for synthetic proband generation',
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
