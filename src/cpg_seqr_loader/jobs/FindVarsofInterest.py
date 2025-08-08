from cpg_utils import Path, config, hail_batch


def find_vars_of_interest_job(
    input_mt: str,
    output_tsv: Path,
    job_attrs: dict,
):
    """Take the single-dataset MT, and subset by AF<0.01 to find rare 3' and 5' UTRS OR splice_region_variants AND RETURN A tsv of vars of interest"""

    subset_j = hail_batch.get_batch().new_bash_job('tsv of vars of interest from dataset MT', job_attrs | {'tool': 'hail query'})
    subset_j.image(config.config_retrieve(['workflow', 'driver_image']))
    subset_j.command(f'python -m cpg_seqr_loader.scripts.subset_mt_for_vars_of_interest_alphagenome.py -i {input_mt} --out {output_tsv!s} ')
    return subset_j
