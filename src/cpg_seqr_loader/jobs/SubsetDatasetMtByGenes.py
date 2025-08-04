from cpg_utils import Path, config, hail_batch


def subset_to_genes_job(
    input_mt: str,
    output_mt: Path,
    gene_list: list[str],
    job_attrs: dict,
):
    """Take the single-dataset MT, and subset by genes, return subsetted MT."""

    subset_j = hail_batch.get_batch().new_bash_job('Subsetted MT from dataset MT', job_attrs | {'tool': 'hail query'})
    subset_j.image(config.config_retrieve(['workflow', 'driver_image']))
    gene_list_string = ' '.join(gene_list)
    subset_j.command(f'python -m cpg_seqr_loader.scripts.subset_mt_to_genes.py -i {input_mt} --out {output_mt!s} --genes {gene_list_string}')
    return subset_j
