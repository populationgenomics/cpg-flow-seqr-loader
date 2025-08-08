from cpg_utils import Path, config, hail_batch


def alphagenome_var_scanner(
    input_tsv: str,
    output_table: Path,
    output_dir: Path,
    output_summary: Path,
    organ_ontologies: list,
    threshold:int,
    min_length:int,
    merge_distance:int,
    window_size:int,
    scan_span:int,
    apikey: str,
    job_attrs: dict,
):
    """Take the vars of interest as a tsv, run alphagenome scanner, and output a tsv with a scan summary, a tsv with an in-depth output on variants that passed, and a directory with plots containing all predictions for vras of interest."""

    alphagenome_scan = hail_batch.get_batch().new_job('Subsetted MT from dataset MT', job_attrs | {'tool': "alphagenome"})
    alphagenome_scan.image(config.config_retrieve(['workflow', 'driver_image']))
    localised_tsv=hail_batch.get_batch().read_input(input_tsv)


    alphagenome_scan.command(f"""
    python -m cpg_seqr_loader.scripts.alphagenome_splice_variant_scanner.py   --variants {localised_tsv}  --organs {" ".join(organ_ontologies) }   --threshold {threshold} --min-length {min_length} --merge-distance {merge_distance} \
    --window-size {window_size}  --scan-span {scan_span}  --output-table-sum{alphagenome_scan.output_summary}  --output-table {alphagenome_scan.output_table} --output-dir {alphagenome_scan.output_dir}   --api-key {apikey}   --plot-non-sig --scan-all-tracks 
    """)
    hail_batch.get_batch().write_output(alphagenome_scan.output_table,output_table)
    hail_batch.get_batch().write_output(alphagenome_scan.output_summary, output_summary)
    hail_batch.get_batch().write_output(alphagenome_scan.output_dir,output_dir)
    return alphagenome_scan
