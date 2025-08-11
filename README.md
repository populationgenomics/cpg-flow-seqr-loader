# cpg-flow-seqr-loader

## Purpose

This repository is a CPG-Flow migration of the [RD_Combiner](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/stages/rd_combiner.py) workflow, which is an evolution of the original [Seqr Loader](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/stages/seqr_loader.py) workflow.

The intention of these workflows is to continue on from the single-sample workflow (alignment, variant calling), and carry out the multi-sample steps of the analysis, including:

- Using the Hail gVCF Combiner to create a VDS representing the whole joint callset.
- Densifying the VDS representation to create a MatrixTable and a VCF
- Joining the fragments of the VCF into a single genome/exome-wide VCF
- Annotating the VCF with gene and transcript consequences (VEP)
- Carrying our Variant quality recalibration (VQSR) on the VCF
- Combining VEP and VQSR annotations into a full MatrixTable
- Preparing that MatrixTable for loading into Seqr (as an ElasticSearch Index)

## Structure

This repository attempts to follow a strict segregation of concerns, with the following structure:

```commandline
src
├── cpg_seqr_loader
│   ├── config_template.toml
│   ├── hail_scripts
│   │   ├── README.md
│   │   ├── annotations.py
│   │   ├── ...
│   ├── jobs
│   │   ├── AnnotateCohort.py
│   │   ├── ...
│   ├── run_workflow.py
│   ├── scripts
│   │   ├── annotate_cohort.py
│   │   ├── ...
│   ├── stages.py
│   └── utils.py
```

- `config_template.toml` is a template for the configuration file, which is used to set up the workflow. It contains all the settings and references required to run the workflow, and can be copied and modified to create a specific configuration for a given run.
- `stages.py` contains the stages of the workflow, which are the individual steps that make up the workflow. Each stage is defined as a class, and the logic for each stage is imported from the `jobs` directory. A Stage defined in this file states the inputs, outputs, and how the Stage is positioned in the overall workflow graph. All other business logic should be imported from a corresponding file in the `jobs` directory.
- `jobs` contains a directory of Python files, each of which contains the logic for a specific stage of the workflow. These are imported into `stages.py` to define the stages of the workflow.
- `scripts` contains the scripts that are run by the workflow. Some Stages run command-line tools (e.g. [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html)), so the Stage's operations are triggered by interacting with a CLI in a specific container. Where the Stage's logic is entirely bespoke code (e.g. a Python script), it is placed in the `scripts` directory, and baked into the container image.
- `hail_scripts` contains the Hail/Gnomad/Seqr-specific scripts used in the workflow. A more detailed README is provided in that directory.
- `run_workflow.py` is the entry point for the workflow. It is responsible for setting up the workflow, and running the stages in the correct order.

## Operation

This is a brief README, and is a stand-in for a more comprehensive SOP.

CPG-Flow workflows are operated entirely by defining input Cohorts (see [here](https://github.com/populationgenomics/team-docs/blob/main/metamist/cohorts.md)). Once a collection of input Cohorts has been defined, the workflow can be run in two steps as follows:

1. Run the first part of the workflow. This combines new gVCFs into the existing VDS, and densifies the resulting data structure into a MatrixTable, which is exported as a collection of sites-only VCFs to be processed an annotated by tools which require a VCF input. This is done by running the `first_workflow` command:

```bash
analysis-runner \
    --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-seqr-loader:0.1.6 \
    --config src/cpg_seqr_loader/config_template.toml \
    --config cohorts.toml \  # containing the inputs_cohorts and sequencing_type
    --dataset seqr \
    --description 'seqr_loader' \
    --access-level full \
    --output-dir seqr_loader \
  first_workflow
```

2. Run the remainder of the workflow, using the same group of Cohorts. This begins from the exported VCF fragments, recombines all single VCFs into a whole-genome/exome VCF, then runs all annotation and variant quality recalibration steps. This is done by running the `full_workflow` command:

```bash
analysis-runner \
    --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-seqr-loader:0.1.6 \
    --config src/cpg_seqr_loader/config_template.toml \
    --config cohorts.toml \  # containing the inputs_cohorts and sequencing_type
    --dataset seqr \
    --description 'seqr_loader' \
    --access-level full \
    --output-dir seqr_loader \
  full_workflow
```

> **NOTE**: if you run the `full_workflow` without first running the `first_workflow`, it will fail. This is because the `full_workflow` command expects the VCF fragments to already exist in the output directory. This is a little frustrating, but Hail Batch's workflow planning has no way of reactively altering the number of jobs in a workflow graph once the workflow has been initiated. The exact number of fragments the VDS will be exported as is not known until the first workflow has been run, so the second workflow cannot be planned until the first has completed. Forcing the data into a specific number of partitions could control this (`config.workflow.densify_partitions`, called during the VDS -> MT densification process), but this isn't exactly respected when VCFs are exported, so it is not a reliable solution.

Once the workflow has run to completion, there will be an AnnotateCohort MatrixTable in the tmp directory. Action of further stages is controlled via configuration options:

* `workflow.write_mt_for_datasets`: the `SubsetMtToDatasetWithHail` and `AnnotateDataset` Stages will only run for Datasets which are specified in this list. If the list is empty, no further MatrixTable will be written.
* `workflow.write_vcf`: the `AnnotatedDatasetMtToVcf` Stage will only run for datasets in this list. This generates a multisample VCF for each required dataset, which can be shared with collaborators or used for further analysis.
* `workflow.create_es_index_for_datasets`: the `ExportMtAsEsIndex` Stage will only run for datasets in this list. This generates an ElasticSearch index for each required dataset, which can be used to load the data into Seqr.

If a dataset appears in any one of these 3 lists, `AnnotateDataset` will run. Control of the final two stages is specific to the controlling config entry for each.

By default, none of these stages will run, which gives us granular control over the data being generated and pushed into long-term storage. Execution of `ExportMtAsEsIndex` is also conditional on the ElasticSearch credentials being successfully accessed as GCP Secrets at runtime, so if the service account running the workflow does not have access to the secrets, or the required configuration entries are absent, this stage will not run.
