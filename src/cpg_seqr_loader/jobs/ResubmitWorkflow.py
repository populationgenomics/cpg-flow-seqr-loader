from typing import TYPE_CHECKING

import toml

from cpg_utils import config, hail_batch, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob


def resubmit_full_workflow(output_path: str) -> 'PythonJob':
    """Chain the second phase of the workflow using analysis-runner."""
    batch_instance = hail_batch.get_batch()
    j = batch_instance.new_bash_job('Submit second Seqr-Loader workflow phase.')
    hail_batch.authenticate_cloud_credentials_in_job(j)
    j.image(config.config_retrieve(['workflow', 'driver_image']))
    config_dict = dict(config._config)  # noqa: SLF001
    _ar = config_dict['workflow'].pop('ar-guid')
    with to_path(output_path).open('w') as output_file:
        output_file.write(toml.dumps(config_dict))

    j.command(f"""
        python -m cpg_seqr_loader.scripts.resubmit_workflow \\
            --dataset {config_dict['workflow']['dataset']} \\
            --image {config_dict['workflow']['driver_image']} \\
            --output {output_path} --script full_workflow \\
            --description '' \\
            --access-level {config_dict['workflow']['access_level']} \\
            --config {output_path} \\
                --skip-repo-checkout
    """)

    return j
