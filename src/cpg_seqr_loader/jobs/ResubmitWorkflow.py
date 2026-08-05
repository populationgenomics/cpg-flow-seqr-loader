from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def resubmit_full_workflow(output_path: str) -> 'BashJob':
    """Chain the second phase of the workflow using analysis-runner."""
    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Submit second Seqr-Loader workflow phase.')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f'yes | python -m cpg_seqr_loader.scripts.resubmit_workflow {output_path}')
    return job
