from typing import TYPE_CHECKING

from cpg_utils import hail_batch, config


if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob


def submission_python_job(output_path: str):
    """"""
    import toml
    from analysis_runner.cli_analysisrunner import run_analysis_runner
    from cpg_utils import config, to_path

    # make sure the config is loaded
    _ar = config.config_retrieve(['workflow', 'ar-guid'])

    # pop off any workflow-unique attributes
    config_dict = dict(config._config)
    _ar = config_dict['workflow'].pop('ar-guid')

    with open('config.toml', 'w') as config_file:
        config_file.write(toml.dumps(config_dict))

    run_analysis_runner(
        dataset=config_dict['workflow']['dataset'],
        image=config_dict['workflow']['driver_image'],
        output_dir='',
        script=['full_workflow'],
        description='',
        access_level=config_dict['workflow']['access_level'],
        config=['config.toml'],
        skip_repo_checkout=True,
    )

    with to_path(output_path).open('w') as output_file:
        output_file.write(toml.dumps(config_dict))


def resubmit_full_workflow(output_path: str) -> 'PythonJob':
    """Chain the second phase of the workflow using analysis-runner."""
    batch_instance = hail_batch.get_batch()
    npj = batch_instance.new_python_job('Submit second Seqr-Loader workflow phase.')
    npj.image(config.config_retrieve(['workflow', 'driver_image']))
    npj.call(submission_python_job, output_path)
    return npj
