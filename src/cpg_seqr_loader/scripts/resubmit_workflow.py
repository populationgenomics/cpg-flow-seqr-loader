import toml
from argparse import ArgumentParser

from analysis_runner.cli_analysisrunner import run_analysis_runner
from cpg_utils import config, to_path


def submission_python_job(output_path: str):
    """"""

    # make sure the config is loaded
    _ar = config.config_retrieve(['workflow', 'ar-guid'])

    # pop off any workflow-unique attributes
    config_dict = dict(config._config)  # noqa: SLF001
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


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('--output', type=str, required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    submission_python_job(args.output)
