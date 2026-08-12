from argparse import ArgumentParser

from analysis_runner.cli_analysisrunner import run_analysis_runner
from loguru import logger


def submit_jobs(
    dataset: str,
    image: str,
    output_dir: str,
    script: str,
    description: str,
    access_level: str,
    config: str,
    skip_repo_checkout: bool = True,
):
    """
    Wrapper script for submitting to analysis-runner
    """

    logger.info(f"Submitting job for dataset {dataset} with description '{description}'")

    run_analysis_runner(
        dataset=dataset,
        image=image,
        output_dir=output_dir,
        script=[script],
        description=description,
        access_level=access_level,
        config=[config],
        skip_repo_checkout=skip_repo_checkout,
    )


if __name__ == '__main__':
    parser = ArgumentParser(description='Submit the second phase of the workflow to analysis-runner')
    parser.add_argument('--dataset', required=True, help='Dataset to submit to analysis-runner')
    parser.add_argument('--image', required=True, help='Docker image to use for the job')
    parser.add_argument('--output', required=True, help='Output directory for the job')
    parser.add_argument('--script', required=True, help='Script to run in the job')
    parser.add_argument('--description', required=True, help='Description of the job')
    parser.add_argument('--access-level', required=True, help='Access level for the job')
    parser.add_argument('--config', required=True, help='Configuration file for the job')
    parser.add_argument('--skip-repo-checkout', action='store_true', help='Skip checking out the repository in the job')
    args = parser.parse_args()

    submit_jobs(
        dataset=args.dataset,
        image=args.image,
        output_dir=args.output,
        script=args.script,
        description=args.description,
        access_level=args.access_level,
        config=args.config,
        skip_repo_checkout=args.skip_repo_checkout,
    )
