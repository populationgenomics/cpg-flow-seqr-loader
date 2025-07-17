import zoneinfo
from datetime import datetime
from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def delete_temp_data_recursive(temp_path: Path, logfile: Path) -> 'BashJob | None':
    """
    takes the path used to write combiner intermediatesm, runs a recursive deletion on that path
    """

    deletion_job = hail_batch.get_batch().new_bash_job('Delete Combiner Intermediates')
    deletion_job.image(config.config_retrieve(['workflow', 'driver_image']))
    deletion_job.command(f'gcloud storage rm -r {temp_path!s}')

    # generate a log file with the deletion timestamp
    # this is generated within the job so that we can ensure the deletion was successful, or no log is written
    timestamp = datetime.now(tz=zoneinfo.ZoneInfo('Australia/Brisbane')).strftime('%Y-%m-%d %H:%M:%S')
    deletion_job.command(f'echo "Deleted {temp_path!s}, date {timestamp}" > {deletion_job.output}')
    hail_batch.get_batch().write_output(deletion_job.output, logfile)

    return deletion_job
