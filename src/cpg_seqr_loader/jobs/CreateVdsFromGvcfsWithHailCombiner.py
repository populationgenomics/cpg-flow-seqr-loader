from typing import TYPE_CHECKING

from cpg_utils import hail_batch, config
from cpg_flow import targets

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_combiner_jobs() -> 'BashJob':
    ...

