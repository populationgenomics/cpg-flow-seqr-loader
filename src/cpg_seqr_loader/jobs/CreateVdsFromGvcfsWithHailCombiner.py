from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_combiner_jobs() -> 'BashJob': ...
