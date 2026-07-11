"""RQ/Redis implementation of the durable execution queue boundary."""

from __future__ import annotations

from typing import Any

from redis import Redis
from redis.exceptions import RedisError
from rq import Queue
from rq.exceptions import DuplicateJobError, InvalidJobOperation
from rq.job import JobStatus
from rq.serializers import JSONSerializer

from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.services.run_queue import (
    RunQueue,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
    RunQueueUnavailableError,
)
from encode_pipeline.workers.settings import WorkerSettings, load_worker_settings


RESULT_TTL_SECONDS = 86_400
FAILURE_TTL_SECONDS = 604_800
RQ_JOB_CLEANUP_GRACE_SECONDS = 30
REUSABLE_JOB_STATUSES = frozenset(
    {
        JobStatus.QUEUED,
        JobStatus.STARTED,
        JobStatus.DEFERRED,
        JobStatus.SCHEDULED,
    }
)


def rq_job_timeout_seconds(workflow_timeout_seconds: int) -> int:
    """Leave a fixed outer window for durable failure mapping and cleanup."""
    return workflow_timeout_seconds + RQ_JOB_CLEANUP_GRACE_SECONDS


def create_redis_connection(settings: WorkerSettings) -> Redis:
    """Create a Redis client from validated shared worker settings."""
    if not isinstance(settings, WorkerSettings):
        raise ValueError("settings must be a WorkerSettings instance")
    return Redis.from_url(settings.redis_url)


def create_rq_queue(
    settings: WorkerSettings,
    *,
    connection: Any | None = None,
) -> Queue:
    """Create the named RQ queue using the safe JSON serializer."""
    redis_connection = (
        create_redis_connection(settings) if connection is None else connection
    )
    return Queue(
        name=settings.queue_name,
        connection=redis_connection,
        serializer=JSONSerializer,
    )


class RqRunQueue(RunQueue):
    """Enqueue execution identities without process-local dependency objects."""

    def __init__(
        self,
        settings: WorkerSettings | None = None,
        *,
        connection: Any | None = None,
    ) -> None:
        self._settings = settings if settings is not None else load_worker_settings()
        if not isinstance(self._settings, WorkerSettings):
            raise ValueError("settings must be a WorkerSettings instance or None")
        self._owns_connection = connection is None
        self._queue = create_rq_queue(self._settings, connection=connection)

    @property
    def backend(self) -> str:
        """Return the durable backend identity."""
        return "rq"

    @property
    def queue_name(self) -> str:
        """Return the configured backend queue name."""
        return self._settings.queue_name

    def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
        """Enqueue only ``run_id`` under the canonical durable ``job_id``."""
        if not isinstance(assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if assignment.backend != "rq" or assignment.queue_name != self.queue_name:
            raise RunQueueIdentityError(
                "execution assignment does not match the configured RQ queue"
            )
        run_id = _stable_identifier(assignment.run_id, "run_id")
        job_id = _stable_identifier(assignment.job_id, "job_id")

        # Import lazily so the stored callable remains a module-level worker entry
        # point while importing this queue adapter has no worker-runtime side effects.
        from encode_pipeline.workers.jobs import run_execution_job

        try:
            try:
                job = self._queue.enqueue(
                    run_execution_job,
                    args=(run_id,),
                    kwargs={},
                    job_id=job_id,
                    job_timeout=rq_job_timeout_seconds(
                        self._settings.job_timeout_seconds
                    ),
                    result_ttl=RESULT_TTL_SECONDS,
                    failure_ttl=FAILURE_TTL_SECONDS,
                    unique=True,
                )
            except DuplicateJobError:
                job = self._queue.fetch_job(job_id)
                if (
                    job is None
                    or job.func_name != "encode_pipeline.workers.jobs.run_execution_job"
                    or tuple(job.args) != (run_id,)
                    or job.kwargs != {}
                    or job.origin != self.queue_name
                ):
                    raise RunQueueIdentityError(
                        "existing RQ job does not match durable execution identity"
                    ) from None
                try:
                    existing_status = job.get_status(refresh=True)
                except InvalidJobOperation:
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    ) from None
                if existing_status not in REUSABLE_JOB_STATUSES:
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    ) from None
        except RedisError as exc:
            raise RunQueueUnavailableError(
                "RQ could not confirm durable execution submission"
            ) from exc
        return job.id

    def close(self) -> None:
        """Release Redis connection-pool resources owned by this adapter."""
        if self._owns_connection:
            self._queue.connection.close()


def _stable_identifier(value: str, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")
    return value.strip()
