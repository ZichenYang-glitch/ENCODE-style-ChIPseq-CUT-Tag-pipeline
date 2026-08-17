"""RQ/Redis implementation of the durable execution queue boundary."""

from __future__ import annotations

import os
import socket
from datetime import datetime, timezone
from typing import Any

from redis import ConnectionPool, Redis
from redis.connection import parse_url
from redis.exceptions import RedisError
from rq import Queue, Worker
from rq.command import send_stop_job_command
from rq.exceptions import (
    DeserializationError,
    DuplicateJobError,
    InvalidJobOperation,
    NoSuchJobError,
)
from rq.job import Callback, Job, JobStatus
from rq.registry import DeferredJobRegistry, ScheduledJobRegistry
from rq.serializers import JSONSerializer

from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.run_recovery import (
    ExecutionQueueEvidenceState,
    RunExecutionQueueEvidence,
)
from encode_pipeline.services.run_queue import (
    RunQueue,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
    RunRecoveryQueue,
    RunQueueUnavailableError,
    RunQueueStopUnavailableError,
    RunStopQueue,
)
from encode_pipeline.workers.settings import WorkerSettings, load_worker_settings


RESULT_TTL_SECONDS = 86_400
FAILURE_TTL_SECONDS = 604_800
RQ_JOB_STARTUP_ALLOWANCE_SECONDS = 300
RQ_JOB_CLEANUP_GRACE_SECONDS = 30
STOPPED_CALLBACK_PATH = "encode_pipeline.workers.jobs.handle_execution_stopped"
STOPPED_CALLBACK_TIMEOUT_SECONDS = 30
REQUEUE_REQUEST_META_KEY = "helixweave_requeue_requested_at"
REUSABLE_JOB_STATUSES = frozenset(
    {
        JobStatus.QUEUED,
        JobStatus.STARTED,
        JobStatus.DEFERRED,
        JobStatus.SCHEDULED,
    }
)
TERMINAL_JOB_STATUSES = frozenset(
    {
        JobStatus.FINISHED,
        JobStatus.FAILED,
        JobStatus.STOPPED,
        JobStatus.CANCELED,
    }
)
EVIDENCE_STATE_BY_JOB_STATUS = {
    JobStatus.CREATED: ExecutionQueueEvidenceState.UNKNOWN,
    JobStatus.QUEUED: ExecutionQueueEvidenceState.QUEUED,
    JobStatus.DEFERRED: ExecutionQueueEvidenceState.DEFERRED,
    JobStatus.SCHEDULED: ExecutionQueueEvidenceState.SCHEDULED,
    JobStatus.FINISHED: ExecutionQueueEvidenceState.FINISHED,
    JobStatus.FAILED: ExecutionQueueEvidenceState.FAILED,
    JobStatus.STOPPED: ExecutionQueueEvidenceState.STOPPED,
    JobStatus.CANCELED: ExecutionQueueEvidenceState.CANCELED,
}


def rq_job_timeout_seconds(workflow_timeout_seconds: int) -> int:
    """Bound startup, workflow execution, and durable cleanup independently.

    RQ starts this outer deadline when it enters the job's execution body,
    before the platform has resolved the runtime and spawned ``ProcessRunner``.
    The fixed startup allowance prevents normal bounded pre-spawn work from
    consuming the configured workflow deadline or its cleanup grace while the
    complete sum remains a finite worker backstop.
    """
    return (
        RQ_JOB_STARTUP_ALLOWANCE_SECONDS
        + workflow_timeout_seconds
        + RQ_JOB_CLEANUP_GRACE_SECONDS
    )


def create_api_redis_connection(settings: WorkerSettings) -> Redis:
    """Create a bounded-latency Redis client for synchronous API commands."""
    if not isinstance(settings, WorkerSettings):
        raise ValueError("settings must be a WorkerSettings instance")
    return _create_redis_connection(
        settings,
        socket_timeout=settings.redis_api_read_timeout_seconds,
    )


def create_worker_redis_connection(settings: WorkerSettings) -> Redis:
    """Create a worker client without imposing the API command read timeout.

    RQ configures its own socket read timeout from the blocking dequeue interval
    when the worker is constructed. The finite connection timeout remains shared
    so an unavailable Redis endpoint cannot stall worker startup indefinitely.
    """
    if not isinstance(settings, WorkerSettings):
        raise ValueError("settings must be a WorkerSettings instance")
    return _create_redis_connection(settings, socket_timeout=None)


def _create_redis_connection(
    settings: WorkerSettings,
    *,
    socket_timeout: float | None,
) -> Redis:
    options = parse_url(settings.redis_url)
    # Explicit application settings win over optional URL query parameters.
    options.update(
        socket_connect_timeout=settings.redis_connect_timeout_seconds,
        socket_timeout=socket_timeout,
        retry_on_timeout=False,
    )
    return Redis.from_pool(ConnectionPool(**options))


def create_rq_queue(
    settings: WorkerSettings,
    *,
    connection: Any | None = None,
) -> Queue:
    """Create the named RQ queue using the safe JSON serializer."""
    redis_connection = (
        create_api_redis_connection(settings) if connection is None else connection
    )
    return Queue(
        name=settings.queue_name,
        connection=redis_connection,
        serializer=JSONSerializer,
    )


def _fetch_job_read_only(queue: Queue, job_id: str) -> Job | None:
    """Fetch one job hash without allowing RQ to prune a stale queue entry."""
    try:
        return Job.fetch(
            job_id,
            connection=queue.connection,
            serializer=JSONSerializer,
        )
    except NoSuchJobError:
        return None


class RqRunQueue(RunQueue, RunStopQueue, RunRecoveryQueue):
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
        job_meta = _requeue_delivery_meta(assignment)

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
                    on_stopped=Callback(
                        STOPPED_CALLBACK_PATH,
                        timeout=STOPPED_CALLBACK_TIMEOUT_SECONDS,
                    ),
                    meta=job_meta,
                    unique=True,
                )
            except DuplicateJobError:
                job = self._queue.fetch_job(job_id)
                if job is None or not _job_matches_execution_identity(
                    job,
                    run_id=run_id,
                    job_id=job_id,
                    queue_name=self.queue_name,
                ):
                    raise RunQueueIdentityError(
                        "existing RQ job does not match durable execution identity"
                    ) from None
                delivery_match = _requeue_delivery_matches_assignment(
                    job,
                    assignment,
                )
                if delivery_match is None:
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
                if existing_status in {
                    JobStatus.QUEUED,
                    JobStatus.DEFERRED,
                    JobStatus.SCHEDULED,
                } and not _scheduling_status_is_queue_owned(
                    job,
                    status=existing_status,
                    queue=self._queue,
                ):
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    ) from None
                if (
                    existing_status is JobStatus.STARTED
                    and not _has_live_local_started_owner(job)
                ):
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    ) from None
                if assignment.requeue_requested_at is not None:
                    if assignment.requeue_confirmed_at is not None:
                        if not delivery_match:
                            raise RunQueueIdentityError(
                                "existing RQ job does not match durable execution identity"
                            ) from None
                    elif not delivery_match:
                        self._bind_pending_requeue_delivery(assignment)
        except (RunQueueIdentityError, RunQueueJobUnavailableError):
            raise
        except (
            DeserializationError,
            InvalidJobOperation,
            NoSuchJobError,
            OSError,
            RedisError,
            TypeError,
            ValueError,
        ) as exc:
            raise RunQueueUnavailableError(
                "RQ could not confirm durable execution submission"
            ) from exc
        return job.id

    def _bind_pending_requeue_delivery(
        self,
        assignment: RunExecutionAssignment,
    ) -> None:
        """Atomically bind one exact active duplicate to its durable request."""
        job_key = Job.key_for(assignment.job_id)
        with self._queue.connection.pipeline() as pipeline:
            pipeline.watch(job_key)
            job = self._queue.fetch_job(assignment.job_id)
            if job is None or not _job_matches_execution_identity(
                job,
                run_id=assignment.run_id,
                job_id=assignment.job_id,
                queue_name=assignment.queue_name,
            ):
                raise RunQueueIdentityError(
                    "existing RQ job does not match durable execution identity"
                )
            status = job.get_status(refresh=True)
            delivery_match = _requeue_delivery_matches_assignment(job, assignment)
            if delivery_match is None:
                raise RunQueueIdentityError(
                    "existing RQ job does not match durable execution identity"
                )
            if delivery_match:
                return
            if status not in REUSABLE_JOB_STATUSES:
                raise RunQueueJobUnavailableError(
                    "existing RQ job is not in a reusable scheduling state"
                )
            if status is JobStatus.STARTED:
                worker_key = _started_owner_key(job)
                if worker_key is None:
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    )
                pipeline.watch(worker_key)
                if not _has_live_local_started_owner(job):
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    )
            else:
                owner_key = _scheduling_owner_key(
                    job,
                    status=status,
                    queue=self._queue,
                )
                pipeline.watch(owner_key)
                if not _scheduling_status_is_queue_owned(
                    job,
                    status=status,
                    queue=self._queue,
                ):
                    raise RunQueueJobUnavailableError(
                        "existing RQ job is not in a reusable scheduling state"
                    )
            job.meta.update(_requeue_delivery_meta(assignment))
            pipeline.multi()
            job.save(pipeline=pipeline, include_meta=True)
            pipeline.exists(job_key)
            if status is JobStatus.STARTED:
                pipeline.exists(worker_key)
            else:
                pipeline.exists(owner_key)
            pipeline.execute()

    def request_stop(self, assignment: RunExecutionAssignment) -> None:
        """Publish RQ's stop-job command after strict durable identity checks."""
        if not isinstance(assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if (
            assignment.backend != self.backend
            or assignment.queue_name != self.queue_name
        ):
            raise RunQueueIdentityError(
                "execution assignment does not match the configured RQ queue"
            )
        if (
            assignment.claimed_at is None
            or assignment.cancellation_requested_at is None
            or assignment.cancellation_reason is None
        ):
            raise RunQueueStopUnavailableError(
                "RQ could not confirm execution cancellation"
            )

        try:
            job = self._queue.fetch_job(assignment.job_id)
            if job is None or not _job_matches_execution_identity(
                job,
                run_id=assignment.run_id,
                job_id=assignment.job_id,
                queue_name=assignment.queue_name,
            ):
                raise RunQueueStopUnavailableError(
                    "RQ could not confirm execution cancellation"
                )
            status = job.get_status(refresh=True)
            if status is not JobStatus.STARTED:
                raise RunQueueStopUnavailableError(
                    "RQ could not confirm execution cancellation"
                )
            if not isinstance(job.worker_name, str) or not job.worker_name.strip():
                raise RunQueueStopUnavailableError(
                    "RQ could not confirm execution cancellation"
                )
            send_stop_job_command(
                self._queue.connection,
                assignment.job_id,
                serializer=JSONSerializer,
            )
        except RunQueueStopUnavailableError:
            raise
        except (
            DeserializationError,
            InvalidJobOperation,
            NoSuchJobError,
            OSError,
            RedisError,
        ) as exc:
            raise RunQueueStopUnavailableError(
                "RQ could not confirm execution cancellation"
            ) from exc

    def inspect_execution(
        self,
        assignment: RunExecutionAssignment,
    ) -> RunExecutionQueueEvidence:
        """Project RQ state into bounded evidence for one durable identity."""
        if not isinstance(assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if (
            assignment.backend != self.backend
            or assignment.queue_name != self.queue_name
        ):
            return _queue_evidence(ExecutionQueueEvidenceState.IDENTITY_DRIFT)

        try:
            job = _fetch_job_read_only(self._queue, assignment.job_id)
            if job is None:
                return _queue_evidence(ExecutionQueueEvidenceState.MISSING)
            if not _job_matches_execution_identity(
                job,
                run_id=assignment.run_id,
                job_id=assignment.job_id,
                queue_name=assignment.queue_name,
            ):
                return _queue_evidence(ExecutionQueueEvidenceState.IDENTITY_DRIFT)

            requeue_delivery_matches = _requeue_delivery_matches_assignment(
                job,
                assignment,
            )
            if requeue_delivery_matches is None:
                return _queue_evidence(ExecutionQueueEvidenceState.IDENTITY_DRIFT)

            status = job.get_status(refresh=True)
            if status is JobStatus.STARTED:
                state = (
                    ExecutionQueueEvidenceState.STARTED_LIVE
                    if _has_live_local_started_owner(job)
                    else ExecutionQueueEvidenceState.STARTED_UNPROVEN
                )
                return _queue_evidence(
                    state,
                    requeue_delivery_matches_request=requeue_delivery_matches,
                )
            if status in {
                JobStatus.QUEUED,
                JobStatus.DEFERRED,
                JobStatus.SCHEDULED,
            } and not _scheduling_status_is_queue_owned(
                job,
                status=status,
                queue=self._queue,
            ):
                return _queue_evidence(
                    ExecutionQueueEvidenceState.UNKNOWN,
                    requeue_delivery_matches_request=requeue_delivery_matches,
                )
            return _queue_evidence(
                EVIDENCE_STATE_BY_JOB_STATUS.get(
                    status,
                    ExecutionQueueEvidenceState.UNKNOWN,
                ),
                requeue_delivery_matches_request=requeue_delivery_matches,
            )
        except (
            DeserializationError,
            InvalidJobOperation,
            NoSuchJobError,
            OSError,
            RedisError,
            TypeError,
            ValueError,
        ):
            return _queue_evidence(ExecutionQueueEvidenceState.UNAVAILABLE)

    def requeue_execution(self, assignment: RunExecutionAssignment) -> str:
        """Replace a missing or exact terminal job after durable authorization."""
        if not isinstance(assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if (
            assignment.backend != self.backend
            or assignment.queue_name != self.queue_name
        ):
            raise RunQueueIdentityError(
                "execution assignment does not match the configured RQ queue"
            )
        if (
            assignment.dispatched_at is None
            or assignment.claimed_at is not None
            or assignment.cancellation_requested_at is not None
            or assignment.requeue_requested_at is None
            or assignment.requeue_confirmed_at is not None
        ):
            raise RunQueueJobUnavailableError(
                "RQ job is not available for durable execution requeue"
            )

        run_id = _stable_identifier(assignment.run_id, "run_id")
        job_id = _stable_identifier(assignment.job_id, "job_id")
        job_key = Job.key_for(job_id)
        enqueue_after_read = False

        try:
            with self._queue.connection.pipeline() as pipeline:
                pipeline.watch(job_key)
                job = self._queue.fetch_job(job_id)
                if job is None:
                    pipeline.unwatch()
                    enqueue_after_read = True
                else:
                    if not _job_matches_execution_identity(
                        job,
                        run_id=run_id,
                        job_id=job_id,
                        queue_name=self.queue_name,
                    ):
                        raise RunQueueIdentityError(
                            "existing RQ job does not match durable execution identity"
                        )
                    status = job.get_status(refresh=True)
                    delivery_match = _requeue_delivery_matches_assignment(
                        job,
                        assignment,
                    )
                    if delivery_match is None:
                        raise RunQueueIdentityError(
                            "existing RQ job does not match durable execution identity"
                        )
                    if status in REUSABLE_JOB_STATUSES:
                        if status is JobStatus.STARTED:
                            worker_key = _started_owner_key(job)
                            if worker_key is None:
                                raise RunQueueJobUnavailableError(
                                    "RQ job is not available for durable execution requeue"
                                )
                            pipeline.watch(worker_key)
                            if not _has_live_local_started_owner(job):
                                raise RunQueueJobUnavailableError(
                                    "RQ job is not available for durable execution requeue"
                                )
                        else:
                            pipeline.watch(
                                _scheduling_owner_key(
                                    job,
                                    status=status,
                                    queue=self._queue,
                                )
                            )
                        if status is not JobStatus.STARTED and not (
                            _scheduling_status_is_queue_owned(
                                job,
                                status=status,
                                queue=self._queue,
                            )
                        ):
                            raise RunQueueJobUnavailableError(
                                "RQ job is not available for durable execution requeue"
                            )
                        # Execute a watched no-op read to define the point at which
                        # this retry proved that the exact active job still existed.
                        pipeline.multi()
                        if not delivery_match:
                            job.meta.update(_requeue_delivery_meta(assignment))
                            job.save(pipeline=pipeline, include_meta=True)
                        pipeline.exists(job_key)
                        if status is JobStatus.STARTED:
                            pipeline.exists(worker_key)
                        pipeline.execute()
                        return job_id
                    if status not in TERMINAL_JOB_STATUSES:
                        raise RunQueueJobUnavailableError(
                            "RQ job is not available for durable execution requeue"
                        )
                    if delivery_match:
                        pipeline.multi()
                        pipeline.exists(job_key)
                        pipeline.execute()
                        return job_id

                    # The watched job hash covers identity and status changes
                    # between the strict read above and deletion. Registry cleanup
                    # and hash removal commit together, or WatchError preserves it.
                    pipeline.multi()
                    job.delete(pipeline=pipeline)
                    pipeline.execute()
                    enqueue_after_read = True
        except (RunQueueIdentityError, RunQueueJobUnavailableError):
            raise
        except (
            DeserializationError,
            InvalidJobOperation,
            NoSuchJobError,
            OSError,
            RedisError,
            TypeError,
            ValueError,
        ) as exc:
            raise RunQueueUnavailableError(
                "RQ could not confirm durable execution requeue"
            ) from exc

        if not enqueue_after_read:
            raise RunQueueUnavailableError(
                "RQ could not confirm durable execution requeue"
            )
        return self.enqueue_execution(assignment)

    def close(self) -> None:
        """Release Redis connection-pool resources owned by this adapter."""
        if self._owns_connection:
            self._queue.connection.close()


def _stable_identifier(value: str, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")
    return value.strip()


def _queue_evidence(
    state: ExecutionQueueEvidenceState,
    *,
    requeue_delivery_matches_request: bool = False,
) -> RunExecutionQueueEvidence:
    return RunExecutionQueueEvidence(
        state=state,
        requeue_delivery_matches_request=requeue_delivery_matches_request,
    )


def _requeue_delivery_meta(assignment: RunExecutionAssignment) -> dict[str, str]:
    if (
        assignment.requeue_requested_at is None
        or assignment.requeue_confirmed_at is not None
    ):
        return {}
    return {
        REQUEUE_REQUEST_META_KEY: _canonical_requeue_timestamp(
            assignment.requeue_requested_at
        )
    }


def _requeue_delivery_matches_assignment(
    job,
    assignment: RunExecutionAssignment,
) -> bool | None:
    meta = getattr(job, "meta", None)
    if not isinstance(meta, dict):
        return None
    marker = meta.get(REQUEUE_REQUEST_META_KEY)
    if marker is None:
        return False
    if assignment.requeue_requested_at is None or not isinstance(marker, str):
        return None
    try:
        expected = _canonical_requeue_timestamp(assignment.requeue_requested_at)
    except ValueError:
        return None
    return True if marker == expected else None


def _canonical_requeue_timestamp(value: datetime) -> str:
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError("requeue request timestamp must be timezone-aware")
    return (
        value.astimezone(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _has_live_local_started_owner(job) -> bool:
    worker_name = getattr(job, "worker_name", None)
    if not isinstance(worker_name, str) or not worker_name.strip():
        return False
    worker_key = f"{Worker.redis_worker_namespace_prefix}{worker_name}"
    worker_fields = job.connection.hgetall(worker_key)
    if not isinstance(worker_fields, dict) or not worker_fields:
        return False
    current_job = _redis_hash_text(worker_fields, "current_job")
    hostname = _redis_hash_text(worker_fields, "hostname")
    pid_text = _redis_hash_text(worker_fields, "pid")
    if current_job != job.id:
        return False
    if hostname != socket.gethostname() or pid_text is None:
        return False
    try:
        pid = int(pid_text)
    except ValueError:
        return False
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def _started_owner_key(job) -> str | None:
    worker_name = getattr(job, "worker_name", None)
    if not isinstance(worker_name, str) or not worker_name.strip():
        return None
    return f"{Worker.redis_worker_namespace_prefix}{worker_name}"


def _redis_hash_text(values: dict, field_name: str) -> str | None:
    """Read one RQ worker field without invoking registry cleanup helpers."""
    raw = values.get(field_name)
    if raw is None:
        raw = values.get(field_name.encode("utf-8"))
    if isinstance(raw, str):
        return raw
    if isinstance(raw, bytes):
        try:
            return raw.decode("utf-8")
        except UnicodeDecodeError:
            return None
    return None


def _scheduling_status_is_queue_owned(job, *, status: JobStatus, queue: Queue) -> bool:
    if status is JobStatus.QUEUED:
        return job.id in queue.get_job_ids()
    return job.id in _scheduling_registry(job, status=status)


def _scheduling_owner_key(job, *, status: JobStatus, queue: Queue) -> str:
    if status is JobStatus.QUEUED:
        return queue.key
    return _scheduling_registry(job, status=status).key


def _scheduling_registry(job, *, status: JobStatus):
    registry_type = {
        JobStatus.DEFERRED: DeferredJobRegistry,
        JobStatus.SCHEDULED: ScheduledJobRegistry,
    }.get(status)
    if registry_type is None:
        raise ValueError("status does not have an RQ scheduling registry")
    return registry_type(
        name=job.origin,
        connection=job.connection,
        job_class=job.__class__,
        serializer=job.serializer,
    )


def _job_matches_execution_identity(
    job,
    *,
    run_id: str,
    job_id: str,
    queue_name: str,
) -> bool:
    try:
        return (
            job.id == job_id
            and job.func_name == "encode_pipeline.workers.jobs.run_execution_job"
            and tuple(job.args) == (run_id,)
            and job.kwargs == {}
            and job.origin == queue_name
            and getattr(job, "_stopped_callback_name", None) == STOPPED_CALLBACK_PATH
            and job.stopped_callback_timeout == STOPPED_CALLBACK_TIMEOUT_SECONDS
        )
    except (AttributeError, DeserializationError, TypeError, ValueError):
        return False
