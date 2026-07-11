"""Command-line entry point for the local RQ worker process."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from rq import Worker
from rq.serializers import JSONSerializer

from encode_pipeline.workers.rq_queue import (
    create_redis_connection,
    create_rq_queue,
)
from encode_pipeline.workers.settings import load_worker_settings


def build_parser() -> argparse.ArgumentParser:
    """Build the stable worker command-line parser."""
    parser = argparse.ArgumentParser(description="Run the ENCODE platform worker.")
    parser.add_argument(
        "--burst",
        action="store_true",
        help="Process currently queued jobs and exit.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Start one RQ worker using the shared environment configuration."""
    args = build_parser().parse_args(argv)
    settings = load_worker_settings()
    connection = create_redis_connection(settings)
    try:
        queue = create_rq_queue(settings, connection=connection)
        worker = Worker(
            [queue],
            connection=connection,
            serializer=JSONSerializer,
        )
        worker.work(burst=args.burst)
    finally:
        connection.close()
    return 0


if __name__ == "__main__":  # pragma: no cover - console-script fallback
    raise SystemExit(main())
