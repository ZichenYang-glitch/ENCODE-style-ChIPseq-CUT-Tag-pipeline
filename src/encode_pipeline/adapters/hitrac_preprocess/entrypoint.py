"""Consume the adapter's identity-bound plan in a dedicated qualification process."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import threading

from encode_pipeline.platform.adapters import WorkflowInputs

from .admission import Sample, _regular
from .calls import strict_json
from .execution import (
    ATTEMPT,
    REQUEST,
    admit_runtime,
    bind_inputs,
    encoded,
    workspace_path,
)
from .qualification import qualify


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    cancelled = threading.Event()
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda *_: cancelled.set())
    result = {"status": "rejected", "reason_code": "execution_binding_invalid"}
    try:
        data = _regular(args.request).read_bytes()
        if hashlib.sha256(data).hexdigest() != args.sha256:
            raise ValueError
        request = strict_json(data.decode())
        if request["schema_version"] != "hitrac-private-plan-v1":
            raise ValueError
        workspace = workspace_path(request["workspace"])
        if Path.cwd() != workspace or args.request != workspace / REQUEST:
            raise ValueError
        runtime = admit_runtime(
            Path(request["runtime_binding"]),
            request["runtime_sha256"],
            timeout=request["timeout"],
        )
        binding = bind_inputs(
            runtime,
            WorkflowInputs(**request["inputs"]),
            {
                "schema_version": "hitrac-reference-profile-v1",
                "binding": request["reference_binding"],
                "sha256": request["reference_sha256"],
            },
        )
        if binding.identity() != request["identity"] or binding.samples_json != encoded(
            request["samples"]
        ):
            raise ValueError
        # Only this fresh owned workspace is restricted; original inputs untouched.
        os.chmod(workspace, 0o700)
        os.chmod(args.request, 0o600)
        inputs = WorkflowInputs(**request["inputs"])
        runtime_checked = runtime.verify()
        expected = {
            "runtime_binding_sha256": runtime_checked.binding_sha256,
            "runtime_lock_sha256": runtime_checked.lock_sha256,
            "reference_binding_sha256": request["reference_sha256"],
            "samples": request["samples"],
        }
        result = qualify(
            runtime_binding=runtime.binding,
            reference_binding=binding.reference,
            reference_sha256=binding.reference_sha256,
            samples=[
                Sample(row["sample_id"], Path(row["fastq_1"]), Path(row["fastq_2"]))
                for row in inputs.samples
            ],
            attempt=workspace / ATTEMPT,
            threads=inputs.options["threads"],
            mapq=inputs.options["mapq"],
            timeout=runtime.timeout,
            cancelled=cancelled.is_set,
            expected_input_identity=expected,
            expected_implementation_sha256=request["identity"]["implementation"],
        )
    except (OSError, ValueError, TypeError, KeyError):
        pass
    # No paths, raw exception text or scientific logs enter platform log streams.
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key in {"status", "reason_code", "sample", "stage", "collection"}
            },
            sort_keys=True,
        )
    )
    return 0 if result["status"] == "complete" else 1
