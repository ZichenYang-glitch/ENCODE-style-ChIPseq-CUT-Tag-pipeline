"""Private, offline qualification of fixed tracPre2; no registry or API binding."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import sys
import threading
import traceback
from collections.abc import Callable

from .admission import (
    AdmissionError,
    Sample,
    load_reference_binding,
    load_runtime_binding,
    prepare_inputs,
)
from .calls import (
    CallFailure,
    digest,
    plan_calls,
    prepare_shims,
    strict_json,
    verify_calls,
    write_exclusive,
)
from .outputs import OutputRejected, iter_bedpe, verify_outputs
from .pairs import PairError, verify_pairs
from .process import execute


class QualificationError(Exception):
    def __init__(self, code: str, sample=None, stage=None, collection=None):
        self.code, self.sample, self.stage, self.collection = (
            code,
            sample,
            stage,
            collection,
        )
        super().__init__(code)


def implementation_identity() -> dict:
    package = Path(__file__).parent
    root = package.parents[3]
    files = sorted(package.glob("*.py")) + [
        root / "config/hitrac_preprocess/tools.lock.json",
        root / "config/hitrac_preprocess/conda-explicit.lock",
        root / "scripts/qualify_hitrac_preprocess.py",
        root / "scripts/prepare_hitrac_bindings.py",
        root / "scripts/run_hitrac_preprocess.py",
    ]
    # Bind actual shared execution and H4 publication/download consumers. This
    # closure is rechecked after science and again before result publication.
    shared = (
        "platform/adapters.py",
        "platform/builds.py",
        "platform/registry.py",
        "platform/planning.py",
        "platform/snapshots.py",
        "platform/reference_profiles.py",
        "platform/result_generations.py",
        "platform/artifact_publications.py",
        "platform/results.py",
        "platform/runs.py",
        "platform/execution.py",
        "services/defaults.py",
        "services/validated_inputs.py",
        "services/planning.py",
        "services/materialization.py",
        "services/command_builder.py",
        "services/workflow_builds.py",
        "services/reference_profile_runtime.py",
        "services/private_reference_profiles.py",
        "services/workflow_info.py",
        "services/validation.py",
        "services/local_execution.py",
        "services/local_run_driver.py",
        "services/preflight.py",
        "services/run_submission.py",
        "services/run_cancellation.py",
        "services/runs.py",
        "services/run_repositories.py",
        "services/process_runner.py",
        "services/artifact_extraction.py",
        "services/qc_summary_indexing.py",
        "services/artifact_downloads.py",
        "services/artifact_publications.py",
        "persistence/repositories.py",
        "persistence/models.py",
        "persistence/runtime.py",
        "api/models.py",
        "workers/timeouts.py",
        "workers/jobs.py",
        "workers/runtime.py",
        "workers/settings.py",
        "workers/rq_queue.py",
    )
    files.extend(root / "src/encode_pipeline" / path for path in shared)
    files.extend(
        root / "scripts" / path
        for path in ("checkout_bootstrap.py", "source_provenance.py")
    )
    values = {str(p.relative_to(root)): digest(p) for p in files}
    encoded = json.dumps(values, sort_keys=True, separators=(",", ":")).encode()
    return {"files": values, "sha256": hashlib.sha256(encoded).hexdigest()}


def qualify(
    *,
    runtime_binding: Path,
    reference_binding: Path,
    reference_sha256: str,
    samples: list[Sample],
    attempt: Path,
    threads: int = 2,
    mapq: int = 10,
    timeout: float = 300,
    cancelled: Callable[[], bool] = lambda: False,
    expected_input_identity: dict | None = None,
    expected_implementation_sha256: str | None = None,
) -> dict:
    """Create one private attempt; never reuse, publish partially or delete outputs.

    timeout bounds the original scientific process. Admission/postcondition work
    is streamed but is not a hard CPU, RSS or end-to-end wall-clock limit.
    """
    if (
        not attempt.is_absolute()
        or attempt.resolve() != attempt
        or not attempt.parent.is_dir()
    ):
        return {"status": "rejected", "reason_code": "attempt_path_invalid"}
    try:
        attempt.mkdir(mode=0o700)
    except FileExistsError:
        return {"status": "rejected", "reason_code": "attempt_exists"}
    private = attempt / "private"
    private.mkdir(mode=0o700)
    outcome = None
    try:
        if (
            type(threads) is not int
            or not 2 <= threads <= 8
            or type(mapq) is not int
            or not 0 <= mapq <= 255
            or isinstance(timeout, bool)
            or not isinstance(timeout, (float, int))
            or not math.isfinite(timeout)
            or timeout <= 0
        ):
            raise QualificationError("parameters_invalid")
        runtime = load_runtime_binding(runtime_binding)
        reference = load_reference_binding(reference_binding, reference_sha256)
        prepared = prepare_inputs(samples, runtime, reference, attempt)
        identity = implementation_identity()
        if (
            expected_input_identity is not None
            and prepared.input_identity != expected_input_identity
        ):
            raise QualificationError("input_binding_changed")
        if (
            expected_implementation_sha256 is not None
            and identity["sha256"] != expected_implementation_sha256
        ):
            raise QualificationError("implementation_binding_changed")
        identity.update(
            runtime_lock_sha256=runtime.lock_sha256,
            runtime_binding_sha256=runtime.binding_sha256,
            reference_binding_sha256=reference.binding_sha256,
        )
        write_exclusive(private / "identity.json", identity)
        write_exclusive(private / "inputs.json", prepared.input_identity)
        if cancelled():
            raise QualificationError("cancelled")
        output = attempt / "output"
        output.mkdir(mode=0o700)
        plan = plan_calls(
            output, prepared.reference_prefix, list(prepared.samples), threads, mapq
        )
        science_tools = {
            name: runtime.tools[name]
            for name in ("bowtie2", "samtools", "bamToBed", "gzip", "cLoops2", "rm")
        }
        shim_dir, config = prepare_shims(
            attempt, science_tools, plan, identity["sha256"]
        )
        config_sha = digest(config)
        argv = [
            str(runtime.python),
            "-I",
            "-B",
            str(runtime.script),
            "-fqd",
            str(prepared.staged_fastq_dir),
            "-o",
            str(output),
            "-ref",
            str(prepared.reference_prefix),
            "-n",
            "1",
            "-p",
            str(threads),
            "-mapq",
            str(mapq),
        ]
        # Inherit neither caller PATH nor Python injection; no HOME reassignment.
        env = {
            "PATH": f"{shim_dir}:{runtime.prefix}/bin:/usr/bin:/bin",
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONNOUSERSITE": "1",
            "LANG": "C.UTF-8",
            "LC_ALL": "C.UTF-8",
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "MPLCONFIGDIR": str(private / "mpl"),
            "XDG_CACHE_HOME": str(private / "cache"),
            "TMPDIR": str(private / "tmp"),
        }
        (private / "tmp").mkdir(mode=0o700)
        write_exclusive(
            private / "command.json",
            {"argv": argv, "cwd": str(attempt), "environment": env},
        )
        execution = execute(argv, attempt, env, private, timeout, cancelled)
        write_exclusive(private / "execution.json", asdict(execution))
        if execution.reason:
            raise QualificationError(execution.reason)
        if digest(config) != config_sha:
            raise QualificationError("call_plan_changed")
        # Prefer a precise child failure over the top-level parser consequence.
        verify_calls(config)
        if execution.returncode:
            raise QualificationError("upstream_failed")
        if cancelled():
            raise QualificationError("cancelled")
        counts = {k: v["raw_pairs"] for k, v in prepared.samples.items()}
        (private / "outputs-check").mkdir(mode=0o700)
        results = verify_outputs(
            output, counts, prepared.contigs, mapq, private / "outputs-check"
        )
        for token in prepared.samples:
            if cancelled():
                raise QualificationError("cancelled")
            directory = private / ("pair-" + token)
            directory.mkdir(mode=0o700)
            bam = output / token / (token + ".bam")
            command = [str(runtime.tools["samtools"]), "view", "-h", str(bam)]
            viewed = execute(command, attempt, env, directory, timeout, cancelled)
            write_exclusive(
                directory / "execution.json", asdict(viewed) | {"argv": command}
            )
            if viewed.reason or viewed.returncode:
                raise QualificationError(
                    viewed.reason or "bam_read_failed", token, "pairing"
                )
            try:
                with (directory / "upstream.stdout").open() as sam:
                    verify_pairs(
                        sam,
                        iter_bedpe(
                            output / token / (token + "_all.bedpe.gz"), prepared.contigs
                        ),
                        directory / "origins.sqlite",
                    )
            except PairError as error:
                raise QualificationError(
                    error.reason_code, token, "pairing", "all"
                ) from None
        if cancelled():
            raise QualificationError("cancelled")
        # Recheck pinned bytes and staged inputs before committing a completion fact.
        if load_runtime_binding(runtime_binding) != runtime:
            raise QualificationError("runtime_changed")
        if load_reference_binding(reference_binding, reference_sha256) != reference:
            raise QualificationError("reference_changed")
        if implementation_identity()["sha256"] != identity["sha256"]:
            raise QualificationError("qualification_changed")
        for token, sample in prepared.samples.items():
            for mate in ("r1", "r2"):
                if (
                    digest(
                        prepared.staged_fastq_dir
                        / (token + "_" + mate.upper() + ".fastq.gz")
                    )
                    != sample["sha256"][mate]
                ):
                    raise QualificationError("staged_input_changed", token)
        if cancelled():
            raise QualificationError("cancelled")
        complete = {
            "schema_version": "hitrac-qualification-complete-v1",
            "identity": identity,
            "results": results,
            "mapq": mapq,
            "threads": threads,
        }
        temporary = private / "completion.pending.json"
        write_exclusive(temporary, complete)
        outcome = {
            "status": "complete",
            "reason_code": None,
            "qualification_sha256": identity["sha256"],
            "samples": list(prepared.samples),
        }
        # All diagnostic writes precede the atomic completion fact.
        try:
            write_exclusive(private / "outcome.json", outcome)
        except OSError:
            raise QualificationError("diagnostic_write_failed") from None
        if cancelled():
            raise QualificationError("cancelled")
        os.rename(temporary, attempt / "complete.json")
        return outcome
    except (QualificationError, AdmissionError, CallFailure, OutputRejected) as error:
        code = getattr(
            error, "reason_code", getattr(error, "code", "qualification_rejected")
        )
        outcome = {"status": "rejected", "reason_code": code}
        for name in ("sample", "stage", "collection"):
            value = getattr(error, name, None)
            if value is not None:
                outcome[name] = value
    except Exception:
        # Never expose exception text, argv or paths through the CLI result.
        try:
            (private / "qualification-exception.txt").write_text(traceback.format_exc())
        except OSError:
            return {"status": "rejected", "reason_code": "diagnostic_write_failed"}
        outcome = {"status": "rejected", "reason_code": "qualification_internal_error"}
    try:
        target = private / "outcome.json"
        if target.exists():
            # Rename failure after a prepared outcome must not retain a success
            # diagnosis; completion still requires complete.json.
            pending = private / "rejection.pending.json"
            write_exclusive(pending, outcome)
            os.replace(pending, target)
        else:
            write_exclusive(target, outcome)
    except OSError:
        return {"status": "rejected", "reason_code": "diagnostic_write_failed"}
    return outcome


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Private fixed tracPre2 qualification; not public execution."
    )
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--attempt", type=Path, required=True)
    args = parser.parse_args()
    cancelled = threading.Event()
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda *_: cancelled.set())
    try:
        request = strict_json(args.request.read_text())
        result = qualify(
            runtime_binding=Path(request["runtime_binding"]),
            reference_binding=Path(request["reference_binding"]),
            reference_sha256=request["reference_sha256"],
            samples=[
                Sample(s["id"], Path(s["r1"]), Path(s["r2"]))
                for s in request["samples"]
            ],
            attempt=args.attempt,
            threads=request.get("threads", 2),
            mapq=request.get("mapq", 10),
            timeout=request.get("timeout", 300),
            cancelled=cancelled.is_set,
        )
    except (OSError, ValueError, TypeError, KeyError):
        result = {"status": "rejected", "reason_code": "request_invalid"}
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    sys.exit(main())
