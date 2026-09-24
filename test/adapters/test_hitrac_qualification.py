"""Private qualification-entry contracts before any scientific execution."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
from threading import Event
from types import SimpleNamespace

import pytest

from encode_pipeline.adapters.hitrac_preprocess.qualification import qualify
from encode_pipeline.adapters.hitrac_preprocess import qualification as entry
from encode_pipeline.adapters.hitrac_preprocess.process import Execution


def _arguments(tmp_path):
    return {
        "runtime_binding": tmp_path / "runtime.json",
        "reference_binding": tmp_path / "reference.json",
        "reference_sha256": "1" * 64,
        "samples": [],
        "attempt": tmp_path / "attempt",
    }


def test_previous_attempt_is_rejected_without_touching_its_evidence(tmp_path):
    arguments = _arguments(tmp_path)
    attempt = arguments["attempt"]
    attempt.mkdir()
    evidence = attempt / "private-diagnostic"
    evidence.write_bytes(b"prior outputs and failure evidence")
    before = hashlib.sha256(evidence.read_bytes()).hexdigest()
    assert qualify(**arguments) == {
        "status": "rejected",
        "reason_code": "attempt_exists",
    }
    assert list(attempt.iterdir()) == [evidence]
    assert hashlib.sha256(evidence.read_bytes()).hexdigest() == before


@pytest.mark.parametrize(
    "change",
    [
        {"threads": 1},
        {"threads": True},
        {"mapq": "10"},
        {"mapq": -1},
        {"timeout": 0},
        {"timeout": float("nan")},
    ],
)
def test_invalid_execution_parameters_fail_before_tool_or_binding_access(
    tmp_path, change
):
    arguments = _arguments(tmp_path) | change
    result = qualify(**arguments)
    assert result == {"status": "rejected", "reason_code": "parameters_invalid"}
    attempt = arguments["attempt"]
    assert not (attempt / "complete.json").exists()
    assert not (attempt / "output").exists()
    assert json.loads((attempt / "private/outcome.json").read_text()) == result
    assert (attempt.stat().st_mode & 0o777) == 0o700


def test_symlink_attempt_does_not_create_output_at_target(tmp_path):
    arguments = _arguments(tmp_path)
    target = tmp_path / "target"
    target.mkdir()
    arguments["attempt"].symlink_to(target, target_is_directory=True)
    assert qualify(**arguments) == {
        "status": "rejected",
        "reason_code": "attempt_path_invalid",
    }
    assert list(target.iterdir()) == []


@pytest.mark.parametrize(
    "request_text",
    [
        '{"runtime_binding":"private-payload-marker", "samples": NaN}',
        '{"runtime_binding":"private-payload-marker", "runtime_binding":"duplicate"}',
        '{"samples":[{"id":"private-payload-marker"}]}',
    ],
)
def test_original_private_cli_bad_request_is_controlled_and_path_free(
    tmp_path, request_text
):
    source = Path(__file__).resolve().parents[2]
    path = tmp_path / "private-path-marker.json"
    path.write_text(request_text)
    attempt = tmp_path / "attempt"
    command = [
        sys.executable,
        "-I",
        "-S",
        str(source / "scripts/qualify_hitrac_preprocess.py"),
        "--request",
        str(path),
        "--attempt",
        str(attempt),
    ]
    process = subprocess.run(
        command, cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert process.returncode == 1
    assert json.loads(process.stdout) == {
        "status": "rejected",
        "reason_code": "request_invalid",
    }
    assert process.stderr == ""
    assert "private-payload-marker" not in process.stdout
    assert "private-path-marker" not in process.stdout
    assert str(tmp_path) not in process.stdout
    assert not attempt.exists()


def _completed_science_stubs(tmp_path, monkeypatch):
    """Isolate final persistence; scientific prerequisites are explicit stubs.

    This verifies the real qualify control flow and atomic completion boundary,
    not a scientific execution, child plan or runtime-admission success.
    """
    arguments = _arguments(tmp_path)
    staged = tmp_path / "staged"
    staged.mkdir()
    for mate in ("R1", "R2"):
        (staged / f"s000001_{mate}.fastq.gz").write_bytes(b"synthetic-staged-input")
    tool = Path(sys.executable)
    runtime = SimpleNamespace(
        prefix=tmp_path,
        python=tool,
        script=tmp_path / "original-script",
        tools={
            name: tool
            for name in ("bowtie2", "samtools", "bamToBed", "gzip", "cLoops2", "rm")
        },
        lock_sha256="b" * 64,
        binding_sha256="c" * 64,
    )
    reference = SimpleNamespace(binding_sha256="d" * 64)
    prepared = SimpleNamespace(
        staged_fastq_dir=staged,
        reference_prefix=tmp_path / "reference",
        samples={
            "s000001": {
                "raw_pairs": 1,
                "sha256": {
                    mate: entry.digest(staged / f"s000001_{mate.upper()}.fastq.gz")
                    for mate in ("r1", "r2")
                },
            }
        },
        contigs={"chrA": 100},
        input_identity={"test": "isolated completion persistence"},
    )
    monkeypatch.setattr(entry, "load_runtime_binding", lambda _: runtime)
    monkeypatch.setattr(entry, "load_reference_binding", lambda *_: reference)
    monkeypatch.setattr(entry, "prepare_inputs", lambda *_: prepared)
    monkeypatch.setattr(
        entry, "implementation_identity", lambda: {"sha256": "a" * 64, "files": {}}
    )
    monkeypatch.setattr(entry, "verify_calls", lambda _: [])
    monkeypatch.setattr(entry, "verify_pairs", lambda *_: None)

    def simulated_execution(argv, cwd, env, private, timeout, cancelled):
        (private / "upstream.stdout").write_text("synthetic, not scientific evidence\n")
        (private / "upstream.stderr").write_text("")
        return Execution(0, None, 12345, 0)

    def accepted_outputs(output, *_):
        sample = output / "s000001"
        sample.mkdir()
        (sample / "retained-private-output").write_bytes(
            b"retain on persistence failure"
        )
        return {"samples": {"s000001": {"all": 1, "noBg": 1}}}

    monkeypatch.setattr(entry, "execute", simulated_execution)
    monkeypatch.setattr(entry, "verify_outputs", accepted_outputs)
    return arguments


def test_outcome_write_failure_cannot_leave_a_complete_marker(tmp_path, monkeypatch):
    arguments = _completed_science_stubs(tmp_path, monkeypatch)
    original_write = entry.write_exclusive

    def fail_outcome(path, payload):
        if path.name == "outcome.json":
            raise PermissionError("synthetic-private-marker")
        return original_write(path, payload)

    monkeypatch.setattr(entry, "write_exclusive", fail_outcome)
    result = qualify(**arguments)
    assert result == {"status": "rejected", "reason_code": "diagnostic_write_failed"}
    assert not (arguments["attempt"] / "complete.json").exists()
    assert (
        arguments["attempt"] / "output/s000001/retained-private-output"
    ).read_bytes() == b"retain on persistence failure"


def test_cancel_during_final_identity_recheck_refuses_completion(tmp_path, monkeypatch):
    arguments = _completed_science_stubs(tmp_path, monkeypatch)
    cancelled = Event()
    original_load = entry.load_runtime_binding
    calls = []

    def cancel_during_second_check(binding):
        calls.append(binding)
        if len(calls) == 2:
            cancelled.set()
        return original_load(binding)

    monkeypatch.setattr(entry, "load_runtime_binding", cancel_during_second_check)
    result = qualify(**arguments, cancelled=cancelled.is_set)
    assert len(calls) == 2 and cancelled.is_set()
    assert result == {"status": "rejected", "reason_code": "cancelled"}
    assert not (arguments["attempt"] / "complete.json").exists()
    assert (
        arguments["attempt"] / "output/s000001/retained-private-output"
    ).read_bytes() == b"retain on persistence failure"
