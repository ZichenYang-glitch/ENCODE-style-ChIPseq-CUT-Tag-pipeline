"""Fixed call plan, real transparent child process, and fail-closed receipts."""

import concurrent.futures
import json
from pathlib import Path
import subprocess

import pytest

from encode_pipeline.adapters.hitrac_preprocess.calls import (
    CallFailure,
    digest,
    plan_calls,
    prepare_shims,
    verify_calls,
)


def setup(tmp_path, argv=None, body="import sys;sys.exit(0)"):
    attempt = tmp_path / "attempt"
    attempt.mkdir()
    (attempt / "private").mkdir()
    tool = tmp_path / "tool"
    import sys

    tool.write_text("#!" + sys.executable + "\n" + body + "\n")
    tool.chmod(0o700)
    entry = dict(tool="samtools", stage="view", sample="s000001", argv=argv or ["view"])
    shim, config = prepare_shims(
        attempt, {"samtools": tool}, [entry], "pinned-test-identity"
    )
    return attempt, shim, config


def run(attempt, shim, argv=("view",)):
    return subprocess.run(
        [str(shim / "samtools"), *argv], cwd=attempt, capture_output=True, text=True
    )


def test_transparent_argv_cwd_and_streams(tmp_path):
    attempt, shim, config = setup(
        tmp_path,
        ["view", "space value"],
        "import os,sys;print(repr(sys.argv[1:]));print(os.getcwd(),file=sys.stderr)",
    )
    done = run(attempt, shim, ("view", "space value"))
    assert done.returncode == 0
    assert done.stdout == "['view', 'space value']\n"
    assert done.stderr == str(attempt) + "\n"
    receipts = verify_calls(config)
    assert receipts[0]["argv"] == ["view", "space value"]


@pytest.mark.parametrize("exitcode", [1, 73, 127])
def test_nonzero_is_preserved_and_rejected(tmp_path, exitcode):
    attempt, shim, config = setup(tmp_path, body=f"import sys;sys.exit({exitcode})")
    assert run(attempt, shim).returncode == exitcode
    with pytest.raises(CallFailure) as caught:
        verify_calls(config)
    assert (caught.value.code, caught.value.sample, caught.value.stage) == (
        "child_failed",
        "s000001",
        "view",
    )


def test_child_signal_is_preserved(tmp_path):
    attempt, shim, config = setup(
        tmp_path, body="import os,signal;os.kill(os.getpid(),signal.SIGTERM)"
    )
    assert run(attempt, shim).returncode == -15
    with pytest.raises(CallFailure, match="child_failed"):
        verify_calls(config)


@pytest.mark.parametrize(
    "damage", ["missing_end", "malformed", "extra", "wrong_sample", "duplicate_key"]
)
def test_damaged_journal(tmp_path, damage):
    attempt, shim, config = setup(tmp_path)
    assert run(attempt, shim).returncode == 0
    directory = attempt / "private/calls"
    start = next(directory.glob("*.start.json"))
    if damage == "missing_end":
        next(directory.glob("*.end.json")).unlink()
    elif damage == "malformed":
        start.write_text("{")
    elif damage == "extra":
        (directory / "unexpected").write_text("{}")
    elif damage == "wrong_sample":
        d = json.loads(start.read_text())
        d["sample"] = "s000002"
        start.write_text(json.dumps(d))
    else:
        start.write_text('{"x":1,"x":2}')
    with pytest.raises(CallFailure):
        verify_calls(config)


def test_unrecognized_call_is_not_executed(tmp_path):
    marker = tmp_path / "must-not-exist"
    attempt, shim, config = setup(
        tmp_path, body=f"from pathlib import Path;Path({str(marker)!r}).touch()"
    )
    assert run(attempt, shim, ("unexpected",)).returncode == 126
    assert not marker.exists()
    with pytest.raises(CallFailure, match="call_unrecognized"):
        verify_calls(config)


def test_concurrent_unique_records_reject_duplicate_stage(tmp_path):
    attempt, shim, config = setup(tmp_path)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda _: run(attempt, shim), range(4)))
    assert all(r.returncode == 0 for r in results)
    assert len(list((attempt / "private/calls").glob("*.json"))) == 8
    with pytest.raises(CallFailure, match="call_duplicate"):
        verify_calls(config)


def test_journal_write_failure_leaves_missing_call(tmp_path):
    attempt, shim, config = setup(tmp_path)
    directory = attempt / "private/calls"
    directory.rmdir()
    directory.write_text("not a directory")
    assert run(attempt, shim).returncode == 125
    with pytest.raises(CallFailure, match="call_records_invalid"):
        verify_calls(config)


def test_redirection_prevents_tool_and_plan_detects_missing(tmp_path):
    attempt, shim, config = setup(tmp_path)
    out = attempt / "nonexistent" / "result"
    import shlex

    cmd = shlex.join([str(shim / "samtools"), "view"]) + " > " + shlex.quote(str(out))
    result = subprocess.run(["/bin/sh", "-c", cmd], cwd=attempt, capture_output=True)
    assert result.returncode != 0
    assert not list((attempt / "private/calls").iterdir())
    with pytest.raises(CallFailure, match="call_missing"):
        verify_calls(config)


def test_plan_includes_only_original_cleanup_and_dynamic_mapq(tmp_path):
    plan = plan_calls(tmp_path / "out", tmp_path / "ref", ["s000001", "s000002"], 2, 17)
    assert len(plan) == 15
    assert [p["stage"] for p in plan].count("remove_sam") == 2
    assert [p["stage"] for p in plan].count("remove_qc") == 1
    view = next(p for p in plan if p["stage"] == "view")
    assert view["argv"][view["argv"].index("-q") + 1] == "17"
    removed = [arg for p in plan if p["tool"] == "rm" for arg in p["argv"]]
    assert all(p.endswith((".sam", "_bedpeQc.txt")) for p in removed)


def test_config_and_tool_drift_refused(tmp_path):
    attempt, shim, config = setup(tmp_path)
    doc = json.loads(config.read_text())
    Path(doc["tools"]["samtools"]["path"]).write_text("drift")
    assert run(attempt, shim).returncode == 125
    with pytest.raises(CallFailure, match="child_failed"):
        verify_calls(config)


@pytest.mark.parametrize(
    "field,kind,value",
    [
        ("pid", "start", None),
        ("pid", "start", True),
        ("pid", "start", 0),
        ("time_ns", "start", None),
        ("time_ns", "end", None),
        ("time_ns", "end", False),
        ("unexpected", "end", "unplanned"),
    ],
)
def test_incomplete_or_malformed_receipt_metadata_rejected(
    tmp_path, field, kind, value
):
    attempt, shim, config = setup(tmp_path)
    assert run(attempt, shim).returncode == 0
    directory = attempt / "private/calls"
    path = next(directory.glob(f"*.{kind}.json"))
    data = json.loads(path.read_text())
    if value is None:
        data.pop(field)
    else:
        data[field] = value
    path.write_text(json.dumps(data))
    if kind == "start":
        # Keep the start digest consistent to isolate the metadata schema check.
        end = next(directory.glob("*.end.json"))
        entry = json.loads(end.read_text())
        entry["start_sha256"] = digest(path)
        end.write_text(json.dumps(entry))
    with pytest.raises(CallFailure, match="call_records_invalid"):
        verify_calls(config)
