"""ENCODE command equivalence and fail-closed private execution configuration."""

from dataclasses import replace
import json
from pathlib import Path

import pytest

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.adapters.encode_execution import (
    EXECUTION_CONFIG_PATH,
    EncodeExecutionBinding,
)
from encode_pipeline.platform.adapters import CommandSpec, WorkflowInputs
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.command_builder import CommandBuilder


ROOT = Path(__file__).resolve().parents[2]
# Captured by running the pre-PR-5b CommandBuilder before modifying it.
BEFORE_COMMANDS = json.loads(
    (ROOT / "test/fixtures/encode-command-before.json").read_text()
)


def _inputs(tmp_path, options):
    return WorkflowInputs(
        config={"threads": 8, "use_control": False},
        samples=[
            {
                "sample": "S1",
                "fastq_1": str(tmp_path / "reads.fastq.gz"),
                "layout": "SE",
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
                "genome": "hs",
                "bowtie2_index": str(tmp_path / "index"),
            }
        ],
        options=options,
    )


def _runtime(root):
    for name in (
        "runner/bin/snakemake",
        "runner/bin/conda",
        "runner/libexec/micromamba",
        "mamba-root/bin/activate",
    ):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"#!/bin/sh\n")
        path.chmod(0o555)
    (root / "conda-envs").mkdir()
    for path in [root, *[p for p in root.rglob("*") if p.is_dir()]]:
        path.chmod(0o755)
    return EncodeExecutionBinding(
        project_root=ROOT,
        snakemake_executable=root / "runner/bin/snakemake",
        conda_prefix=root / "conda-envs",
    )


@pytest.mark.parametrize("case", BEFORE_COMMANDS)
def test_adapter_and_service_match_pre_migration_commands_byte_for_byte(tmp_path, case):
    runtime = tmp_path / "runtime"
    binding = (
        _runtime(runtime)
        if case["managed_runtime"]
        else EncodeExecutionBinding(project_root=ROOT)
    )
    adapter = EncodeStyleWorkflowAdapter(execution=binding)
    workspace = tmp_path / "workspace space α"
    inputs = _inputs(tmp_path, case["options"])
    planned = adapter.plan_workspace(inputs, workspace)
    assert planned.is_success, planned.issues
    result = adapter.build_command(planned.value, workspace)
    assert result.is_success, result.issues
    assert isinstance(result.value, CommandSpec)
    expected = json.dumps(case["command"], ensure_ascii=False)
    for marker, path in [
        ("{workspace}", workspace),
        ("{runtime}", runtime),
        ("{project_root}", ROOT),
    ]:
        expected = expected.replace(marker, str(path))
    expected = json.loads(expected)
    command = result.value
    assert [s.encode() for s in command.argv] == [s.encode() for s in expected["argv"]]
    assert {k.encode(): v.encode() for k, v in command.env.items()} == {
        k.encode(): v.encode() for k, v in expected["env"].items()
    }
    assert [s.encode() for s in command.preflight_argv] == [
        s.encode() for s in expected["preflight_argv"]
    ]
    assert command.to_dict() == expected
    assert command.cwd is None
    assert not workspace.exists()
    registry = WorkflowRegistry([adapter], legacy_execution_fallbacks=(adapter,))
    delegated = CommandBuilder(registry).build_command(
        ExecutionPlan(
            plan_id="p",
            run_id="r",
            workflow_id=adapter.metadata.workflow_id,
            status=PlanStatus.PENDING,
            inputs_snapshot=inputs.to_dict(),
            workspace_plan=planned.value,
        ),
        workspace,
    )
    assert delegated.is_success, delegated.issues
    assert delegated.value.command_spec == command


def _plan(tmp_path):
    adapter = EncodeStyleWorkflowAdapter()
    result = adapter.plan_workspace(
        _inputs(tmp_path, {"cores": 8}), tmp_path / "workspace"
    )
    assert result.is_success
    return adapter, result.value


def test_workspace_file_contract_includes_canonical_execution_json(tmp_path):
    _adapter, plan = _plan(tmp_path)
    assert [path for path, _ in plan.files] == [
        "config/config.yaml",
        "config/samples.tsv",
        EXECUTION_CONFIG_PATH,
    ]
    raw = dict(plan.files)[EXECUTION_CONFIG_PATH]
    value = json.loads(raw)
    assert set(value) == {"schema_version", "cores", "workspace_contract_sha256"}
    assert value["schema_version"] == "1.0.0"
    assert value["cores"] == 8
    assert len(value["workspace_contract_sha256"]) == 64
    assert (
        raw
        == (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()
    )


@pytest.mark.parametrize(
    "case",
    [
        "missing",
        "extra-key",
        "missing-key",
        "duplicate-key",
        "version",
        "hash",
        "valid-cores-changed",
        "whitespace",
        "invalid-json",
        "oversized",
        "non-object",
        "other-file-changed",
        "directories-changed",
    ],
)
def test_private_execution_configuration_tampering_is_rejected(tmp_path, case):
    adapter, plan = _plan(tmp_path)
    raw = dict(plan.files)[EXECUTION_CONFIG_PATH]
    value = json.loads(raw)
    if case == "extra-key":
        value["extra"] = "private"
    elif case == "missing-key":
        value.pop("cores")
    elif case == "version":
        value["schema_version"] = "2.0.0"
    elif case == "hash":
        value["workspace_contract_sha256"] = "0" * 64
    elif case == "valid-cores-changed":
        value["cores"] = 9
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()
    if case == "duplicate-key":
        raw = raw.replace(b'"cores":8', b'"cores":8,"cores":8')
    elif case == "whitespace":
        raw += b"\n"
    elif case == "invalid-json":
        raw = b"{"
    elif case == "oversized":
        raw = b" " * 1025
    elif case == "non-object":
        raw = b"[]"
    files = tuple(
        (path, raw if path == EXECUTION_CONFIG_PATH else content)
        for path, content in plan.files
        if case != "missing" or path != EXECUTION_CONFIG_PATH
    )
    if case == "other-file-changed":
        files = tuple(
            (path, content + b"\n" if path == "config/config.yaml" else content)
            for path, content in files
        )
    altered = replace(plan, files=files)
    if case == "directories-changed":
        altered = replace(altered, directories=("other",))
    result = adapter.build_command(altered, tmp_path / "workspace")
    assert result.is_failure
    assert result.issues[0].code == "ENCODE_EXECUTION_CONFIG_INVALID"
    assert str(tmp_path) not in repr(result.issues)


@pytest.mark.parametrize("cores", [None, True, False, "8", 8.0, [], 0, -1, 1025])
def test_private_execution_cores_type_and_range_are_revalidated(tmp_path, cores):
    adapter, plan = _plan(tmp_path)
    value = json.loads(dict(plan.files)[EXECUTION_CONFIG_PATH])
    value["cores"] = cores
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()
    plan = replace(
        plan,
        files=tuple(
            (path, raw if path == EXECUTION_CONFIG_PATH else data)
            for path, data in plan.files
        ),
    )
    result = adapter.build_command(plan, tmp_path / "workspace")
    assert result.is_failure
    assert result.issues[0].code == "COMMAND_BUILD_INVALID_CORES"


@pytest.mark.parametrize("case", ["match", "changed", "symlink", "directory"])
def test_materialized_execution_file_must_match_planned_bytes(tmp_path, case):
    adapter, plan = _plan(tmp_path)
    workspace = tmp_path / "workspace"
    target = workspace / EXECUTION_CONFIG_PATH
    target.parent.mkdir(parents=True)
    raw = dict(plan.files)[EXECUTION_CONFIG_PATH]
    if case == "symlink":
        source = workspace / "source.json"
        source.write_bytes(raw)
        target.symlink_to(source)
    elif case == "directory":
        target.mkdir()
    else:
        target.write_bytes(
            raw if case == "match" else raw.replace(b'"cores":8', b'"cores":9')
        )
    result = adapter.build_command(plan, workspace)
    if case == "match":
        assert result.is_success, result.issues
    else:
        assert result.is_failure
        assert result.issues[0].code == "ENCODE_EXECUTION_CONFIG_INVALID"
