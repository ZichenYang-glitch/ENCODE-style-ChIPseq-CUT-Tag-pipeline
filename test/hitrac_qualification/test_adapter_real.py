"""H3 private adapter/runner qualification, never public execution or publication."""

from collections import Counter
import json
from pathlib import Path

import pytest

from encode_pipeline.adapters.hitrac_preprocess.adapter import HiTracPreprocessAdapter
from encode_pipeline.adapters.hitrac_preprocess.admission import sha256_file
from encode_pipeline.adapters.hitrac_preprocess.execution import admit_runtime
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.defaults import create_default_process_runner
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import WorkspacePlanner
from encode_pipeline.workers.settings import WorkerSettings
from test_qualification_real import (
    coordinates as shared_coordinates,
    expected_rows,
    rows,
)

coordinates = shared_coordinates

pytestmark = pytest.mark.real_execution


def composed(coordinates, tmp_path, scenario="positive", *, timeout=120):
    directory = Path(coordinates["TINY_INPUTS"]) / scenario / "fastq"
    samples = [
        {
            "sample_id": f"sample {i}",
            "fastq_1": str(p),
            "fastq_2": str(p.with_name(p.name.replace("_R1", "_R2"))),
        }
        for i, p in enumerate(sorted(directory.glob("*_R1.fastq.gz")), 1)
    ]
    inputs = WorkflowInputs(
        config={}, samples=samples, options={"threads": 2, "mapq": 10}
    )
    rb = Path(coordinates["RUNTIME_BINDING"])
    runtime = admit_runtime(rb, sha256_file(rb), timeout=timeout)
    adapter = HiTracPreprocessAdapter(runtime=runtime)
    result = adapter.bind_reference_profile(
        inputs,
        {
            "schema_version": "hitrac-reference-profile-v1",
            "binding": coordinates["REFERENCE_BINDING"],
            "sha256": coordinates["REFERENCE_SHA256"],
        },
    )
    assert result.is_success, result.issues
    adapter = result.value.adapter
    assert adapter.execution_availability().execution == "not_configured"
    registry = WorkflowRegistry([adapter])
    workspace = tmp_path / "workspace"
    plan = ExecutionPlan(
        plan_id="qualification",
        run_id="private",
        workflow_id=adapter.metadata.workflow_id,
        status=PlanStatus.PENDING,
        inputs_snapshot=inputs.to_dict(),
    )
    planned = WorkspacePlanner(registry).plan_workspace(plan, workspace)
    assert planned.is_success, planned.issues
    built = CommandBuilder(registry).build_command(planned.value, workspace)
    assert built.is_success, built.issues
    materialized = WorkspaceMaterializer().materialize(
        built.value.workspace_plan, workspace
    )
    assert materialized.is_success, materialized.issues
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'unused.db'}",
        redis_url="redis://127.0.0.1:1/0",
        queue_name="unused",
        workspace_root=tmp_path,
        job_timeout_seconds=180,
    )
    runner = create_default_process_runner(registry=registry, settings=settings)
    assert runner._managed_container_cleaner is None
    return adapter, runner, built.value.command_spec, workspace, inputs


def record(tmp_path, result, spec):
    (tmp_path / "runner-result.json").write_text(
        json.dumps(
            {
                "argv": spec.argv,
                "cwd": spec.cwd,
                "issues": [i.to_dict() for i in result.issues],
                "result": None if result.value is None else vars(result.value),
            },
            indent=2,
        )
    )


def test_real_adapter_workspace_command_and_runner_double_sample(coordinates, tmp_path):
    adapter, runner, spec, workspace, inputs = composed(coordinates, tmp_path)
    before = {
        p: sha256_file(Path(p))
        for row in inputs.samples
        for k, p in row.items()
        if k.startswith("fastq")
    }
    result = runner.run(spec)
    record(tmp_path, result, spec)
    assert result.is_success, result.issues
    assert result.value.exit_code == 0, result.value
    assert json.loads(result.value.stdout)["status"] == "complete"
    complete = json.loads((workspace / "hitrac-attempt/complete.json").read_text())
    expected = expected_rows(coordinates["design"])
    wanted = coordinates["design"]["scenarios"]["positive"][
        "conditional_expected_summary"
    ]
    assert set(complete["results"]["samples"]) == {"s000001", "s000002"}
    for token, sample in complete["results"]["samples"].items():
        output = workspace / "hitrac-attempt/output" / token
        assert sample["all"] == 8 and sample["noBg"] == 5
        assert sample["metrics"] == pytest.approx(wanted, rel=0, abs=1e-12)
        assert Counter(map(tuple, rows(output / f"{token}_all.bedpe.gz"))) == Counter(
            expected
        )
        assert Counter(
            map(tuple, rows(output / f"{token}_unique.bedpe.gz"))
        ) == Counter(expected[i] for i in (0, 3, 4, 5, 6))
        assert (output / f"{token}.bam").is_file()
    assert {p: sha256_file(Path(p)) for p in before} == before
    assert adapter.execution_availability().execution == "not_configured"
    repeat = runner.run(spec)
    assert repeat.is_success and repeat.value.exit_code == 1
    assert json.loads(repeat.value.stdout)["reason_code"] == "execution_binding_invalid"
    assert (
        json.loads((workspace / "hitrac-attempt/complete.json").read_text()) == complete
    )


@pytest.mark.parametrize(
    "scenario,reason",
    [
        ("singlemate_only", "pair_support_missing"),
        ("noBg_empty", "empty_pet_set"),
        ("all_trans", "no_cis_denominator"),
        ("multisample_one_empty", "empty_pet_set"),
    ],
)
def test_real_adapter_retains_h2_refusals(coordinates, tmp_path, scenario, reason):
    _, runner, spec, workspace, _ = composed(coordinates, tmp_path, scenario)
    result = runner.run(spec)
    record(tmp_path, result, spec)
    assert result.is_success and result.value.exit_code == 1, result
    assert json.loads(result.value.stdout)["reason_code"] == reason
    assert not (workspace / "hitrac-attempt/complete.json").exists()
    assert (workspace / "hitrac-attempt/output/tracPre_summary.txt").is_file()
