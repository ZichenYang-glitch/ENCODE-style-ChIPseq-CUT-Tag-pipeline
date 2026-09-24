"""Real current-identity Hi-TrAC result consumers; no API/RQ product chain."""

from collections import Counter
from hashlib import sha256
import json
import gzip
from pathlib import Path

import pytest

from encode_pipeline.adapters.hitrac_preprocess.admission import sha256_file
from encode_pipeline.adapters.hitrac_preprocess.execution import admit_runtime
from encode_pipeline.adapters.hitrac_preprocess.results import (
    HiTracPreprocessResultsAdapter,
)
from encode_pipeline.platform.adapters import (
    WorkflowInputs,
    QcSourceArtifact,
    QcSourceDocument,
)
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


def compose_results(coordinates, tmp_path, *, scenario="positive", mapq=10):
    directory = Path(coordinates["TINY_INPUTS"]) / scenario / "fastq"
    inputs = WorkflowInputs(
        config={},
        samples=[
            {
                "sample_id": f"sample {index}",
                "fastq_1": str(path),
                "fastq_2": str(path.with_name(path.name.replace("_R1", "_R2"))),
            }
            for index, path in enumerate(sorted(directory.glob("*_R1.fastq.gz")), 1)
        ],
        options={"threads": 2, "mapq": mapq},
    )
    runtime = Path(coordinates["RUNTIME_BINDING"])
    base = HiTracPreprocessResultsAdapter(
        runtime=admit_runtime(runtime, sha256_file(runtime), timeout=120)
    )
    bound = base.bind_reference_profile(
        inputs,
        {
            "schema_version": "hitrac-reference-profile-v1",
            "binding": coordinates["REFERENCE_BINDING"],
            "sha256": coordinates["REFERENCE_SHA256"],
        },
    )
    assert bound.is_success, bound.issues
    adapter = bound.value.adapter
    assert adapter.execution_availability().execution == "available"
    registry = WorkflowRegistry([adapter])
    workspace = tmp_path / "workspace"
    plan = ExecutionPlan(
        plan_id="results-qualification",
        run_id="private-results",
        workflow_id=adapter.metadata.workflow_id,
        status=PlanStatus.PENDING,
        inputs_snapshot=inputs.to_dict(),
    )
    planned = WorkspacePlanner(registry).plan_workspace(plan, workspace)
    assert planned.is_success, planned.issues
    built = CommandBuilder(registry).build_command(planned.value, workspace)
    assert built.is_success, built.issues
    assert (
        WorkspaceMaterializer()
        .materialize(built.value.workspace_plan, workspace)
        .is_success
    )
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'unused.db'}",
        redis_url="redis://127.0.0.1:1/0",
        queue_name="unused",
        workspace_root=tmp_path,
        job_timeout_seconds=180,
    )
    runner = create_default_process_runner(registry=registry, settings=settings)
    return adapter, inputs, workspace, runner, built.value.command_spec


def source_document(workspace, artifacts):
    item = next(value for value in artifacts if value.output_type == "hitrac_summary")
    return QcSourceDocument(
        QcSourceArtifact(
            artifact_id="summary-source",
            output_type=item.output_type,
            relative_path=item.relative_path,
            metadata={"scope": "run"},
        ),
        (workspace / item.relative_path).read_bytes(),
    )


@pytest.mark.parametrize("mapq", [10, 30])
def test_real_results_double_sample_original_bytes_all_metrics_and_mapping(
    coordinates, tmp_path, mapq
):
    adapter, inputs, workspace, runner, spec = compose_results(
        coordinates, tmp_path, mapq=mapq
    )
    before = {
        row[key]: sha256_file(Path(row[key]))
        for row in inputs.samples
        for key in ("fastq_1", "fastq_2")
    }
    result = runner.run(spec)
    (tmp_path / "runner.json").write_text(
        json.dumps(
            {
                "result": None if result.value is None else vars(result.value),
                "issues": [i.to_dict() for i in result.issues],
            }
        )
    )
    assert result.is_success and result.value.exit_code == 0, result
    extracted = adapter.extract_artifacts(inputs, workspace)
    assert extracted.is_success, extracted.issues
    assert len(extracted.value) == 5
    assert {item.output_type for item in extracted.value} == {
        "hitrac_summary",
        "hitrac_bedpe_all",
        "hitrac_bedpe_no_bg",
    }
    expected = expected_rows(coordinates["design"])
    for item in extracted.value:
        projected = workspace / item.relative_path
        original = (
            workspace
            / "hitrac-attempt/output"
            / Path(item.relative_path).relative_to(
                Path(*Path(item.relative_path).parts[:3])
            )
        )
        assert projected.read_bytes() == original.read_bytes()
        assert sha256(projected.read_bytes()).hexdigest() == item.metadata["sha256"]
        if item.output_type == "hitrac_summary":
            continue
        assert item.metadata["sample_id"] in {"sample 1", "sample 2"}
        wanted = (
            expected
            if item.output_type == "hitrac_bedpe_all"
            else [expected[i] for i in (0, 3, 4, 5, 6)]
        )
        assert Counter(map(tuple, rows(projected))) == Counter(wanted)
    metrics = adapter.extract_qc_metrics(
        inputs, (source_document(workspace, extracted.value),)
    )
    assert metrics.is_success, metrics.issues
    assert len(metrics.value) == 30
    expected_values = coordinates["design"]["scenarios"]["positive"][
        "conditional_expected_summary"
    ]
    for sample in ("sample 1", "sample 2"):
        values = [m for m in metrics.value if m.sample_id == sample]
        assert [float(m.value) for m in values] == pytest.approx(
            expected_values, rel=0, abs=1e-12
        )
        assert values[3].value == 8 and values[9].value == 5
        assert values[3].display_name == f"total mapped PETs (mapq>={mapq})"
        assert all(m.qc_flag is None for m in values)
    repeat = adapter.extract_artifacts(inputs, workspace)
    assert repeat.is_success and repeat.value == extracted.value
    assert {path: sha256_file(Path(path)) for path in before} == before
    assert all(
        not any(
            private in item.relative_path
            for private in (".bam", ".fastq", "private", ".log")
        )
        for item in extracted.value
    )


def test_real_result_consumer_refuses_completion_and_output_drift(
    coordinates, tmp_path
):
    adapter, inputs, workspace, runner, spec = compose_results(coordinates, tmp_path)
    result = runner.run(spec)
    assert result.is_success and result.value.exit_code == 0, result
    attempt = workspace / "hitrac-attempt"
    changes = [
        (attempt / "complete.json", lambda _: b"{}"),
        (
            workspace / "hitrac-request.json",
            lambda data: data.replace(str(workspace).encode(), b"/wrong-attempt"),
        ),
        (attempt / "private/inputs.json", lambda _: b"{}"),
        (attempt / "private/call-plan.json", lambda _: b"{}"),
        (
            attempt / "output/tracPre_summary.txt",
            lambda data: data.replace(b"s000002", b"s999999"),
        ),
        (attempt / "output/s000001/s000001_all.bedpe.gz", lambda data: data[:-8]),
        (attempt / "output/s000002/s000002_unique.bedpe.gz", lambda _: b"bad-gzip"),
        (attempt / "output/s000001/s000001.bam", lambda _: b"bad-bam"),
    ]
    summary_path = attempt / "output/tracPre_summary.txt"
    changes.extend(
        [
            (summary_path, lambda data: b"\n".join(data.splitlines()[:-1]) + b"\n"),
            (summary_path, lambda data: data + data.splitlines()[1] + b"\n"),
            (
                summary_path,
                lambda data: (
                    data + data.splitlines()[1].replace(b"s000001", b"s000003") + b"\n"
                ),
            ),
            (
                summary_path,
                lambda data: data.replace(b"mapping ratio", b"incorrect column"),
            ),
            (
                attempt / "output/s000001/s000001_all.bedpe.gz",
                lambda data: gzip.compress(
                    gzip.decompress(data).replace(b"\t", b" ", 1)
                ),
            ),
        ]
    )
    for path, mutate in changes:
        original = path.read_bytes()
        path.write_bytes(mutate(original))
        try:
            rejected = adapter.extract_artifacts(inputs, workspace)
            assert rejected.is_failure
            assert {issue.code for issue in rejected.issues} == {
                "HITRAC_RESULTS_INVALID"
            }
            assert not (workspace / "results").exists()
        finally:
            path.write_bytes(original)
    extra = attempt / "output/s999999"
    extra.mkdir()
    try:
        assert adapter.extract_artifacts(inputs, workspace).is_failure
    finally:
        extra.rmdir()
    partial = attempt / "output/s000002/s000002_unique.bedpe.gz"
    old_partial = partial.read_bytes()
    partial.unlink()
    try:
        assert adapter.extract_artifacts(inputs, workspace).is_failure
        assert not (workspace / "results").exists()
    finally:
        partial.write_bytes(old_partial)
    receipt = next((attempt / "private/calls").glob("*.end.json"))
    old_receipt = receipt.read_bytes()
    for damaged in (b"broken-json", b"{}"):
        receipt.write_bytes(damaged)
        try:
            assert adapter.extract_artifacts(inputs, workspace).is_failure
        finally:
            receipt.write_bytes(old_receipt)
    receipt.unlink()
    try:
        assert adapter.extract_artifacts(inputs, workspace).is_failure
    finally:
        receipt.write_bytes(old_receipt)
    marker = attempt / "complete.json"
    old = marker.read_bytes()
    marker.unlink()
    try:
        assert adapter.extract_artifacts(inputs, workspace).is_failure
    finally:
        marker.write_bytes(old)
    assert adapter.extract_artifacts(inputs, workspace).is_success


@pytest.mark.parametrize(
    "scenario,reason",
    [
        ("singlemate_only", "pair_support_missing"),
        ("multisample_one_empty", "empty_pet_set"),
    ],
)
def test_real_failed_qualification_cannot_be_published(
    coordinates, tmp_path, scenario, reason
):
    adapter, inputs, workspace, runner, spec = compose_results(
        coordinates, tmp_path, scenario=scenario
    )
    result = runner.run(spec)
    assert result.is_success and result.value.exit_code == 1, result
    assert json.loads(result.value.stdout)["reason_code"] == reason
    assert not (workspace / "hitrac-attempt/complete.json").exists()
    rejected = adapter.extract_artifacts(inputs, workspace)
    assert rejected.is_failure
    assert not (workspace / "results").exists()
