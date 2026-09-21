"""Served schema -> real workbench defaults -> adapter -> tiny execution.

Run in the existing browser CI tier, where locked Node and Python dependencies
are both installed. These checks deliberately do not substitute YAML fixtures
for the bulk form defaults or mock the ENCODE command builder / Snakemake.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from api_test_client import ApiTestClient
from encode_pipeline.cli.results_visibility_fixture import (
    prepare_results_visibility_fixture,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.services.defaults import create_default_command_builder
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import WorkspacePlanner
from encode_pipeline.services.process_runner import ProcessRunner


pytestmark = pytest.mark.platform_real_execution
ROOT = Path(__file__).resolve().parents[2]
ENCODE_ID = "encode-style-chipseq-cuttag-atac-mnase"


@pytest.fixture
def client(reference_ready_app):
    try:
        with ApiTestClient(reference_ready_app) as client:
            yield client
    finally:
        reference_ready_app.state.persistence.close()


def _materialized_defaults(client, workflow_id):
    response = client.get(f"/api/v1/workflows/{workflow_id}/schema")
    assert response.status_code == 200
    assert response.json()["ok"] is True
    completed = subprocess.run(
        [
            "node",
            "node_modules/vite-node/vite-node.mjs",
            "e2e/materialize-authoring-schema.ts",
        ],
        cwd=ROOT / "frontend",
        input=json.dumps(response.json()["schema"]),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout)


def test_served_bulk_defaults_pass_real_adapter_validation(client, tmp_path):
    defaults = _materialized_defaults(client, "bulk-rnaseq")
    config = defaults["config"]
    assert "reference" not in config["standard"]
    # The platform supplies this operator-owned reference after authoring. Only
    # that required binding is added; every gated-section default is untouched.
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text('chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")
    config["standard"]["reference"] = {
        "reference_id": "smoke-reference",
        "fasta": str(fasta),
        "fasta_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        "gtf": str(gtf),
        "gtf_sha256": hashlib.sha256(gtf.read_bytes()).hexdigest(),
        "annotation_style": "ensembl",
    }
    reads = tmp_path / "S1.fastq.gz"
    reads.write_bytes(gzip.compress(b"@read1\nACGT\n+\nIIII\n", mtime=0))
    inputs = WorkflowInputs(
        config=config,
        samples=[
            {
                "sample": "S1",
                "library": "lib1",
                "lane": "L001",
                "layout": "SE",
                "fastq_1": str(reads),
                "strandedness": "auto",
                "platform": "ILLUMINA",
            }
        ],
        options=defaults["options"],
    )
    result = client.app.state.registry.get("bulk-rnaseq").validate(inputs)
    assert result.is_success, result.issues
    for section in ("umi", "ribosomal_rna_removal"):
        assert config["standard"][section] == {"enabled": False}


def test_encode_cores_reach_command_and_real_snakemake_allocation(
    client, tmp_path, monkeypatch
):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    options = _materialized_defaults(client, ENCODE_ID)["options"]
    assert options["cores"] == 1
    options["cores"] = 2
    project = tmp_path / "runtime" / "project"
    inputs = prepare_results_visibility_fixture(project, repository_root=ROOT)
    snakefile = project / "workflow" / "Snakefile"
    workflow = snakefile.read_text(encoding="utf-8")
    original_target = "rule all:\n    input:\n        RESULT_OUTPUTS\n"
    assert workflow.count(original_target) == 1
    # Extend only this task-owned fixture. The existing tiny sample/results
    # task still executes; this extra job records Snakemake's allocated threads.
    snakefile.write_text(
        workflow.replace(
            original_target,
            original_target.rstrip() + ' + ["result/allocated_threads.txt"]\n',
        )
        + "\nrule allocated_threads:\n"
        '    output: "result/allocated_threads.txt"\n'
        "    threads: 8\n"
        "    shell: \"printf '%s\\\\n' {threads} > {output:q}\"\n",
        encoding="utf-8",
    )
    with inputs.samples_path.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    registry = client.app.state.registry
    builder = create_default_command_builder(registry, project_root=project)
    pending = ExecutionPlan(
        plan_id="submission-smoke-plan",
        run_id="submission-smoke-run",
        workflow_id=ENCODE_ID,
        status=PlanStatus.PENDING,
        inputs_snapshot={
            "config": inputs.results_config,
            "samples": samples,
            "options": options,
        },
    )
    workspace = tmp_path / "workspace"
    planned = WorkspacePlanner(registry).plan_workspace(pending, workspace)
    assert planned.is_success, planned.issues
    materialized = WorkspaceMaterializer().materialize(
        planned.value.workspace_plan, workspace
    )
    assert materialized.is_success, materialized.issues
    built = builder.build_command(planned.value, workspace)
    assert built.is_success, built.issues
    command = built.value.command_spec
    executed = ProcessRunner(timeout_seconds=60).run(command)
    assert executed.is_success, executed.issues
    assert executed.value.exit_code == 0, executed.value
    assert (workspace / "result/complete.txt").read_text() == "success\n"
    assert (workspace / "results/C1/01_qc/C1.qc_summary.tsv").read_text() == (
        inputs.expected_qc_summary
    )
    allocated = int((workspace / "result/allocated_threads.txt").read_text())
    assert (command.argv[command.argv.index("--cores") + 1], allocated) == ("2", 2)
