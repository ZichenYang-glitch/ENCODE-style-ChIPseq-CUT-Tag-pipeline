"""Tests for the ENCODE-style workflow adapter."""

import json
import os
from pathlib import Path
import subprocess
import sys

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import WorkflowInputs


def test_workspace_config_yaml_is_canonical_across_mapping_order():
    from encode_pipeline.adapters.encode import _render_config_yaml

    first = {
        "tool_parameters": {
            "samtools_filter": {"filter_flags": "", "extra_args": ""},
            "fastqc": {"extra_args": ""},
        },
        "mnase": {"callers": {"danpos3": False, "sem": False}},
    }
    second = {
        "mnase": {"callers": {"sem": False, "danpos3": False}},
        "tool_parameters": {
            "fastqc": {"extra_args": ""},
            "samtools_filter": {"extra_args": "", "filter_flags": ""},
        },
    }

    assert _render_config_yaml(first) == _render_config_yaml(second)


def test_encode_adapter_declares_workspace_plan_capability():
    adapter = EncodeStyleWorkflowAdapter()
    assert "workspace_plan" in adapter.capabilities.supports


def test_adapter_plan_workspace_returns_config_file(tmp_path):
    import yaml

    config = {
        "samples": str(tmp_path / "samples.tsv"),
        "outdir": "custom",
        "threads": 1,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }
    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )

    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config=config, samples=str(samples_tsv), options={})
    result = adapter.plan_workspace(inputs, str(tmp_path / "workspace"))

    assert result.is_success is True
    plan = result.value
    assert plan.directories == ("logs", "results")
    files = dict(plan.files)
    assert "config/config.yaml" in files

    parsed = yaml.safe_load(files["config/config.yaml"])
    assert parsed["samples"] == "config/samples.tsv"
    assert parsed["outdir"] == "results"


def test_adapter_plan_workspace_translates_semantic_config_only_at_materialization(
    tmp_path,
):
    import yaml

    from encode_pipeline.config.validate import validate_config

    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )
    submitted_config = {
        "replicate_analysis": {"enabled": False},
        "chipseq_idr": {"enabled": False},
        "threads": 1,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }

    result = EncodeStyleWorkflowAdapter().plan_workspace(
        WorkflowInputs(
            config=submitted_config,
            samples=str(samples_tsv),
            options={},
        ),
        str(tmp_path / "workspace"),
    )

    assert result.is_success
    rendered = yaml.safe_load(dict(result.value.files)["config/config.yaml"])
    assert rendered["replicate_analysis"] is False
    assert rendered["chipseq_idr"] is False
    assert "stage4b" not in rendered
    assert "stage5" not in rendered
    assert validate_config(rendered)["replicate_analysis"] is False
    assert validate_config(rendered)["chipseq_idr"] is False
    assert submitted_config == {
        "replicate_analysis": {"enabled": False},
        "chipseq_idr": {"enabled": False},
        "threads": 1,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }


def test_adapter_plan_workspace_semantic_config_has_only_completion_issue(tmp_path):
    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )

    result = EncodeStyleWorkflowAdapter().plan_workspace(
        WorkflowInputs(
            config={
                "replicate_analysis": {"enabled": False},
                "chipseq_idr": {"enabled": False},
                "blacklist_mode": "off",
            },
            samples=str(samples_tsv),
        ),
        str(tmp_path / "workspace"),
    )

    assert result.is_success
    assert [issue.code for issue in result.issues] == [
        "ENCODE_WORKSPACE_PLANNING_COMPLETE",
    ]


def test_adapter_plan_workspace_round_trip_pe_control_replicate(tmp_path):
    import csv
    import io

    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\tcontrol_bam\trole\tcontrol_sample\texperiment\tcondition\t"
        "replicate\tbiological_replicate\ttechnical_replicate\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t"
        "/abs/bt2/GRCh38\t\ttreatment\t\texp1\tcond1\t1\t1\t1\n"
        "S2\t/abs/S2_1.fq.gz\t/abs/S2_2.fq.gz\tPE\tchipseq\tINPUT\tnarrow\ths\t"
        "/abs/bt2/GRCh38\t\tcontrol\t\texp1\tcond1\t1\t1\t2\n",
        encoding="utf-8",
    )

    config = {
        "samples": str(samples_tsv),
        "threads": 1,
        "use_control": False,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }

    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config=config, samples=str(samples_tsv), options={})
    result = adapter.plan_workspace(inputs, str(tmp_path / "workspace"))

    assert result.is_success is True
    plan = result.value
    tsv_bytes = dict(plan.files)["config/samples.tsv"]
    reader = csv.DictReader(io.StringIO(tsv_bytes.decode("utf-8")), delimiter="\t")
    rows = list(reader)
    assert len(rows) == 2
    assert rows[0]["sample"] == "S1"
    assert rows[0]["fastq_2"] == "/abs/S1_2.fq.gz"
    assert rows[0]["biological_replicate"] == "1"
    assert rows[1]["role"] == "control"
    assert rows[1]["technical_replicate"] == "2"


def test_adapter_plan_workspace_rejects_relative_fastq_path(tmp_path):
    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\trelative/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )

    config = {
        "samples": str(samples_tsv),
        "threads": 1,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }

    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config=config, samples=str(samples_tsv), options={})
    result = adapter.plan_workspace(inputs, str(tmp_path / "workspace"))

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "ENCODE_WORKSPACE_RELATIVE_EXTERNAL_PATH"
    assert issue.path == "samples[0].fastq_1"
    assert issue.context == {"field": "fastq_1"}
    payload = issue.to_dict()
    assert "relative/S1_1.fq.gz" not in str(payload)


def test_adapter_plan_workspace_rejects_relative_genome_resource_path(
    tmp_path, monkeypatch
):
    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )

    relative_ref = tmp_path / "relative_ref.fa"
    relative_ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    config = {
        "samples": str(samples_tsv),
        "threads": 1,
        "blacklist_mode": "off",
        "genome_resources": {
            "hs": {
                "effective_genome_size": "hs",
                "reference_fasta": "relative_ref.fa",
            }
        },
    }

    monkeypatch.chdir(tmp_path)

    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config=config, samples=str(samples_tsv), options={})
    result = adapter.plan_workspace(inputs, str(tmp_path / "workspace"))

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "ENCODE_WORKSPACE_RELATIVE_EXTERNAL_PATH"
    assert issue.path == "config.genome_resources.hs.reference_fasta"


def test_adapter_plan_workspace_issue_does_not_leak_user_path(tmp_path):
    missing_path = str(tmp_path / "does_not_exist.tsv")
    config = {
        "samples": missing_path,
        "threads": 1,
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }

    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config=config, samples=missing_path, options={})
    result = adapter.plan_workspace(inputs, str(tmp_path / "workspace"))

    assert result.is_failure is True
    issue = result.issues[0]
    payload = issue.to_dict()
    for value in payload.values():
        if isinstance(value, str):
            assert missing_path not in value
    for ctx_value in payload["context"].values():
        if isinstance(ctx_value, str):
            assert missing_path not in ctx_value
    assert issue.hint is None
    assert issue.context == {}


def test_adapter_plan_workspace_accepts_inline_rows_and_is_byte_deterministic(
    tmp_path,
):
    adapter = EncodeStyleWorkflowAdapter()
    row = {
        "sample": "S1",
        "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
        "fastq_2": str((tmp_path / "S1.R2.fastq.gz").resolve()),
        "layout": "PE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((tmp_path / "indices/hs").resolve()),
    }
    inputs = WorkflowInputs(
        config={"blacklist_mode": "off"},
        samples=[row],
        options={},
    )
    first = adapter.plan_workspace(inputs, tmp_path / "workspace-first")
    second = adapter.plan_workspace(inputs, tmp_path / "workspace-second")

    assert first.is_success is True
    assert second.is_success is True
    assert first.value.files == second.value.files
    files = dict(first.value.files)
    assert b"samples: config/samples.tsv" in files["config/config.yaml"]
    assert b"S1" in files["config/samples.tsv"]


_PLAN_HASH_CHILD = """
import hashlib
import json
import sys
from pathlib import Path

base = Path(sys.argv[1])
sys.path.insert(0, str(Path(sys.argv[2]) / "src"))

import encode_pipeline.adapters.encode as encode_module
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import WorkflowInputs

row = {
    "sample": "S1",
    "fastq_1": str(base / "S1.R1.fastq.gz"),
    "fastq_2": str(base / "S1.R2.fastq.gz"),
    "layout": "PE",
    "assay": "chipseq",
    "target": "CTCF",
    "peak_mode": "narrow",
    "genome": "hs",
    "bowtie2_index": str(base / "indices" / "hs"),
}
adapter = EncodeStyleWorkflowAdapter()
result = adapter.plan_workspace(
    WorkflowInputs(config={"blacklist_mode": "off"}, samples=[row]), base / "workspace"
)
if result.is_failure:
    print(json.dumps({"codes": [issue.code for issue in result.issues]}))
    raise SystemExit(1)
hashes = {
    relpath: hashlib.sha256(content).hexdigest()
    for relpath, content in sorted(result.value.files)
}
hashes["__module_file__"] = encode_module.__file__
probe_keys = frozenset(f"probe-key-{index}" for index in range(16))
hashes["__hash_order_fingerprint__"] = hashlib.sha256(
    ",".join(probe_keys).encode("utf-8")
).hexdigest()
print(json.dumps(hashes, sort_keys=True))
"""


def test_adapter_plan_workspace_is_byte_deterministic_across_hash_seeds(tmp_path):
    # No -I here: -I implies -E, which silently ignores PYTHONHASHSEED and
    # would turn this into a comparison of two random seeds.
    repo_root = Path(__file__).resolve().parents[2]
    base_env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(("PYTEST", "PYTHONPATH", "PYTHONHOME"))
    }
    expected_module = str(
        repo_root / "src" / "encode_pipeline" / "adapters" / "encode.py"
    )
    runs = {}
    for seed in ("1", "2"):
        env = dict(base_env)
        env["PYTHONHASHSEED"] = seed
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                "-c",
                _PLAN_HASH_CHILD,
                str(tmp_path),
                str(repo_root),
            ],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr[-500:]
        runs[seed] = json.loads(completed.stdout)

    fingerprints = {}
    for seed, hashes in runs.items():
        module_file = hashes.pop("__module_file__")
        assert module_file == expected_module
        fingerprints[seed] = hashes.pop("__hash_order_fingerprint__")
    # Self-proof: the two children really ran with different hash seeds.
    # Without this, a broken env handoff would silently compare identical seeds.
    assert fingerprints["1"] != fingerprints["2"]
    assert runs["1"] == runs["2"]


def test_inline_and_server_path_modes_render_equivalent_workspace_files(tmp_path):
    row = {
        "sample": "S1",
        "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
        "fastq_2": str((tmp_path / "S1.R2.fastq.gz").resolve()),
        "layout": "PE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((tmp_path / "indices/hs").resolve()),
    }
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(
        "\t".join(row) + "\n" + "\t".join(row.values()) + "\n",
        encoding="utf-8",
    )
    adapter = EncodeStyleWorkflowAdapter()

    inline = adapter.plan_workspace(
        WorkflowInputs(config={"blacklist_mode": "off"}, samples=[row]),
        tmp_path / "workspace-inline",
    )
    server_path = adapter.plan_workspace(
        WorkflowInputs(config={"samples": str(samples_path), "blacklist_mode": "off"}),
        tmp_path / "workspace-path",
    )

    assert inline.is_success is True
    assert server_path.is_success is True
    assert inline.value.files == server_path.value.files


def test_adapter_plan_workspace_semantic_round_trip_via_loader(tmp_path):
    import csv
    import io

    from encode_pipeline.samples.load import load_and_validate_samples

    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\tcontrol_bam\trole\tcontrol_sample\texperiment\tcondition\t"
        "replicate\tbiological_replicate\ttechnical_replicate\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t"
        "/abs/bt2/GRCh38\t\ttreatment\tS2\texp1\tcond1\t1\t1\t1\n"
        "S2\t/abs/S2_1.fq.gz\t/abs/S2_2.fq.gz\tPE\tchipseq\tINPUT\tnarrow\ths\t"
        "/abs/bt2/GRCh38\t\tcontrol\t\texp1\tcond1\t1\t1\t2\n",
        encoding="utf-8",
    )

    config = {
        "samples": str(samples_tsv),
        "threads": 1,
        "use_control": True,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }

    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config=config, samples=str(samples_tsv), options={})
    result = adapter.plan_workspace(inputs, str(tmp_path / "workspace"))

    assert result.is_success is True
    plan = result.value
    tsv_bytes = dict(plan.files)["config/samples.tsv"]

    reader = csv.DictReader(io.StringIO(tsv_bytes.decode("utf-8")), delimiter="\t")
    assert reader.fieldnames == [
        "sample",
        "fastq_1",
        "fastq_2",
        "layout",
        "assay",
        "target",
        "peak_mode",
        "genome",
        "bowtie2_index",
        "control_bam",
        "role",
        "control_sample",
        "experiment",
        "condition",
        "replicate",
        "biological_replicate",
        "technical_replicate",
    ]

    rendered = tmp_path / "rendered_samples.tsv"
    rendered.write_bytes(tsv_bytes)

    loaded = load_and_validate_samples(
        str(rendered),
        use_control=True,
        chipseq_idr_enabled=False,
        strict_inputs=False,
    )

    assert len(loaded) == 2
    s1 = next(r for r in loaded if r["id"] == "S1")
    s2 = next(r for r in loaded if r["id"] == "S2")
    assert s1["layout"] == "PE"
    assert s1["control_sample"] == "S2"
    assert s2["role"] == "control"
    assert s1["biological_replicate"] == 1
    assert s1["technical_replicate"] == 1
    assert s2["technical_replicate"] == 2
