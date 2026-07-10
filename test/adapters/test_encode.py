"""Tests for the ENCODE-style workflow adapter."""

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import WorkflowInputs


def test_encode_adapter_declares_workspace_plan_capability():
    adapter = EncodeStyleWorkflowAdapter()
    assert "workspace_plan" in adapter.capabilities.supports


def test_adapter_plan_workspace_returns_config_file(tmp_path):
    import yaml

    config = {
        "samples": str(tmp_path / "samples.tsv"),
        "outdir": "custom",
        "threads": 1,
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


def test_adapter_plan_workspace_rejects_relative_genome_resource_path(tmp_path, monkeypatch):
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
