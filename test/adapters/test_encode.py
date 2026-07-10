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
