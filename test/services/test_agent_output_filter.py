import pytest

from encode_pipeline.services.agent_output_filter import OutputFilter, _SAFE_MESSAGE


def test_filters_snakemake_command():
    filter_ = OutputFilter()
    result = filter_.filter_text("Run: snakemake --cores 4 workflow.smk")
    assert result.filtered is True
    assert "filtered because it contained execution-like instructions" in result.text
    assert "snakemake_command" in result.matched_patterns


@pytest.mark.parametrize("unsafe", [
    "$ snakemake --cores 4",
    "rm -rf /tmp/data",
    "rm --force /important/data",
    "python script.py --input x",
    "python3 script.py --input x",
    "Run: python3 script.py --input x",
    "    python3 script.py",
    "bash run.sh",
    "Example: bash run.sh",
    "sbatch job.sh",
    "qsub job.sh",
    "Please run this command.",
    "execute this shell command now",
    "delete this file",
])
def test_filters_unsafe_commands(unsafe: str):
    filter_ = OutputFilter()
    result = filter_.filter_text(unsafe)
    assert result.filtered is True
    assert result.text == _SAFE_MESSAGE
    assert len(result.matched_patterns) > 0


def test_allows_snakemake_prose():
    filter_ = OutputFilter()
    text = "This workflow is implemented with Snakemake."
    result = filter_.filter_text(text)
    assert result.filtered is False
    assert result.text == text


def test_allows_normal_explanation():
    filter_ = OutputFilter()
    text = "FRiP measures the fraction of reads in peaks."
    result = filter_.filter_text(text)
    assert result.filtered is False
    assert result.text == text
