"""Snakemake command contracts for validated tool parameters."""

from __future__ import annotations


SAMPLE_COLUMNS = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\n"
)


def _dry_run(tmp_path, tmp_config, run_snakemake, tool_parameters=None):
    read = tmp_path / "reads.fastq"
    read.write_text("", encoding="utf-8")
    samples = SAMPLE_COLUMNS + (f"S1\t{read}\t\tSE\tchipseq\tCTCF\tnarrow\ths\tidx\n")
    config = {
        "outdir": str(tmp_path / "results"),
        "use_control": False,
        "threads": 1,
        "stage4b": False,
        "stage5": False,
    }
    if tool_parameters is not None:
        config["tool_parameters"] = tool_parameters
    _, config_path, _ = tmp_config(config=config, samples=samples)

    result = run_snakemake(
        config_path,
        extra_args=["--printshellcmds"],
        quiet=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr


def test_default_policy_arguments_reach_shell_commands(
    tmp_path, tmp_config, run_snakemake
):
    output = _dry_run(tmp_path, tmp_config, run_snakemake)

    assert "-q 0.01" in output
    assert "-F 0x904" in output
    assert "--quality" not in output
    assert "--length" not in output
    assert "--stringency" not in output
    assert "--very-sensitive" not in output


def test_explicit_tool_parameters_reach_their_shell_commands(
    tmp_path, tmp_config, run_snakemake
):
    output = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        tool_parameters={
            "fastqc": {"extra_args": "--nogroup"},
            "trim_galore": {"quality": 25},
            "macs3": {"extra_args": "--keep-dup all"},
        },
    )

    assert "--nogroup" in output
    assert "--quality 25" in output
    assert "--keep-dup all" in output
