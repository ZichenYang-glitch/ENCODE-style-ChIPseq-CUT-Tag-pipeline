"""Original-rule input and argv contracts, without requiring science tools."""

import re
import shlex

import pytest

from _complexity import project, snakemake


def job(output, rule):
    match = re.search(rf"(?m)^(?:local)?rule {rule}:\n((?:[ \t]+[^\n]*\n)+)", output)
    assert match is not None, output
    return match.group(1)


@pytest.mark.parametrize("layout", ["SE", "PE"])
@pytest.mark.parametrize(
    "policy,mode,effective,mapq",
    [
        ("yes", "narrow", "yes", 17),
        ("no", "narrow", "no", 42),
        ("auto", "narrow", "yes", 17),
        ("auto", "broad", "no", 42),
    ],
)
def test_complexity_uses_filtered_duplicates_independent_of_dedup_policy(
    tmp_path, snakemake_executable, run_validator, layout, policy, mode, effective, mapq
):
    config = project(tmp_path, layout, remove_dup=policy, peak_mode=mode, mapq=mapq)
    assert run_validator(config).returncode == 0
    # A placeholder supplies only DAG existence, never tool/scientific evidence.
    sorted_bam = tmp_path / "results/S1/02_align/S1.sorted.bam"
    sorted_bam.parent.mkdir(parents=True)
    sorted_bam.write_bytes(b"DAG-only placeholder")
    result = snakemake(
        tmp_path,
        snakemake_executable,
        config,
        [
            "results/S1/01_qc/S1.nrf_pbc.tsv",
            "results/S1/05_qc/preseq/S1.preseq.txt",
            "results/S1/02_align/S1.final.bam",
        ],
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    filtered = f"results/S1/02_align/S1.mapq{mapq}.bam"
    for rule in ("nrf_pbc", "preseq_complexity"):
        inputs = re.search(r"(?m)^    input: (.+)$", job(output, rule)).group(1)
        assert inputs == filtered
    command = (
        re.search(r"(?m)^\s*preseq lc_extrap [^\n]+", output).group().rstrip(" \\")
    )
    argv = shlex.split(command)
    assert argv == [
        "preseq",
        "lc_extrap",
        "-B",
        *(["-P"] if layout == "PE" else []),
        "-o",
        "results/S1/05_qc/preseq/S1.preseq.txt",
        filtered,
        "2>&1",
        "|",
        "tee",
        "results/S1/logs/S1.preseq.log",
    ]
    assert f'if [[ "{effective}" == "yes" ]]' in output
    assert f"-q {mapq} -F 0x904" in output


def test_preseq_remains_opt_in_and_nrf_remains_default(tmp_path, snakemake_executable):
    config = project(tmp_path, "PE", preseq=None)
    result = snakemake(tmp_path, snakemake_executable, config, restricted=False)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert job(output, "nrf_pbc")
    assert re.search(r"(?m)^(?:local)?rule preseq_complexity:", output) is None
