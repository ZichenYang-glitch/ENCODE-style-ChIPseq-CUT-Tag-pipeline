"""Real FastQC / MultiQC 1.35 contract checks for staged read QC."""

import os
from pathlib import Path
import shutil
import subprocess

import pytest

from _tool_resolver import require_external_tools, resolve_tool
from test_fastqc_multiqc_stages import (
    REPO_ROOT,
    TSV_HEADER,
    _fingerprint,
    _run_snakemake,
    _write_config,
    _write_fastq,
)


MULTIQC_CONFIG = REPO_ROOT / "workflow" / "multiqc_config.yaml"
MULTIQC = resolve_tool("multiqc", "MULTIQC")

pytestmark = pytest.mark.real_execution


def test_staged_fastqc_and_multiqc_135_report_contract(tmp_path):
    tool_bin = str(Path(MULTIQC).resolve().parent)
    env = os.environ.copy()
    env["PATH"] = tool_bin + os.pathsep + env.get("PATH", "")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    missing = [
        tool
        for tool in ("snakemake", "fastqc", "trim_galore", "multiqc")
        if not shutil.which(tool, path=env["PATH"])
    ]
    require_external_tools(missing, "staged FastQC / MultiQC real execution")

    samples = []
    raw_inputs = []
    for index, sample in enumerate(("sample.one", "sample.two"), start=1):
        input_dir = tmp_path / "inputs" / sample
        r1 = input_dir / "shared_R1.fastq.gz"
        r2 = input_dir / "shared_R2.fastq.gz"
        _write_fastq(r1, mate=1)
        _write_fastq(r2, mate=2)
        raw_inputs.extend((r1, r2))
        samples.append(
            f"{sample}\t{r1}\t{r2}\tPE\tchipseq\tH3K27ac\tnarrow\ths\t"
            f"unused\texp{index}\tcondition\t1\t1\t1\ttreatment\t\t"
        )

    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(
        TSV_HEADER + "\n" + "\n".join(samples) + "\n",
        encoding="utf-8",
    )
    config = _write_config(tmp_path, samples_path, trim=True)
    fingerprints_before = {path: _fingerprint(path) for path in raw_inputs}
    sentinels = [
        f"results/{sample}/logs/{sample}.fastqc.done"
        for sample in ("sample.one", "sample.two")
    ]

    result = _run_snakemake(tmp_path, config, sentinels, env=env)
    assert result.returncode == 0, result.stdout + result.stderr

    for sample in ("sample.one", "sample.two"):
        assert (tmp_path / f"results/{sample}/logs/{sample}.fastqc.done").is_file()
        for stage in ("raw", "trimmed"):
            for read in ("R1", "R2"):
                for extension in ("html", "zip"):
                    report = tmp_path / (
                        f"results/{sample}/01_qc/fastqc/{stage}/"
                        f"{sample}.{stage}.{read}_fastqc.{extension}"
                    )
                    assert report.is_file()
                    assert report.stat().st_size > 0
    assert {path: _fingerprint(path) for path in raw_inputs} == fingerprints_before

    # Concrete reports remain independently rebuildable after the sentinel exists.
    missing_report = tmp_path / (
        "results/sample.one/01_qc/fastqc/raw/sample.one.raw.R1_fastqc.html"
    )
    missing_report.unlink()
    rebuild = _run_snakemake(
        tmp_path,
        config,
        [str(missing_report.relative_to(tmp_path))],
        env=env,
    )
    assert rebuild.returncode == 0, rebuild.stdout + rebuild.stderr
    assert missing_report.stat().st_size > 0

    flagstat = (
        "100 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "0 + 0 secondary\n"
        "0 + 0 supplementary\n"
        "10 + 0 duplicates\n"
        "90 + 0 mapped (90.00% : N/A)\n"
        "100 + 0 paired in sequencing\n"
        "50 + 0 read1\n"
        "50 + 0 read2\n"
        "80 + 0 properly paired (80.00% : N/A)\n"
        "80 + 0 with itself and mate mapped\n"
        "10 + 0 singletons (10.00% : N/A)\n"
        "0 + 0 with mate mapped to a different chr\n"
        "0 + 0 with mate mapped to a different chr (mapQ>=5)\n"
    )
    for sample in ("sample.one", "sample.two"):
        qc_dir = tmp_path / f"results/{sample}/01_qc"
        (qc_dir / f"{sample}.flagstat.txt").write_text(flagstat, encoding="utf-8")
        (qc_dir / f"{sample}.final.flagstat.txt").write_text(
            flagstat.replace("100 + 0 in total", "80 + 0 in total").replace(
                "90 + 0 mapped", "75 + 0 mapped"
            ),
            encoding="utf-8",
        )
        (qc_dir / f"{sample}.idxstats.txt").write_text(
            "chr1\t1000\t90\t0\n*\t0\t0\t0\n",
            encoding="utf-8",
        )

    multiqc_dir = tmp_path / "multiqc"
    multiqc = subprocess.run(
        [
            MULTIQC,
            str(tmp_path / "results/sample.one"),
            str(tmp_path / "results/sample.two"),
            "--config",
            str(MULTIQC_CONFIG),
            "--outdir",
            str(multiqc_dir),
            "--filename",
            "multiqc_report.html",
            "--data-dir",
            "--force",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert multiqc.returncode == 0, multiqc.stdout + multiqc.stderr

    report_html = (multiqc_dir / "multiqc_report.html").read_text(encoding="utf-8")
    assert 'id="fastqc_raw"' in report_html
    assert 'id="fastqc_trimmed"' in report_html
    assert 'id="samtools_pre_filter"' in report_html
    assert 'id="samtools_final"' in report_html
    assert "fastqc_filtered" not in report_html
    assert "Seqs (raw)" in report_html
    assert "Seqs (trimmed)" in report_html

    data_dir = next(
        path for path in multiqc_dir.iterdir() if path.name.endswith("_data")
    )
    raw_data = (data_dir / "multiqc_fastqc_fastqc_raw.txt").read_text(encoding="utf-8")
    trimmed_data = (data_dir / "multiqc_fastqc_fastqc_trimmed.txt").read_text(
        encoding="utf-8"
    )
    for sample in ("sample.one", "sample.two"):
        for read in ("R1", "R2"):
            assert f"{sample}.raw.{read}" in raw_data
            assert f"{sample}.trimmed.{read}" in trimmed_data

    pre_filter = (
        data_dir / "multiqc_samtools_flagstat_samtools_pre_filter.txt"
    ).read_text(encoding="utf-8")
    final = (data_dir / "multiqc_samtools_flagstat_samtools_final.txt").read_text(
        encoding="utf-8"
    )
    pre_idxstats = (
        data_dir / "multiqc_samtools_idxstats_samtools_pre_filter.txt"
    ).read_text(encoding="utf-8")
    for sample in ("sample.one", "sample.two"):
        assert sample in pre_filter
        assert sample in final
        assert sample in pre_idxstats
    assert not (data_dir / "multiqc_samtools_idxstats_samtools_final.txt").exists()
