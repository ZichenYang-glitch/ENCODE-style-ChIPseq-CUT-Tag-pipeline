"""Dry-run contracts for CUT&Tag fragment-size and SEACR sidecars."""

from __future__ import annotations

import re
from pathlib import Path

import pytest


pytestmark = pytest.mark.full_main


SAMPLE_COLUMNS = (
    "sample",
    "fastq_1",
    "fastq_2",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "experiment",
    "condition",
    "replicate",
    "biological_replicate",
    "technical_replicate",
    "role",
    "control_sample",
    "control_bam",
)

SEACR_RULES = {"seacr_bedgraph", "seacr_call"}


def _sample_sheet(tmp_path: Path, sample_specs: list[dict[str, str]]) -> str:
    rows = ["\t".join(SAMPLE_COLUMNS)]
    for index, spec in enumerate(sample_specs, start=1):
        sample = spec["sample"]
        layout = spec.get("layout", "PE")
        fastq_1 = tmp_path / f"{sample}_R1.fastq"
        fastq_1.write_text("", encoding="utf-8")
        fastq_2 = tmp_path / f"{sample}_R2.fastq"
        if layout == "PE":
            fastq_2.write_text("", encoding="utf-8")

        assay = spec.get("assay", "cuttag")
        target = spec.get("target", "H3K4me3")
        row = {
            "sample": sample,
            "fastq_1": str(fastq_1),
            "fastq_2": str(fastq_2) if layout == "PE" else "",
            "layout": layout,
            "assay": assay,
            "target": target,
            "peak_mode": "narrow",
            "genome": "hs",
            "bowtie2_index": "idx",
            "experiment": spec.get("experiment", f"exp_{sample}"),
            "condition": target,
            "replicate": spec.get("replicate", "1"),
            "biological_replicate": spec.get("biological_replicate", "1"),
            "technical_replicate": str(index),
            "role": spec.get("role", "treatment"),
            "control_sample": spec.get("control_sample", ""),
            "control_bam": "",
        }
        rows.append("\t".join(row[column] for column in SAMPLE_COLUMNS))
    return "\n".join(rows) + "\n"


def _base_config(tmp_path: Path, *, use_control: bool = False) -> dict:
    return {
        "outdir": str(tmp_path / "results"),
        "use_control": use_control,
        "threads": 1,
        "trim": False,
        "multiqc": False,
        "stage4b": False,
        "stage5": False,
        "qc": {
            "blacklist_filter": False,
            "frip": False,
            "library_complexity": False,
            "nrf_pbc": False,
            "signal_tracks": False,
            "summary": False,
            "cross_correlation": False,
            "preseq_complexity": False,
            "picard_metrics": False,
            "tss_enrichment": False,
        },
        "genome_resources": {
            "hs": {
                "effective_genome_size": "hs",
                "chrom_sizes": "",
                "blacklist": "",
                "gtf": "",
                "reference_fasta": "",
            }
        },
    }


def _scheduled_rules(output: str) -> set[str]:
    return {
        match.group(1)
        for line in output.splitlines()
        if (match := re.match(r"^\s*(?:local)?rule\s+(\w+)\s*:", line))
    }


def _dry_run(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
    *,
    sample_specs,
    config,
):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / ".cache"))
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    _, config_path, _ = tmp_config(
        config=config,
        samples=_sample_sheet(tmp_path, sample_specs),
    )
    result = run_snakemake(config_path, extra_args=["--printshellcmds"])
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-5000:]
    return output, _scheduled_rules(output)


@pytest.mark.parametrize("layout", ["PE", "SE"])
def test_fragment_size_schedules_only_cuttag_with_the_sample_layout(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
    layout,
):
    config = _base_config(tmp_path)
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        monkeypatch,
        sample_specs=[
            {"sample": "CUTTAG", "assay": "cuttag", "layout": layout},
            {"sample": "CHIP", "assay": "chipseq", "layout": "PE"},
        ],
        config=config,
    )

    cuttag_output = tmp_path / "results/CUTTAG/01_qc/CUTTAG.cuttag_fragment_size.tsv"
    chip_output = tmp_path / "results/CHIP/01_qc/CHIP.cuttag_fragment_size.tsv"
    assert "cuttag_fragment_size" in rules
    assert str(cuttag_output) in output
    assert str(chip_output) not in output
    assert "python3 scripts/calc_cuttag_fragment_size.py" in output
    assert "--sample CUTTAG" in output
    assert f"--layout {layout}" in output


def test_fragment_size_false_omits_the_sidecar(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
):
    config = _base_config(tmp_path)
    config["qc"]["cuttag_fragment_size"] = False
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        monkeypatch,
        sample_specs=[{"sample": "CUTTAG", "assay": "cuttag"}],
        config=config,
    )

    assert "cuttag_fragment_size" not in rules
    assert ".cuttag_fragment_size.tsv" not in output
    assert "calc_cuttag_fragment_size.py" not in output


def test_fragment_size_includes_cuttag_control_but_not_other_assays(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
):
    config = _base_config(tmp_path, use_control=True)
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        monkeypatch,
        sample_specs=[
            {
                "sample": "CUTTAG_T",
                "assay": "cuttag",
                "experiment": "exp_cuttag",
                "control_sample": "CUTTAG_C",
            },
            {
                "sample": "CUTTAG_C",
                "assay": "cuttag",
                "layout": "SE",
                "target": "IgG",
                "experiment": "exp_cuttag",
                "role": "control",
            },
            {"sample": "CHIP", "assay": "chipseq", "layout": "SE"},
        ],
        config=config,
    )

    assert "cuttag_fragment_size" in rules
    for sample in ("CUTTAG_T", "CUTTAG_C"):
        expected = (
            tmp_path / f"results/{sample}/01_qc/{sample}.cuttag_fragment_size.tsv"
        )
        assert str(expected) in output
    excluded = tmp_path / "results/CHIP/01_qc/CHIP.cuttag_fragment_size.tsv"
    assert str(excluded) not in output


@pytest.mark.parametrize("explicitly_disabled", [False, True])
def test_seacr_default_and_false_do_not_schedule_sidecars(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
    explicitly_disabled,
):
    config = _base_config(tmp_path)
    config["qc"]["cuttag_fragment_size"] = False
    if explicitly_disabled:
        config["cuttag"] = {"seacr": {"enabled": False}}
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        monkeypatch,
        sample_specs=[{"sample": "CUTTAG", "assay": "cuttag"}],
        config=config,
    )

    assert rules.isdisjoint(SEACR_RULES)
    assert "/04_peaks_seacr/" not in output


def test_seacr_enabled_cuttag_pe_wires_bedgraph_and_caller_command(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
):
    config = _base_config(tmp_path)
    config["qc"]["cuttag_fragment_size"] = False
    config["cuttag"] = {
        "seacr": {
            "enabled": True,
            "mode": "relaxed",
            "normalization": "non",
            "threshold": 0.05,
        }
    }
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        monkeypatch,
        sample_specs=[{"sample": "CUTTAG", "assay": "cuttag", "layout": "PE"}],
        config=config,
    )
    shell = " ".join(output.replace("\\", " ").split())

    bedgraph = tmp_path / "results/CUTTAG/04_peaks_seacr/CUTTAG.bedgraph"
    final_bam = tmp_path / "results/CUTTAG/02_align/CUTTAG.final.bam"
    peak = tmp_path / "results/CUTTAG/04_peaks_seacr/CUTTAG/CUTTAG.seacr.relaxed.bed"
    assert SEACR_RULES <= rules
    assert f"bedtools genomecov -ibam {final_bam} -bg -pc > {bedgraph}" in shell
    assert "SEACR_SCRIPT=$(command -v SEACR_1.3.sh)" in shell
    assert f'bash "$SEACR_SCRIPT" {bedgraph} 0.05 non relaxed' in shell
    assert str(peak) in shell
    assert f"$(dirname {peak})/CUTTAG.seacr" in shell


@pytest.mark.parametrize(
    ("assay", "layout"),
    [("chipseq", "PE"), ("cuttag", "SE")],
    ids=["chipseq-excluded", "single-end-cuttag-excluded"],
)
def test_seacr_enabled_excludes_unsupported_samples(
    tmp_path,
    tmp_config,
    run_snakemake,
    monkeypatch,
    assay,
    layout,
):
    config = _base_config(tmp_path)
    config["qc"]["cuttag_fragment_size"] = False
    config["cuttag"] = {"seacr": {"enabled": True}}
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        monkeypatch,
        sample_specs=[{"sample": "SAMPLE", "assay": assay, "layout": layout}],
        config=config,
    )

    assert rules.isdisjoint(SEACR_RULES)
    assert "/04_peaks_seacr/" not in output
