"""Dry-run contracts for the consolidated reproducibility IDR rules."""

from __future__ import annotations

from pathlib import Path
import re

import pytest


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
    "biological_replicate",
)

NARROW_RULES = {
    "idr_macs3_biorep_narrow",
    "idr_true_replicates_narrow",
    "idr_split_pseudoreps_narrow",
    "idr_macs3_pseudorep_narrow",
    "idr_self_pseudoreps_narrow",
    "idr_pooled_pseudoreps_narrow",
}
BROAD_RULES = {
    "idr_macs3_biorep_broad",
    "idr_true_replicates_broad",
    "idr_split_pseudoreps_broad",
    "idr_macs3_pseudorep_broad",
    "idr_self_pseudoreps_broad",
    "idr_pooled_pseudoreps_broad",
}

MODE_CASES = [
    pytest.param(
        "atac",
        "narrow",
        "atac_narrow",
        "idr_summary_atac_narrow",
        NARROW_RULES,
        "atac.macs3.narrow.replicate_validated.idr.narrowPeak",
        id="atac-narrow",
    ),
    pytest.param(
        "cuttag",
        "narrow",
        "cuttag_narrow",
        "idr_summary_cuttag_narrow",
        NARROW_RULES,
        "cuttag.macs3.narrow.replicate_validated.idr.narrowPeak",
        id="cuttag-narrow",
    ),
    pytest.param(
        "chipseq",
        "broad",
        "chipseq_broad_experimental",
        "idr_summary_chipseq_broad",
        BROAD_RULES,
        "chipseq.macs3.broad.replicate_validated.idr.broadPeak",
        id="chipseq-broad",
    ),
    pytest.param(
        "cuttag",
        "broad",
        "cuttag_broad_experimental",
        "idr_summary_cuttag_broad",
        BROAD_RULES,
        "cuttag.macs3.broad.replicate_validated.idr.broadPeak",
        id="cuttag-broad",
    ),
]


def _sample_sheet(
    tmp_path: Path,
    rows: list[tuple[str, str, str, str, str, str]],
) -> str:
    rendered = ["\t".join(SAMPLE_COLUMNS)]
    for sample, assay, peak_mode, experiment, biorep, layout in rows:
        fastq_1 = tmp_path / f"{sample}_R1.fq"
        fastq_1.write_text("", encoding="utf-8")
        fastq_2 = tmp_path / f"{sample}_R2.fq"
        if layout == "PE":
            fastq_2.write_text("", encoding="utf-8")
            fastq_2_value = str(fastq_2)
        else:
            fastq_2_value = ""
        target = "H3K27me3" if peak_mode == "broad" else "CTCF"
        rendered.append(
            "\t".join(
                (
                    sample,
                    str(fastq_1),
                    fastq_2_value,
                    layout,
                    assay,
                    target,
                    peak_mode,
                    "hs",
                    "idx",
                    experiment,
                    biorep,
                )
            )
        )
    return "\n".join(rendered) + "\n"


def _base_config(tmp_path: Path) -> dict:
    return {
        "outdir": str(tmp_path / "results"),
        "use_control": False,
        "threads": 1,
        "trim": False,
        "multiqc": False,
        "stage4b": True,
        "stage5": False,
        "qc": {
            "blacklist_filter": False,
            "frip": False,
            "library_complexity": False,
            "nrf_pbc": False,
            "signal_tracks": False,
            "summary": False,
            "cuttag_fragment_size": False,
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
    *,
    sample_rows,
    reproducibility,
    idr=None,
    stage5=False,
):
    config = _base_config(tmp_path)
    config["stage5"] = stage5
    config["reproducibility"] = reproducibility
    if idr is not None:
        config["idr"] = idr
    _, config_path, _ = tmp_config(
        config=config,
        samples=_sample_sheet(tmp_path, sample_rows),
    )
    result = run_snakemake(
        config_path,
        extra_args=["--printshellcmds"],
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-5000:]
    return output, _scheduled_rules(output)


def test_custom_idr_threshold_rank_and_seed_reach_snakemake_shell_commands(
    tmp_path,
    tmp_config,
    run_snakemake,
):
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("S1", "atac", "narrow", "exp_atac", "1", "PE"),
            ("S2", "atac", "narrow", "exp_atac", "2", "PE"),
        ],
        reproducibility={
            "enabled": True,
            "consensus": {"enabled": False},
            "idr": {"atac_narrow": True},
        },
        idr={
            "threshold": 0.0125,
            "rank": "signal.value",
            "seed": 777,
        },
    )

    assert {
        "idr_true_replicates_narrow",
        "idr_split_pseudoreps_narrow",
        "idr_self_pseudoreps_narrow",
        "idr_pooled_pseudoreps_narrow",
    } <= rules
    assert "--idr-threshold 0.0125" in output
    assert "--rank signal.value" in output
    assert "--seed 777" in output


@pytest.mark.parametrize(
    (
        "assay",
        "peak_mode",
        "flag",
        "summary_rule",
        "shared_rules",
        "final_suffix",
    ),
    MODE_CASES,
)
def test_enabled_modes_schedule_current_unified_idr_rules_and_commands(
    tmp_path,
    tmp_config,
    run_snakemake,
    assay,
    peak_mode,
    flag,
    summary_rule,
    shared_rules,
    final_suffix,
):
    experiment = f"exp_{assay}_{peak_mode}"
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("S1", assay, peak_mode, experiment, "1", "PE"),
            ("S2", assay, peak_mode, experiment, "2", "PE"),
        ],
        reproducibility={
            "enabled": True,
            "consensus": {"enabled": False},
            "idr": {flag: True},
        },
    )

    assert shared_rules | {summary_rule} <= rules
    assert final_suffix in output
    assert "reproducibility_summary.tsv" in output
    assert f"--assay {assay}" in output
    input_type = "broadPeak" if peak_mode == "broad" else "narrowPeak"
    assert f"--input-file-type {input_type}" in output
    if peak_mode == "broad":
        assert "--broad" in output
        assert "--broad-cutoff" in output


@pytest.mark.parametrize(
    (
        "assay",
        "peak_mode",
        "flag",
        "summary_rule",
        "shared_rules",
        "_final_suffix",
    ),
    MODE_CASES,
)
def test_disabled_mode_flags_omit_idr_targets_instead_of_scheduling_generic_rules(
    tmp_path,
    tmp_config,
    run_snakemake,
    assay,
    peak_mode,
    flag,
    summary_rule,
    shared_rules,
    _final_suffix,
):
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("S1", assay, peak_mode, "exp1", "1", "PE"),
            ("S2", assay, peak_mode, "exp1", "2", "PE"),
        ],
        reproducibility={
            "enabled": True,
            "consensus": {"enabled": False},
            "idr": {flag: False},
        },
    )

    assert summary_rule not in rules
    assert shared_rules.isdisjoint(rules)
    assert "/06_reproducibility/idr/" not in output


def test_disabled_reproducibility_block_omits_idr_dag(
    tmp_path,
    tmp_config,
    run_snakemake,
):
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("A1", "atac", "narrow", "atac_exp", "1", "PE"),
            ("A2", "atac", "narrow", "atac_exp", "2", "PE"),
        ],
        reproducibility={
            "enabled": False,
            "idr": {"atac_narrow": True},
        },
    )

    assert NARROW_RULES.isdisjoint(rules)
    assert "idr_summary_atac_narrow" not in rules
    assert "/06_reproducibility/idr/" not in output


def test_cuttag_single_end_layout_uses_the_same_narrow_idr_graph(
    tmp_path,
    tmp_config,
    run_snakemake,
):
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("T1", "cuttag", "narrow", "cuttag_se", "1", "SE"),
            ("T2", "cuttag", "narrow", "cuttag_se", "2", "SE"),
        ],
        reproducibility={
            "enabled": True,
            "consensus": {"enabled": False},
            "idr": {"cuttag_narrow": True},
        },
    )

    assert NARROW_RULES | {"idr_summary_cuttag_narrow"} <= rules
    assert "cuttag.macs3.narrow.replicate_validated.idr.narrowPeak" in output


def test_mixed_atac_and_chipseq_schedules_only_the_enabled_atac_idr_namespace(
    tmp_path,
    tmp_config,
    run_snakemake,
):
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("A1", "atac", "narrow", "atac_exp", "1", "PE"),
            ("A2", "atac", "narrow", "atac_exp", "2", "PE"),
            ("C1", "chipseq", "narrow", "chip_exp", "1", "PE"),
            ("C2", "chipseq", "narrow", "chip_exp", "2", "PE"),
        ],
        reproducibility={
            "enabled": True,
            "consensus": {"enabled": False},
            "idr": {"atac_narrow": True},
        },
    )

    assert "idr_summary_atac_narrow" in rules
    assert "experiments/atac_exp/06_reproducibility/idr/" in output
    assert "experiments/chip_exp/06_reproducibility/idr/" not in output
    assert "experiments/chip_exp/06_idr/" not in output


def test_cuttag_idr_coexists_with_chipseq_narrow_idr(
    tmp_path,
    tmp_config,
    run_snakemake,
):
    output, rules = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        sample_rows=[
            ("T1", "cuttag", "narrow", "cuttag_exp", "1", "PE"),
            ("T2", "cuttag", "narrow", "cuttag_exp", "2", "PE"),
            ("C1", "chipseq", "narrow", "chip_exp", "1", "PE"),
            ("C2", "chipseq", "narrow", "chip_exp", "2", "PE"),
        ],
        reproducibility={
            "enabled": True,
            "consensus": {"enabled": False},
            "idr": {"cuttag_narrow": True},
        },
        stage5=True,
    )

    assert "idr_summary_cuttag_narrow" in rules
    assert "stage5b_summary" in rules
    assert "experiments/cuttag_exp/06_reproducibility/idr/" in output
    assert "experiments/chip_exp/06_idr/" in output
