"""Snakemake dry-run contracts for consensus scheduling and final selection."""

from __future__ import annotations

import pytest


pytestmark = pytest.mark.full_main


SAMPLE_COLUMNS = [
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
]


def _replicates(assay, peak_mode, experiment, *, count=2, layout="PE"):
    return [
        {
            "sample": f"{experiment}_{index}",
            "layout": layout,
            "assay": assay,
            "peak_mode": peak_mode,
            "experiment": experiment,
            "biological_replicate": str(index),
        }
        for index in range(1, count + 1)
    ]


def _config(
    tmp_path,
    *,
    reproducibility=True,
    consensus=True,
    seacr=None,
    stage5=False,
    cuttag_idr=None,
):
    reproducibility_config = {
        "enabled": reproducibility,
        "consensus": {"enabled": consensus},
    }
    if cuttag_idr is not None:
        reproducibility_config["idr"] = {"cuttag_narrow": cuttag_idr}

    config = {
        "outdir": str(tmp_path / "results"),
        "use_control": False,
        "threads": 1,
        "stage4b": True,
        "stage5": stage5,
        "reproducibility": reproducibility_config,
    }
    if seacr is not None:
        config["cuttag"] = {
            "seacr": {
                "enabled": seacr,
                "mode": "stringent",
            }
        }
    return config


def _macs3_consensus_peak(experiment, assay, peak_mode):
    suffix = "narrowPeak" if peak_mode == "narrow" else "broadPeak"
    return (
        f"06_reproducibility/consensus/{experiment}.{assay}.macs3."
        f"{peak_mode}.consensus.{suffix}"
    )


def _macs3_consensus_summary(experiment, assay, peak_mode):
    return (
        f"06_reproducibility/consensus/{experiment}.{assay}.macs3."
        f"{peak_mode}.consensus.summary.tsv"
    )


def _macs3_consensus_final(experiment, assay, peak_mode):
    suffix = "narrowPeak" if peak_mode == "narrow" else "broadPeak"
    return (
        f"06_reproducibility/final/{experiment}.{assay}.macs3.{peak_mode}."
        f"replicate_validated.consensus.{suffix}"
    )


def _seacr_consensus_peak(experiment):
    return (
        f"06_reproducibility/consensus/{experiment}.cuttag.seacr.stringent."
        "consensus.bed"
    )


def _seacr_consensus_summary(experiment):
    return (
        f"06_reproducibility/consensus/{experiment}.cuttag.seacr.stringent."
        "consensus.summary.tsv"
    )


def _seacr_consensus_final(experiment):
    return (
        f"06_reproducibility/final/{experiment}.cuttag.seacr.stringent."
        "replicate_validated.consensus.bed"
    )


@pytest.fixture
def consensus_dry_run(tmp_path, tmp_config, run_snakemake, monkeypatch):
    """Run one isolated, non-executing consensus DAG case."""

    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / ".cache"))
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")

    def _run(sample_specs, config, *, printshellcmds=False):
        read1 = tmp_path / "reads_R1.fastq"
        read2 = tmp_path / "reads_R2.fastq"
        read1.write_text("", encoding="utf-8")
        read2.write_text("", encoding="utf-8")

        rows = ["\t".join(SAMPLE_COLUMNS)]
        for spec in sample_specs:
            row = {
                "sample": spec["sample"],
                "fastq_1": str(read1),
                "fastq_2": str(read2) if spec["layout"] == "PE" else "",
                "layout": spec["layout"],
                "assay": spec["assay"],
                "target": "T",
                "peak_mode": spec["peak_mode"],
                "genome": "hs",
                "bowtie2_index": "idx",
                "experiment": spec["experiment"],
                "biological_replicate": spec["biological_replicate"],
            }
            rows.append("\t".join(row[column] for column in SAMPLE_COLUMNS))

        _, config_path, _ = tmp_config(
            config=config,
            samples="\n".join(rows) + "\n",
        )
        extra_args = ["--printshellcmds"] if printshellcmds else None
        result = run_snakemake(
            config_path,
            extra_args=extra_args,
            quiet=False,
        )
        output = result.stdout + result.stderr
        assert result.returncode == 0, (
            "Consensus DAG dry-run failed:\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
        return output

    return _run


@pytest.mark.parametrize(
    ("reproducibility", "consensus"),
    [(False, True), (True, False)],
    ids=["reproducibility-disabled", "consensus-disabled"],
)
def test_consensus_respects_global_gates(
    consensus_dry_run, tmp_path, reproducibility, consensus
):
    output = consensus_dry_run(
        _replicates("chipseq", "narrow", "exp_cs"),
        _config(
            tmp_path,
            reproducibility=reproducibility,
            consensus=consensus,
        ),
    )
    assert "06_reproducibility/consensus/" not in output


@pytest.mark.parametrize(
    ("assay", "peak_mode", "experiment", "expects_final"),
    [
        ("chipseq", "narrow", "exp_cs", False),
        ("chipseq", "broad", "exp_cb", True),
        ("cuttag", "narrow", "exp_ct", True),
        ("cuttag", "broad", "exp_ctb", True),
        ("atac", "narrow", "exp_at", False),
    ],
    ids=[
        "chipseq-narrow",
        "chipseq-broad",
        "cuttag-narrow",
        "cuttag-broad",
        "atac-narrow",
    ],
)
def test_macs3_consensus_schedules_mode_outputs(
    consensus_dry_run,
    tmp_path,
    assay,
    peak_mode,
    experiment,
    expects_final,
):
    printshellcmds = assay == "chipseq" and peak_mode == "narrow"
    output = consensus_dry_run(
        _replicates(assay, peak_mode, experiment),
        _config(tmp_path),
        printshellcmds=printshellcmds,
    )

    assert _macs3_consensus_peak(experiment, assay, peak_mode) in output
    assert _macs3_consensus_summary(experiment, assay, peak_mode) in output
    final_path = _macs3_consensus_final(experiment, assay, peak_mode)
    if expects_final:
        assert final_path in output
    else:
        assert final_path not in output

    if printshellcmds:
        assert (
            "06_reproducibility/consensus/biorep_peaks/"
            "exp_cs.biorep1.chipseq.macs3.narrow_peaks.narrowPeak"
        ) in output
        assert "06_idr/" not in output

    if assay == "cuttag" and peak_mode == "narrow":
        assert ".cuttag.seacr." not in output
        assert "consensus_compute_seacr" not in output


def test_consensus_requires_multiple_bioreps(consensus_dry_run, tmp_path):
    output = consensus_dry_run(
        _replicates("chipseq", "narrow", "exp_one", count=1),
        _config(tmp_path),
    )
    assert "06_reproducibility/consensus/" not in output


@pytest.mark.parametrize(
    "final_policy",
    ["chipseq-idr", "cuttag-idr", "cuttag-consensus"],
)
def test_idr_takes_precedence_only_for_the_final_peak(
    consensus_dry_run, tmp_path, final_policy
):
    if final_policy == "chipseq-idr":
        assay = "chipseq"
        experiment = "exp_cs"
        config = _config(tmp_path, stage5=True)
        expected_final = "06_idr/final/conservative.narrowPeak"
        rejected_final = _macs3_consensus_final(experiment, assay, "narrow")
    elif final_policy == "cuttag-idr":
        assay = "cuttag"
        experiment = "exp_ct"
        config = _config(tmp_path, cuttag_idr=True)
        expected_final = (
            "06_reproducibility/final/exp_ct.cuttag.macs3.narrow."
            "replicate_validated.idr.narrowPeak"
        )
        rejected_final = _macs3_consensus_final(experiment, assay, "narrow")
    else:
        assay = "cuttag"
        experiment = "exp_ct"
        config = _config(tmp_path, cuttag_idr=False)
        expected_final = _macs3_consensus_final(experiment, assay, "narrow")
        rejected_final = (
            "06_reproducibility/final/exp_ct.cuttag.macs3.narrow."
            "replicate_validated.idr.narrowPeak"
        )

    output = consensus_dry_run(
        _replicates(assay, "narrow", experiment),
        config,
    )
    assert _macs3_consensus_peak(experiment, assay, "narrow") in output
    assert _macs3_consensus_summary(experiment, assay, "narrow") in output
    assert expected_final in output
    assert rejected_final not in output


@pytest.mark.parametrize(
    ("seacr", "reproducibility", "consensus"),
    [(False, True, True), (True, False, True), (True, True, False)],
    ids=["seacr-disabled", "reproducibility-disabled", "consensus-disabled"],
)
def test_seacr_consensus_respects_independent_gates(
    consensus_dry_run,
    tmp_path,
    seacr,
    reproducibility,
    consensus,
):
    output = consensus_dry_run(
        _replicates("cuttag", "narrow", "exp_ct"),
        _config(
            tmp_path,
            seacr=seacr,
            reproducibility=reproducibility,
            consensus=consensus,
        ),
    )
    assert _seacr_consensus_peak("exp_ct") not in output


@pytest.mark.parametrize(
    ("assay", "layout", "replicate_count", "experiment"),
    [
        ("cuttag", "PE", 1, "exp_one"),
        ("cuttag", "SE", 2, "exp_ctse"),
        ("chipseq", "PE", 2, "exp_cs"),
    ],
    ids=["one-biorep", "single-end", "chipseq"],
)
def test_seacr_consensus_requires_eligible_cuttag_pe_replicates(
    consensus_dry_run,
    tmp_path,
    assay,
    layout,
    replicate_count,
    experiment,
):
    output = consensus_dry_run(
        _replicates(
            assay,
            "narrow",
            experiment,
            count=replicate_count,
            layout=layout,
        ),
        _config(tmp_path, seacr=True),
    )
    assert f"{experiment}.cuttag.seacr.stringent.consensus.bed" not in output


def test_seacr_consensus_schedules_outputs_and_command(consensus_dry_run, tmp_path):
    output = consensus_dry_run(
        _replicates("cuttag", "narrow", "exp_ct"),
        _config(tmp_path, seacr=True),
        printshellcmds=True,
    )

    assert _seacr_consensus_peak("exp_ct") in output
    assert _seacr_consensus_summary("exp_ct") in output
    assert _seacr_consensus_final("exp_ct") in output
    assert "--format bed" in output
    assert "--caller seacr" in output
    assert "--final-method consensus" in output

    assert "seacr_bedgraph" in output
    assert "seacr_call" in output
    assert "04_peaks_seacr" in output
    assert _macs3_consensus_peak("exp_ct", "cuttag", "narrow") in output

    lowered = output.lower()
    for forbidden in ("seacr.idr", "seacr_idr", "seacr_experimental"):
        assert forbidden not in lowered
