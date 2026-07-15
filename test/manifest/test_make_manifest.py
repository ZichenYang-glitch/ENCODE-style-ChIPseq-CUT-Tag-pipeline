"""Native behavior tests for scripts/make_manifest.py."""

import csv
import io
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import scripts.make_manifest as make_manifest


def _make_samples(
    sample_rows=None,
    columns=None,
):
    """Return a samples.tsv string."""
    if columns is None:
        columns = [
            "sample",
            "fastq_1",
            "fastq_2",
            "layout",
            "assay",
            "target",
            "peak_mode",
            "genome",
            "bowtie2_index",
        ]
    if sample_rows is None:
        sample_rows = [
            [
                "S1",
                "R1.fq",
                "R2.fq",
                "PE",
                "chipseq",
                "CTCF",
                "narrow",
                "hs",
                "/opt/idx/ref",
            ]
        ]
    lines = ["\t".join(columns)]
    for row in sample_rows:
        lines.append("\t".join(str(v) for v in row))
    return "\n".join(lines) + "\n"


def _base_config():
    return {
        "use_control": False,
        "multiqc": True,
        "stage4b": True,
        "stage5": False,
        "qc": {"signal_tracks": True, "summary": True},
        "genome_resources": {
            "hs": {"effective_genome_size": "hs", "chrom_sizes": ""},
        },
    }


@pytest.fixture
def manifest_case(tmp_config):
    """Return a helper that builds a make_manifest case."""

    def _make(config_updates=None, sample_rows=None, columns=None, use_json=False):
        cfg = _base_config()
        if config_updates:
            for key, value in config_updates.items():
                if isinstance(value, dict) and isinstance(cfg.get(key), dict):
                    cfg[key].update(value)
                else:
                    cfg[key] = value

        workdir, config_path, samples_path = tmp_config(
            config=cfg,
            samples=_make_samples(sample_rows=sample_rows, columns=columns),
        )
        out = workdir / "manifest.tsv"
        return workdir, str(config_path), str(samples_path), str(out)

    return _make


def _run_manifest(config_path, out, use_json=False):
    """Run make_manifest.py directly and return parsed rows."""
    old_argv = sys.argv
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    stdout_capture = io.StringIO()
    stderr_capture = io.StringIO()
    sys.stdout = stdout_capture
    sys.stderr = stderr_capture
    sys.argv = ["make_manifest.py", "--output", out]
    if use_json:
        with open(config_path) as fh:
            sys.argv.extend(["--config-json", fh.read()])
    else:
        sys.argv.extend(["--config", config_path])

    try:
        make_manifest.main()
        returncode = 0
    except SystemExit as exc:
        returncode = exc.code if exc.code is not None else 0
    finally:
        sys.argv = old_argv
        sys.stdout = old_stdout
        sys.stderr = old_stderr

    class _Result:
        def __init__(self, returncode, stdout, stderr):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    result = _Result(returncode, stdout_capture.getvalue(), stderr_capture.getvalue())
    rows = []
    if result.returncode == 0 and os.path.exists(out):
        with open(out, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    return result, rows


# ---------------------------------------------------------------------------
# Default config
# ---------------------------------------------------------------------------


def test_default_expected_output_types(manifest_case):
    _, config_path, _, out = manifest_case()
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    types = {r["output_type"] for r in rows}
    expected = {
        "final_bam",
        "final_bai",
        "cpm_bigwig",
        "macs3_peak",
        "qc_summary",
        "macs3_fe_bdg",
        "macs3_ppois_bdg",
        "macs3_fe_bw",
        "macs3_ppois_bw",
        "stage3_qc_summary",
        "multiqc_report",
    }
    assert expected.issubset(types), f"Missing types: {expected - types}"


def test_signal_tracks_false(manifest_case):
    _, config_path, _, out = manifest_case(
        config_updates={"qc": {"signal_tracks": False, "summary": True}},
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0
    signal_types = {
        "macs3_fe_bdg",
        "macs3_ppois_bdg",
        "macs3_fe_bw",
        "macs3_ppois_bw",
    }
    signal_rows = [r for r in rows if r["output_type"] in signal_types]
    assert len(signal_rows) == len(signal_types)
    assert {r["output_type"] for r in signal_rows} == signal_types
    assert all(r["status"] == "not_applicable" for r in signal_rows)


def test_chrom_sizes_empty_bw_not_applicable(manifest_case):
    _, config_path, _, out = manifest_case()
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0
    bw_rows = [r for r in rows if r["output_type"] in ("macs3_fe_bw", "macs3_ppois_bw")]
    bdg_rows = [
        r for r in rows if r["output_type"] in ("macs3_fe_bdg", "macs3_ppois_bdg")
    ]
    assert len(bw_rows) == 2
    assert len(bdg_rows) == 2
    assert {r["output_type"] for r in bw_rows} == {"macs3_fe_bw", "macs3_ppois_bw"}
    assert {r["output_type"] for r in bdg_rows} == {
        "macs3_fe_bdg",
        "macs3_ppois_bdg",
    }
    assert all(r["status"] == "not_applicable" for r in bw_rows)
    assert all(r["status"] in ("present", "missing") for r in bdg_rows)


def test_disabled_chipseq_narrow_idr_omits_manifest_rows(manifest_case):
    columns = [
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
    sample_rows = [
        [
            "S1",
            "R1.fq",
            "R2.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "1",
        ],
        [
            "S2",
            "R2.fq",
            "R3.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "2",
        ],
    ]
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0
    idr_types = {r["output_type"] for r in rows if r["output_type"].startswith("idr_")}
    assert idr_types == set()


def test_missing_file_status(manifest_case):
    _, config_path, _, out = manifest_case()
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0
    check_rows = [r for r in rows if r["status"] != "not_applicable"]
    assert check_rows
    assert all(r["status"] == "missing" for r in check_rows)


def test_deterministic_ordering(manifest_case):
    sample_rows = [
        [
            "SB",
            "R1.fq",
            "R2.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
        ],
        [
            "SA",
            "R1.fq",
            "R2.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
        ],
    ]
    _, config_path, _, out1 = manifest_case(sample_rows=sample_rows)
    _, config_path2, _, out2 = manifest_case(sample_rows=sample_rows)

    result1, rows1 = _run_manifest(config_path, out1)
    assert result1.returncode == 0
    first_text = Path(out1).read_text(encoding="utf-8")

    result2, _ = _run_manifest(config_path2, out2)
    assert result2.returncode == 0
    assert first_text == Path(out2).read_text(encoding="utf-8")
    sample_order = list(
        dict.fromkeys(row["sample_id"] for row in rows1 if row["sample_id"])
    )
    assert sample_order == ["SA", "SB"]


def test_no_crlf(manifest_case):
    _, config_path, _, out = manifest_case()
    result, _ = _run_manifest(config_path, out)
    assert result.returncode == 0
    assert b"\r\n" not in open(out, "rb").read()


def test_strict_flag_fails_on_missing(manifest_case):
    _, config_path, _, out = manifest_case()
    old_argv = sys.argv
    sys.argv = [
        "make_manifest.py",
        "--config",
        config_path,
        "--output",
        out,
        "--strict",
    ]
    try:
        make_manifest.main()
        returncode = 0
    except SystemExit as exc:
        returncode = exc.code if exc.code is not None else 0
    finally:
        sys.argv = old_argv
    assert returncode != 0


# ---------------------------------------------------------------------------
# Config via JSON
# ---------------------------------------------------------------------------


def test_custom_outdir_config_json(manifest_case):
    custom_outdir = "/tmp/custom_results_phase2b"
    cfg = {
        "samples": "{samples}",
        "outdir": custom_outdir,
        "use_control": False,
        "multiqc": True,
        "stage4b": True,
        "stage5": False,
        "qc": {"signal_tracks": True, "summary": True},
        "genome_resources": {"hs": {"effective_genome_size": "hs", "chrom_sizes": ""}},
    }
    workdir, _, samples_path, out = manifest_case()
    cfg["samples"] = samples_path
    config_path = workdir / "config.json"
    with open(config_path, "w") as fh:
        json.dump(cfg, fh)

    result, rows = _run_manifest(str(config_path), out, use_json=True)
    assert result.returncode == 0, result.stderr[-200:]
    non_na = [r for r in rows if r["status"] != "not_applicable"]
    assert non_na
    assert all(r["path"].startswith(custom_outdir) for r in non_na)


# ---------------------------------------------------------------------------
# Sample-sheet variations
# ---------------------------------------------------------------------------


def test_missing_experiment_defaults(manifest_case):
    columns = [
        "sample",
        "fastq_1",
        "fastq_2",
        "layout",
        "assay",
        "target",
        "peak_mode",
        "genome",
        "bowtie2_index",
        "biological_replicate",
    ]
    sample_rows = [
        [
            "S1",
            "R1.fq",
            "R2.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "1",
        ],
        [
            "S2",
            "R1b.fq",
            "R2b.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "2",
        ],
    ]
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    exp_rows = [r for r in rows if r["experiment_id"]]
    assert len(exp_rows) == 0


def test_chrom_sizes_configured_bw_present(manifest_case, tmp_path):
    chrom_sizes = tmp_path / "chr.sizes"
    chrom_sizes.write_text("", encoding="utf-8")
    _, config_path, _, out = manifest_case(
        config_updates={
            "genome_resources": {
                "hs": {"effective_genome_size": "hs", "chrom_sizes": str(chrom_sizes)},
            },
        },
    )
    # make_manifest only checks existence when chrom_sizes is non-empty
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0
    bw_rows = [r for r in rows if r["output_type"] in ("macs3_fe_bw", "macs3_ppois_bw")]
    assert len(bw_rows) == 2
    assert all(r["status"] != "not_applicable" for r in bw_rows)


# ---------------------------------------------------------------------------
# Technical / biological replicate gating
# ---------------------------------------------------------------------------


def test_techrep_biorep_gating(manifest_case):
    columns = [
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
        "technical_replicate",
    ]
    sample_rows = [
        [
            "T1a",
            "R1.fq",
            "R1b.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "1",
            "1",
        ],
        [
            "T1b",
            "R2.fq",
            "R2b.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "1",
            "2",
        ],
    ]
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    types = {r["output_type"] for r in rows}
    assert "biorep1_final_bam" in types
    assert "pooled_final_bam" not in types
    assert "pooled_macs3_peak" not in types


def test_idr_with_tech_reps(manifest_case):
    columns = [
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
        "technical_replicate",
    ]
    sample_rows = [
        [
            "T1a",
            "R1.fq",
            "R1b.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "1",
            "1",
        ],
        [
            "T1b",
            "R2.fq",
            "R2b.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "1",
            "2",
        ],
        [
            "T2",
            "R3.fq",
            "R3b.fq",
            "PE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx/ref",
            "EXP1",
            "2",
            "1",
        ],
    ]
    _, config_path, _, out = manifest_case(
        config_updates={
            "stage5": True,
            "idr": {"seed": 42, "threshold": 0.05, "rank": "p.value"},
        },
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    types = {r["output_type"] for r in rows}
    expected_idr = {"idr_conservative", "idr_optimal", "idr_reproducibility_summary"}
    assert expected_idr.issubset(types), f"Missing IDR types: {expected_idr - types}"


# ---------------------------------------------------------------------------
# Reproducibility mode and final-selection contracts
# ---------------------------------------------------------------------------


REPRODUCIBILITY_COLUMNS = [
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


def _reproducibility_samples(assay, peak_mode, *, layout="PE", experiment="EXP1"):
    target = "H3K27me3" if peak_mode == "broad" else "CTCF"
    return [
        [
            f"{experiment}_1",
            "R1.fq",
            "R2.fq" if layout == "PE" else "",
            layout,
            assay,
            target,
            peak_mode,
            "hs",
            "/opt/idx/ref",
            experiment,
            "1",
        ],
        [
            f"{experiment}_2",
            "R3.fq",
            "R4.fq" if layout == "PE" else "",
            layout,
            assay,
            target,
            peak_mode,
            "hs",
            "/opt/idx/ref",
            experiment,
            "2",
        ],
    ]


def _run_reproducibility_manifest(
    manifest_case,
    *,
    assay,
    peak_mode,
    config_updates,
    layout="PE",
):
    _, config_path, _, out = manifest_case(
        config_updates=config_updates,
        sample_rows=_reproducibility_samples(assay, peak_mode, layout=layout),
        columns=REPRODUCIBILITY_COLUMNS,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-500:]
    return [row for row in rows if "/06_reproducibility/" in row["path"]]


@pytest.mark.parametrize(
    (
        "assay",
        "peak_mode",
        "flag",
        "consensus_peak_type",
        "consensus_summary_type",
        "idr_final_type",
        "idr_summary_type",
        "final_suffix",
    ),
    [
        pytest.param(
            "atac",
            "narrow",
            "atac_narrow",
            "atac_macs3_narrow_consensus_peak",
            "atac_macs3_narrow_consensus_summary",
            "atac_macs3_narrow_idr_final_peak",
            "atac_macs3_narrow_idr_summary",
            "atac.macs3.narrow.replicate_validated.idr.narrowPeak",
            id="atac-narrow",
        ),
        pytest.param(
            "cuttag",
            "narrow",
            "cuttag_narrow",
            "cuttag_macs3_narrow_consensus_peak",
            "cuttag_macs3_narrow_consensus_summary",
            "cuttag_macs3_narrow_idr_final_peak",
            "cuttag_macs3_narrow_idr_summary",
            "cuttag.macs3.narrow.replicate_validated.idr.narrowPeak",
            id="cuttag-narrow",
        ),
        pytest.param(
            "chipseq",
            "broad",
            "chipseq_broad_experimental",
            "chipseq_macs3_broad_consensus_peak",
            "chipseq_macs3_broad_consensus_summary",
            "chipseq_macs3_broad_idr_final_peak",
            "chipseq_macs3_broad_idr_summary",
            "chipseq.macs3.broad.replicate_validated.idr.broadPeak",
            id="chipseq-broad",
        ),
        pytest.param(
            "cuttag",
            "broad",
            "cuttag_broad_experimental",
            "cuttag_macs3_broad_consensus_peak",
            "cuttag_macs3_broad_consensus_summary",
            "cuttag_macs3_broad_idr_final_peak",
            "cuttag_macs3_broad_idr_summary",
            "cuttag.macs3.broad.replicate_validated.idr.broadPeak",
            id="cuttag-broad",
        ),
    ],
)
def test_enabled_idr_modes_emit_exact_rows_and_supersede_consensus_final(
    manifest_case,
    assay,
    peak_mode,
    flag,
    consensus_peak_type,
    consensus_summary_type,
    idr_final_type,
    idr_summary_type,
    final_suffix,
):
    rows = _run_reproducibility_manifest(
        manifest_case,
        assay=assay,
        peak_mode=peak_mode,
        config_updates={
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {flag: True},
            }
        },
    )
    by_type = {row["output_type"]: row for row in rows}

    assert len(rows) == len(by_type)
    assert set(by_type) == {
        consensus_peak_type,
        consensus_summary_type,
        idr_final_type,
        idr_summary_type,
    }
    assert by_type[idr_final_type]["path"].endswith(f"EXP1.{final_suffix}")
    assert by_type[idr_final_type]["method"] == "idr"
    assert by_type[idr_summary_type]["path"].endswith("reproducibility_summary.tsv")
    assert all(row["status"] == "missing" for row in rows)


@pytest.mark.parametrize(
    (
        "assay",
        "peak_mode",
        "flag",
        "consensus_peak_type",
        "consensus_summary_type",
        "consensus_final_type",
    ),
    [
        pytest.param(
            "atac",
            "narrow",
            "atac_narrow",
            "atac_macs3_narrow_consensus_peak",
            "atac_macs3_narrow_consensus_summary",
            None,
            id="atac-narrow",
        ),
        pytest.param(
            "cuttag",
            "narrow",
            "cuttag_narrow",
            "cuttag_macs3_narrow_consensus_peak",
            "cuttag_macs3_narrow_consensus_summary",
            "cuttag_macs3_narrow_consensus_final_peak",
            id="cuttag-narrow",
        ),
        pytest.param(
            "chipseq",
            "broad",
            "chipseq_broad_experimental",
            "chipseq_macs3_broad_consensus_peak",
            "chipseq_macs3_broad_consensus_summary",
            "chipseq_macs3_broad_consensus_final_peak",
            id="chipseq-broad",
        ),
        pytest.param(
            "cuttag",
            "broad",
            "cuttag_broad_experimental",
            "cuttag_macs3_broad_consensus_peak",
            "cuttag_macs3_broad_consensus_summary",
            "cuttag_macs3_broad_consensus_final_peak",
            id="cuttag-broad",
        ),
    ],
)
def test_disabled_idr_modes_are_omitted_and_consensus_remains_primary(
    manifest_case,
    assay,
    peak_mode,
    flag,
    consensus_peak_type,
    consensus_summary_type,
    consensus_final_type,
):
    rows = _run_reproducibility_manifest(
        manifest_case,
        assay=assay,
        peak_mode=peak_mode,
        config_updates={
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {flag: False},
            }
        },
    )
    types = {row["output_type"] for row in rows}
    expected = {consensus_peak_type, consensus_summary_type}
    if consensus_final_type is not None:
        expected.add(consensus_final_type)

    assert len(rows) == len(types)
    assert types == expected
    assert not any("idr" in output_type for output_type in types)
    assert all(row["status"] != "not_applicable" for row in rows)


@pytest.mark.parametrize(
    "config_updates",
    [
        pytest.param(
            {"reproducibility": {"enabled": False}},
            id="reproducibility-disabled",
        ),
        pytest.param(
            {
                "stage4b": False,
                "reproducibility": {
                    "enabled": True,
                    "consensus": {"enabled": True},
                },
            },
            id="replicate-stage-disabled",
        ),
    ],
)
def test_global_reproducibility_gates_omit_all_rows(manifest_case, config_updates):
    rows = _run_reproducibility_manifest(
        manifest_case,
        assay="cuttag",
        peak_mode="narrow",
        config_updates=config_updates,
    )

    assert rows == []


def test_seacr_consensus_rows_are_mode_specific_and_idr_free(manifest_case):
    rows = _run_reproducibility_manifest(
        manifest_case,
        assay="cuttag",
        peak_mode="narrow",
        config_updates={
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
            },
            "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
        },
    )
    types = {row["output_type"] for row in rows}

    assert len(rows) == len(types)
    assert types == {
        "cuttag_macs3_narrow_consensus_peak",
        "cuttag_macs3_narrow_consensus_summary",
        "cuttag_macs3_narrow_consensus_final_peak",
        "cuttag_seacr_consensus_peak",
        "cuttag_seacr_consensus_summary",
        "cuttag_seacr_consensus_final_peak",
    }
    assert not any("idr" in output_type for output_type in types)


def test_chipseq_narrow_never_uses_the_new_idr_namespace(manifest_case):
    rows = _run_reproducibility_manifest(
        manifest_case,
        assay="chipseq",
        peak_mode="narrow",
        config_updates={
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
            }
        },
    )
    types = {row["output_type"] for row in rows}

    assert len(rows) == len(types)
    assert types == {
        "chipseq_macs3_narrow_consensus_peak",
        "chipseq_macs3_narrow_consensus_summary",
    }
    assert "chipseq_macs3_narrow_idr_final_peak" not in types
    assert not any("replicate_validated.idr" in row["path"] for row in rows)


# ---------------------------------------------------------------------------
# MNase-specific rows
# ---------------------------------------------------------------------------


def _mnase_samples(experiment=None, bioreps=None):
    columns = [
        "sample",
        "fastq_1",
        "fastq_2",
        "layout",
        "assay",
        "target",
        "peak_mode",
        "genome",
        "bowtie2_index",
    ]
    rows = [
        [
            "M1",
            "R1.fq",
            "R2.fq",
            "PE",
            "mnase",
            "H3",
            "nucleosome",
            "hs",
            "/opt/idx/ref",
        ]
    ]
    if experiment:
        columns.extend(["experiment", "biological_replicate"])
        for idx, br in enumerate(bioreps or ["1"]):
            row = [
                f"M{idx + 1}",
                f"R{idx + 1}.fq",
                f"R{idx + 1}b.fq",
                "PE",
                "mnase",
                "H3",
                "nucleosome",
                "hs",
                "/opt/idx/ref",
                experiment,
                br,
            ]
            rows.append(row)
        rows = rows[1:]
    return columns, rows


def test_mnase_sample_rows(manifest_case):
    columns, sample_rows = _mnase_samples()
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    types = {r["output_type"] for r in rows}
    expected = {
        "mnase_mono_bam",
        "mnase_mono_bai",
        "mnase_dyad_bigwig",
        "mnase_mono_bigwig",
        "mnase_sub_bam",
        "mnase_sub_bai",
        "mnase_di_bam",
        "mnase_di_bai",
        "mnase_qc_summary",
    }
    assert expected.issubset(types), f"Missing MNase types: {expected - types}"


def test_mnase_excludes_peak_outputs(manifest_case):
    columns, sample_rows = _mnase_samples()
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    peak_only = {
        "macs3_peak",
        "qc_summary",
        "macs3_fe_bdg",
        "macs3_ppois_bdg",
        "macs3_fe_bw",
        "macs3_ppois_bw",
    }
    for pt in peak_only:
        matching = [r for r in rows if r["output_type"] == pt]
        if matching:
            assert matching[0]["status"] == "not_applicable", (
                f"{pt} status={matching[0]['status']}, expected not_applicable"
            )


def test_mnase_pooled_rows(manifest_case):
    columns, sample_rows = _mnase_samples(
        experiment="EXP1",
        bioreps=["1", "2"],
    )
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    types = {r["output_type"] for r in rows}
    expected = {
        "pooled_mnase_mono_bam",
        "pooled_mnase_mono_bai",
        "pooled_mnase_dyad_bigwig",
        "pooled_mnase_mono_bigwig",
    }
    assert expected.issubset(types), f"Missing pooled MNase types: {expected - types}"
    excluded = {
        "pooled_macs3_peak",
        "pooled_qc_summary",
        "pooled_fe_bdg",
        "pooled_ppois_bdg",
        "pooled_fe_bw",
        "pooled_ppois_bw",
    }
    for pt in excluded:
        matching = [r for r in rows if r["output_type"] == pt]
        if matching:
            assert matching[0]["status"] == "not_applicable", (
                f"{pt} status={matching[0]['status']}, expected not_applicable"
            )


def test_mnase_pure_no_stage3_qc_summary(manifest_case):
    columns, sample_rows = _mnase_samples()
    _, config_path, _, out = manifest_case(
        sample_rows=sample_rows,
        columns=columns,
    )
    result, rows = _run_manifest(config_path, out)
    assert result.returncode == 0, result.stderr[-200:]
    matching = [r for r in rows if r["output_type"] == "stage3_qc_summary"]
    assert len(matching) == 1
    assert matching[0]["status"] == "not_applicable"


# ---------------------------------------------------------------------------
# Dependency targets
# ---------------------------------------------------------------------------


def test_manifest_dependency_targets_include_downstream(manifest_case, run_snakemake):
    columns = [
        "sample",
        "fastq_1",
        "layout",
        "assay",
        "target",
        "peak_mode",
        "genome",
        "bowtie2_index",
        "experiment",
        "biological_replicate",
    ]
    workdir, config_path, samples_path, _ = manifest_case(
        sample_rows=[
            [
                "S1",
                "S1_R1.fq",
                "SE",
                "chipseq",
                "CTCF",
                "narrow",
                "hs",
                "/opt/idx",
                "EXP1",
                "1",
            ]
        ],
        columns=columns,
    )
    # Placeholder FASTQ referenced by absolute path in the sample sheet
    fastq_path = workdir / "S1_R1.fq"
    fastq_path.write_text("", encoding="utf-8")
    sample_rows = [
        [
            "S1",
            str(fastq_path),
            "SE",
            "chipseq",
            "CTCF",
            "narrow",
            "hs",
            "/opt/idx",
            "EXP1",
            "1",
        ]
    ]
    samples_path = Path(samples_path)
    samples_path.write_text(
        "\t".join(columns) + "\n" + "\t".join(sample_rows[0]) + "\n", encoding="utf-8"
    )

    result = run_snakemake(config_path, quiet=False)
    assert result.returncode == 0, f"Dry-run failed:\n{result.stderr}"
    assert "result_manifest" in (result.stdout + result.stderr)
