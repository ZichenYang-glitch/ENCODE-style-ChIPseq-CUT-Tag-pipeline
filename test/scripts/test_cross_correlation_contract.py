"""Behavior contracts for SCI-6 input, interpretation and detail consumers."""

import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

from encode_pipeline.config.qc import validate_qc_config
from encode_pipeline.qc.contract import CROSS_CORRELATION_QC_V2_COLUMNS
from scripts import make_cross_correlation_tagalign as tags
from scripts import qc_common
from scripts.assemble_qc_summary import _cross_correlation

ROOT = Path(__file__).resolve().parents[2]


def invoke(name, *args, check=True):
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
        capture_output=True,
        text=True,
    )
    if check:
        assert result.returncode == 0, result.stderr
    return result


@pytest.mark.parametrize(
    "lengths,expected", [([50, 51], 50), ([51, 52], 52), ([50] * 500 + [100], 50)]
)
def test_first_500_estimate_and_all_input_variability(tmp_path, lengths, expected):
    path = tmp_path / "tags.tsv"
    path.write_text(
        "".join(f"chr1\t{i}\t{i + n}\tN\t1000\t+\n" for i, n in enumerate(lengths))
    )
    length, distribution, variable = tags._read_stats(path)
    assert length == expected
    assert sum(distribution.values()) == len(lengths)
    assert variable is (len(set(lengths)) > 1)


@pytest.mark.parametrize(
    "assay,target,requested,resolved,upper",
    [
        ("chipseq", "H3K27ac", "auto", "histone", 100),
        ("chipseq", "CTCF", "auto", "unknown", 100),
        ("chipseq", "CTCF", "tf", "tf", 60),
        ("chipseq", "CTCF", "histone", "histone", 100),
        ("atac", "CTCF", "tf", "unknown", 100),
        ("cuttag", "H3K27ac", "auto", "unknown", 100),
        ("mnase", "H3", "histone", "unknown", 100),
    ],
)
def test_target_resolution_and_window(assay, target, requested, resolved, upper):
    kind, _, _ = tags._resolve_target(assay, target, requested)
    assert kind == resolved
    assert tags._window(50, kind)[:2] == (-500, upper)
    assert tags._window(150, kind)[1] == 160


@pytest.mark.parametrize("value", ["tf", "histone", "auto"])
def test_config_target_type_is_stable(value):
    once = validate_qc_config({"cross_correlation_target_type": value})
    assert validate_qc_config(once) == once


@pytest.mark.parametrize("value", [[], {}, None, True, "TF", "bad"])
def test_invalid_target_type_is_a_validation_error(value):
    with pytest.raises(ValueError, match="cross_correlation_target_type"):
        validate_qc_config({"cross_correlation_target_type": value})


@pytest.mark.parametrize("layout,expected", [("PE", 1), ("SE", 2)])
def test_tagalign_real_script_respects_layout_and_policy(tmp_path, layout, expected):
    # Includes a marked duplicate: input analysis policy retains it here.
    source = tmp_path / "input.bam"
    source.write_text(
        "p\t1123\tchr1\t100\t60\t50M\t=\t250\t200\t*\t*\n"
        "p\t1171\tchr1\t250\t60\t50M\t=\t100\t-200\t*\t*\n"
        "supp\t2048\tchr1\t500\t60\t50M\t*\t0\t0\t*\t*\n"
    )
    tool = tmp_path / "samtools"
    tool.write_text('#!/bin/sh\ncat "$2"\n')
    tool.chmod(0o755)
    meta = tmp_path / "meta.json"
    out = tmp_path / "out.tagAlign"
    invoke(
        "make_cross_correlation_tagalign.py",
        "--analysis-bam",
        source,
        "--samtools",
        tool,
        "--layout",
        layout,
        "--sample",
        "S1",
        "--assay",
        "chipseq",
        "--target",
        "CTCF",
        "--target-type",
        "tf",
        "--tagalign",
        out,
        "--metadata",
        meta,
    )
    rows = out.read_text().splitlines()
    assert len(rows) == expected
    assert rows[0] == "chr1\t99\t149\tN\t1000\t+"
    metadata = json.loads(meta.read_text())
    assert metadata["input_stage"] == (
        "analysis_r1" if layout == "PE" else "analysis_se"
    )
    assert metadata["read_length"] == 50 and metadata["window_upper"] == 60


@pytest.fixture
def detail(tmp_path):
    def build(
        *,
        peak="200",
        scores="0.9",
        nsc="1.4",
        rsc="2",
        phantom="50",
        phantom_score="0.5",
        **overrides,
    ):
        metadata = {
            "schema_version": 1,
            "status": "ready",
            "assay": "chipseq",
            "layout": "SE",
            "input_stage": "analysis_se",
            "target_type": "tf",
            "read_length": 50,
            "read_length_method": "median_first_500_tag_spans_r_round",
            "variable_read_length": False,
            "window_lower": -500,
            "window_upper": 60,
        }
        metadata.update(overrides)
        meta = tmp_path / "meta.json"
        meta.write_text(json.dumps(metadata))
        raw = tmp_path / "S1.cc.qc"
        raw.write_text(
            f"S1\t1000\t{peak}\t{scores}\t{phantom}\t{phantom_score}\t1500\t0.1\t{nsc}\t{rsc}\t1\n"
        )
        out = tmp_path / "detail.tsv"
        invoke(
            "write_cross_correlation_qc.py",
            "--sample",
            "S1",
            "--assay",
            metadata["assay"],
            "--layout",
            "SE",
            "--metadata",
            meta,
            "--cc-qc",
            raw,
            "--output",
            out,
        )
        version, row = qc_common.read_cross_correlation_detail(str(out))
        assert version == "2" and tuple(row) == CROSS_CORRELATION_QC_V2_COLUMNS
        return out, row

    return build


@pytest.mark.parametrize(
    "overrides,status,hint",
    [
        ({}, "ok", "local_pass"),
        ({"rsc": "-1332.5"}, "low_rsc", "local_fail"),
        (
            {
                "peak": "50",
                "scores": "0.9265894939787",
                "phantom_score": "0.9265895",
                "rsc": "1",
            },
            "phantom_selected",
            "not_evaluated",
        ),
        (
            {"peak": "200,50", "scores": "0.9,0.9", "phantom_score": "0.9"},
            "phantom_selected",
            "not_evaluated",
        ),
        (
            {"peak": "200,50", "scores": "0.9,0.8", "phantom_score": "0.8"},
            "ok",
            "local_pass",
        ),
        ({"target_type": "unknown"}, "unknown_target", "not_evaluated"),
        ({"assay": "atac"}, "non_chip_diagnostic", "not_evaluated"),
        ({"variable_read_length": True}, "variable_read_length", "not_evaluated"),
        ({"rsc": "NaN"}, "missing_values", "not_evaluated"),
        ({"peak": "200,300", "scores": "0.9"}, "missing_values", "not_evaluated"),
    ],
)
def test_detail_writer_and_sample_consumer(detail, overrides, status, hint):
    path, row = detail(**overrides)
    assert row["execution_status"] == "reported"
    assert row["interpretation_status"] == status and row["threshold_hint"] == hint
    nsc, rsc, interpreted = _cross_correlation("", True, str(path))
    assert interpreted == status
    assert nsc == row["nsc"] and rsc == row["rsc"]
    if "rsc" in overrides and overrides["rsc"].startswith("-"):
        assert float(rsc) < 0


def test_no_data_is_not_an_applicability_claim(detail):
    path, row = detail(status="no_data")
    assert row["execution_status"] == row["interpretation_status"] == "no_data"
    assert row["nsc"] == row["rsc"] == "NA" and row["candidate_count"] == "0"
    assert _cross_correlation("", True, str(path)) == ("NA", "NA", "no_data")


@pytest.mark.parametrize(
    "key,value",
    [
        ("schema_version", "999"),
        ("rsc", "Inf"),
        ("candidate_count", "3"),
        ("threshold_hint", "bogus"),
        ("variable_read_length", "1"),
    ],
)
def test_corrupt_detail_is_rejected_by_both_consumers(detail, key, value, tmp_path):
    path, row = detail()
    row[key] = value
    qc_common.write_exact_tsv(path, CROSS_CORRELATION_QC_V2_COLUMNS, row)
    with pytest.raises(SystemExit):
        _cross_correlation("", True, str(path))
    result = invoke(
        "summarize_cross_correlation.py",
        "--input",
        path,
        "--output",
        tmp_path / "summary.tsv",
        check=False,
    )
    assert result.returncode != 0 and not (tmp_path / "summary.tsv").exists()


def test_project_summary_preserves_diagnostic_status(detail, tmp_path):
    path, row = detail(target_type="unknown", rsc="-2")
    out = tmp_path / "summary.tsv"
    invoke("summarize_cross_correlation.py", "--input", path, "--output", out)
    with out.open() as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert summary["interpretation_status"] == "unknown_target"
    assert summary["threshold_hint"] == "not_evaluated" and summary["rsc"] == "-2"


@pytest.mark.parametrize(
    "metadata,ok",
    [
        ({"status": "ready", "window_lower": -500, "window_upper": 60}, True),
        ({"status": "no_data", "window_lower": -500, "window_upper": 100}, False),
        ({"status": "ready", "window_lower": -500, "window_upper": False}, False),
        ({"status": "ready", "window_lower": 100, "window_upper": 60}, False),
        ({"status": "ready", "window_lower": -500, "window_upper": "60"}, False),
    ],
)
def test_window_cli_distinguishes_ready_and_unusable_metadata(tmp_path, metadata, ok):
    path = tmp_path / "meta.json"
    path.write_text(json.dumps({"schema_version": 1, **metadata}))
    result = invoke("cross_correlation_spp_window.py", "--metadata", path, check=False)
    assert (result.returncode == 0) is ok
    assert result.stdout == ("-500:60\n" if ok else "")


@pytest.mark.parametrize("version", [999, True, "1"])
def test_unknown_metadata_version_is_parse_failure_not_a_quality_result(
    tmp_path, version
):
    metadata = tmp_path / "meta.json"
    metadata.write_text(json.dumps({"schema_version": version, "status": "ready"}))
    detail = tmp_path / "detail.tsv"
    invoke(
        "write_cross_correlation_qc.py",
        "--sample",
        "S1",
        "--assay",
        "chipseq",
        "--layout",
        "SE",
        "--metadata",
        metadata,
        "--output",
        detail,
    )
    _, row = qc_common.read_cross_correlation_detail(str(detail))
    assert (
        row["execution_status"] == "parse_failed"
        and row["threshold_hint"] == "not_evaluated"
    )
    assert row["nsc"] == row["rsc"] == "NA"
    window = invoke(
        "cross_correlation_spp_window.py", "--metadata", metadata, check=False
    )
    assert window.returncode != 0 and not window.stdout
