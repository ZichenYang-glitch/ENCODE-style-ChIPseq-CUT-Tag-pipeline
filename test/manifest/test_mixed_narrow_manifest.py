"""Manifest rows for explicit mixed narrow-peak routes.

Route selection stays in encode_pipeline.samples.mixed_replicates. These tests
supply the route object and check the manifest records it.
"""

from types import SimpleNamespace

import pytest

from encode_pipeline.manifest.make import (
    MANIFEST_COLUMNS,
    _build_idr_rows,
    _build_reproducibility_rows,
    _resolve_mixed_replicate_routes,
    build_manifest_rows,
)


def _route(method, n_bioreps, run_macs3_consensus):
    return SimpleNamespace(
        method=method,
        n_bioreps=n_bioreps,
        run_macs3_consensus=run_macs3_consensus,
    )


def _sample(experiment, assay, peak_mode, biorep, layout="PE"):
    return {
        "id": f"{experiment}_{biorep}",
        "role": "treatment",
        "experiment": experiment,
        "assay": assay,
        "peak_mode": peak_mode,
        "biological_replicate": biorep,
        "target": "CTCF",
        "genome": "hs",
        "layout": layout,
    }


def _groups():
    samples = []
    samples.extend(_sample("GRP2", "chipseq", "narrow", biorep) for biorep in (1, 2))
    samples.extend(
        _sample("GRP3", "chipseq", "narrow", biorep) for biorep in (1, 2, 3)
    )
    samples.extend(_sample("ATAC3", "atac", "narrow", biorep) for biorep in (1, 2, 3))
    samples.extend(_sample("CT2", "cuttag", "narrow", biorep) for biorep in (1, 2))
    samples.extend(_sample("CT3", "cuttag", "narrow", biorep) for biorep in (1, 2, 3))
    return samples


def _config():
    return {
        "replicate_analysis": True,
        "reproducibility": {
            "enabled": True,
            "consensus": {"enabled": True, "min_replicates": 2},
            "idr": {
                "atac_narrow": True,
                "cuttag_narrow": True,
            },
        },
        "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
        "chipseq_idr": True,
    }


def _flags():
    return {
        "chipseq_narrow": True,
        "atac_narrow": True,
        "cuttag_narrow": True,
        "chipseq_broad": False,
        "cuttag_broad": False,
    }


def _pairs(rows):
    return {(row["experiment_id"], row["output_type"]) for row in rows}


def _repro_rows(routes):
    return _build_reproducibility_rows(
        _groups(), _config(), "/tmp/results", _flags(), routes
    )


def test_manifest_keeps_ten_columns():
    assert len(MANIFEST_COLUMNS) == 10
    assert MANIFEST_COLUMNS == [
        "sample_id",
        "experiment_id",
        "assay",
        "target",
        "genome",
        "output_type",
        "method",
        "path",
        "status",
        "qc_flag",
    ]


def test_strict_and_unresolved_routes_match_and_omit_new_finals():
    explicit = _repro_rows({})
    resolved = _build_reproducibility_rows(
        _groups(), _config(), "/tmp/results", _flags()
    )
    assert explicit == resolved
    pairs = _pairs(explicit)
    assert ("GRP3", "chipseq_macs3_narrow_consensus_final_peak") not in pairs
    assert ("ATAC3", "atac_macs3_narrow_consensus_final_peak") not in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_peak") in pairs
    assert ("CT2", "cuttag_macs3_narrow_consensus_peak") in pairs
    assert ("CT2", "cuttag_seacr_consensus_final_peak") in pairs
    assert ("CT3", "cuttag_seacr_consensus_final_peak") in pairs
    assert ("CT3", "cuttag_macs3_narrow_consensus_final_peak") in pairs
    assert ("CT2", "cuttag_macs3_narrow_idr_final_peak") in pairs
    idr_rows = _build_idr_rows(_groups(), "/tmp/results", True, {})
    assert {row["output_type"] for row in idr_rows} >= {
        "idr_conservative",
        "idr_optimal",
    }
    assert all(row["experiment_id"] == "GRP2" for row in idr_rows)


def test_three_replicate_consensus_is_primary_and_two_replicate_idr_stays():
    routes = {
        "GRP2": _route("idr", 2, False),
        "GRP3": _route("consensus", 3, True),
        "ATAC3": _route("consensus", 3, True),
        "CT2": _route("idr", 2, False),
        "CT3": _route("consensus", 3, True),
    }
    rows = _repro_rows(routes)
    pairs = _pairs(rows)
    by_key = {(row["experiment_id"], row["output_type"]): row for row in rows}

    assert ("GRP2", "chipseq_macs3_narrow_consensus_peak") not in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_summary") not in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_final_peak") not in pairs
    assert ("GRP3", "chipseq_macs3_narrow_consensus_peak") in pairs
    assert ("GRP3", "chipseq_macs3_narrow_consensus_summary") in pairs
    final = by_key[("GRP3", "chipseq_macs3_narrow_consensus_final_peak")]
    assert final["method"] == "cp"
    assert final["path"].endswith(
        "GRP3.chipseq.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    assert "GRP2" not in final["path"]

    atac = by_key[("ATAC3", "atac_macs3_narrow_consensus_final_peak")]
    assert atac["method"] == "cp"
    assert atac["path"].endswith(
        "ATAC3.atac.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    assert ("ATAC3", "atac_macs3_narrow_idr_final_peak") not in pairs

    assert ("CT2", "cuttag_macs3_narrow_consensus_peak") not in pairs
    assert ("CT2", "cuttag_macs3_narrow_consensus_summary") not in pairs
    assert ("CT2", "cuttag_macs3_narrow_consensus_final_peak") not in pairs
    assert ("CT2", "cuttag_macs3_narrow_idr_final_peak") in pairs
    assert ("CT2", "cuttag_seacr_consensus_peak") in pairs
    assert ("CT2", "cuttag_seacr_consensus_final_peak") in pairs
    assert ("CT3", "cuttag_macs3_narrow_consensus_final_peak") in pairs
    assert ("CT3", "cuttag_macs3_narrow_idr_final_peak") not in pairs
    assert ("CT3", "cuttag_seacr_consensus_final_peak") in pairs

    idr_rows = _build_idr_rows(_groups(), "/tmp/results", True, routes)
    assert {row["experiment_id"] for row in idr_rows} == {"GRP2"}
    assert ("GRP3", "idr_conservative") not in {
        (row["experiment_id"], row["output_type"]) for row in idr_rows
    }


def test_two_replicate_auxiliary_consensus_stays_when_route_runs_it():
    routes = {"GRP2": _route("idr", 2, True)}
    pairs = _pairs(_repro_rows(routes))
    assert ("GRP2", "chipseq_macs3_narrow_consensus_peak") in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_summary") in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_final_peak") not in pairs
    idr_rows = _build_idr_rows(_groups(), "/tmp/results", True, routes)
    assert any(row["output_type"] == "idr_optimal" for row in idr_rows)


def test_manifest_passes_config_and_idr_flags_to_the_shared_helpers(monkeypatch, tmp_path):
    seen = {}

    def loader(path, **kwargs):
        del path
        seen["loader"] = kwargs
        return [_sample("EXP1", "chipseq", "narrow", 1)]

    def routes(samples, config, idr_flags):
        seen["samples"] = samples
        seen["config"] = config
        seen["idr_flags"] = idr_flags
        return {}

    monkeypatch.setattr(
        "encode_pipeline.manifest.make.load_and_validate_samples", loader
    )
    monkeypatch.setattr(
        "encode_pipeline.manifest.make._resolve_mixed_replicate_routes", routes
    )
    config = {
        "outdir": str(tmp_path / "results"),
        "samples": str(tmp_path / "samples.tsv"),
        "use_control": False,
        "replicate_analysis": False,
        "chipseq_idr": False,
        "multiqc": False,
        "qc": {},
    }
    build_manifest_rows(config)
    assert seen["loader"]["mixed_replicate_config"] is config
    assert seen["loader"]["chipseq_idr_enabled"] is False
    assert seen["loader"]["reproducibility_idr_atac_narrow"] is False
    assert seen["loader"]["reproducibility_idr_cuttag_narrow"] is False
    assert seen["loader"]["reproducibility_idr_chipseq_broad"] is False
    assert seen["loader"]["reproducibility_idr_cuttag_broad"] is False
    assert seen["config"] is config
    assert seen["idr_flags"]["chipseq_narrow"] is False
    assert seen["samples"][0]["experiment"] == "EXP1"


def test_missing_route_module_keeps_strict_empty_and_refuses_mixed_policy(monkeypatch):
    monkeypatch.setitem(
        __import__("sys").modules,
        "encode_pipeline.samples.mixed_replicates",
        None,
    )
    samples = [_sample("EXP1", "chipseq", "narrow", 1)]
    flags = _flags()
    assert (
        _resolve_mixed_replicate_routes(
            samples, {"reproducibility": {"enabled": True}}, flags
        )
        == {}
    )
    assert (
        _resolve_mixed_replicate_routes(
            samples,
            {
                "reproducibility": {
                    "enabled": True,
                    "mixed_replicate_policy": "strict",
                }
            },
            flags,
        )
        == {}
    )
    assert (
        _resolve_mixed_replicate_routes(
            samples,
            {
                "reproducibility": {
                    "enabled": False,
                    "mixed_replicate_policy": "per_group_fallback",
                }
            },
            flags,
        )
        == {}
    )
    with pytest.raises(ImportError, match="mixed_replicates"):
        _resolve_mixed_replicate_routes(
            samples,
            {
                "reproducibility": {
                    "enabled": True,
                    "mixed_replicate_policy": "per_group_fallback",
                }
            },
            flags,
        )
