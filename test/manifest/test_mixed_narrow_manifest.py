"""Mixed reproducibility rows generated from real sample sheets."""

from pathlib import Path

import pytest

from encode_pipeline.errors import ValidationError
from encode_pipeline.manifest import build_manifest_rows, write_manifest_tsv


ASSAYS = ("chipseq", "atac", "cuttag")


def _config(tmp_path, *, counts=(2, 3), minimum=2, cuttag_pair_layout="PE"):
    samples = tmp_path / "samples.tsv"
    lines = [
        "sample\tlayout\tfastq_1\tfastq_2\tassay\tpeak_mode\tgenome\tbowtie2_index"
        "\ttarget\texperiment\tbiological_replicate"
    ]
    for assay in ASSAYS:
        for count in counts:
            exp = f"{assay}{count}"
            layout = cuttag_pair_layout if exp == "cuttag2" else "PE"
            for rep in range(1, count + 1):
                sample = f"{exp}_r{rep}"
                r2 = f"{sample}_R2.fq" if layout == "PE" else ""
                lines.append(
                    f"{sample}\t{layout}\t{sample}_R1.fq\t{r2}\t{assay}\tnarrow"
                    f"\ths\t/idx\tCTCF\t{exp}\t{rep}"
                )
    samples.write_text("\n".join(lines) + "\n")
    return {
        "outdir": str(tmp_path / "results"),
        "samples": str(samples),
        "blacklist_mode": "off",
        "chipseq_idr": True,
        "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
        "reproducibility": {
            "enabled": True,
            "mixed_replicate_policy": "per_group_fallback",
            "consensus": {"enabled": True, "min_replicates": minimum},
            "idr": {"atac_narrow": True, "cuttag_narrow": True},
        },
    }


@pytest.mark.parametrize("enable_idr", [True, False])
def test_strict_and_omitted_policy_write_identical_manifest_bytes(tmp_path, enable_idr):
    config = _config(tmp_path, counts=(2,) if enable_idr else (2, 3))
    if not enable_idr:
        config["chipseq_idr"] = False
        config["reproducibility"]["idr"] = {}
    del config["reproducibility"]["mixed_replicate_policy"]
    rows, _, _ = build_manifest_rows(config)
    legacy = tmp_path / "omitted.tsv"
    write_manifest_tsv(rows, legacy)
    config["reproducibility"]["mixed_replicate_policy"] = "strict"
    rows, _, _ = build_manifest_rows(config)
    explicit = tmp_path / "strict.tsv"
    write_manifest_tsv(rows, explicit)
    assert explicit.read_bytes() == legacy.read_bytes()
    assert explicit.read_bytes().splitlines()[0].split(b"\t") == [
        b"sample_id",
        b"experiment_id",
        b"assay",
        b"target",
        b"genome",
        b"output_type",
        b"method",
        b"path",
        b"status",
        b"qc_flag",
    ]
    types = {row["output_type"] for row in rows}
    assert "cuttag_seacr_consensus_final_peak" in types
    assert ("idr_optimal" in types) == enable_idr
    assert ("idr_conservative" in types) == enable_idr
    assert ("cuttag_macs3_narrow_consensus_final_peak" in types) == (not enable_idr)
    assert "chipseq_macs3_narrow_consensus_final_peak" not in types
    assert "atac_macs3_narrow_consensus_final_peak" not in types


@pytest.mark.parametrize("minimum", [2, 3])
def test_mixed_manifest_selects_finals_by_experiment_and_caller(tmp_path, minimum):
    config = _config(
        tmp_path, minimum=minimum, cuttag_pair_layout="SE" if minimum == 3 else "PE"
    )
    experiments = Path(config["outdir"]) / "experiments"
    # Exercise the same existing auxiliary files under each support threshold.
    auxiliary = experiments / "chipseq2/06_reproducibility/consensus"
    auxiliary.mkdir(parents=True)
    for suffix in ("narrowPeak", "summary.tsv"):
        (auxiliary / f"chipseq2.chipseq.macs3.narrow.consensus.{suffix}").touch()

    rows, _, _ = build_manifest_rows(config)
    records = {(row["experiment_id"], row["output_type"]): row for row in rows}
    present = {Path(row["path"]) for row in rows if row["status"] == "present"}
    assert present == (set(auxiliary.iterdir()) if minimum == 2 else set())
    for assay in ASSAYS:
        pair, triple = f"{assay}2", f"{assay}3"
        prefix = f"{assay}_macs3_narrow"
        for kind in ("peak", "summary"):
            assert ((pair, f"{prefix}_consensus_{kind}") in records) == (minimum == 2)
            assert (triple, f"{prefix}_consensus_{kind}") in records

        idr = (
            ("idr_conservative", "idr_optimal")
            if assay == "chipseq"
            else (f"{prefix}_idr_final_peak",)
        )
        for output_type in idr:
            assert records[(pair, output_type)]["method"] == "idr"
            assert (triple, output_type) not in records

        final_type = f"{prefix}_consensus_final_peak"
        assert (pair, final_type) not in records
        final = records[(triple, final_type)]
        assert final["method"] == "cp"
        assert Path(final["path"]) == (
            experiments
            / triple
            / "06_reproducibility/final"
            / f"{triple}.{assay}.macs3.narrow.replicate_validated.consensus.narrowPeak"
        )

    assert ("cuttag3", "cuttag_seacr_consensus_final_peak") in records
    assert (("cuttag2", "cuttag_seacr_consensus_final_peak") in records) == (
        minimum == 2
    )


def test_three_replicate_atac_needs_its_own_idr_opt_in_for_a_consensus_final(tmp_path):
    config = _config(tmp_path)
    config["reproducibility"]["idr"]["atac_narrow"] = False
    rows, _, _ = build_manifest_rows(config)
    types = {row["output_type"] for row in rows if row["experiment_id"] == "atac3"}
    assert "atac_macs3_narrow_consensus_peak" in types
    assert "atac_macs3_narrow_consensus_final_peak" not in types


def test_three_of_three_rejects_an_under_supported_seacr_final(tmp_path):
    with pytest.raises(ValidationError, match="SEACR"):
        build_manifest_rows(_config(tmp_path, minimum=3))
