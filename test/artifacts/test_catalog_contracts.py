"""Contract tests for the artifact catalog."""

import ast
import csv
import io
import os
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import scripts.make_manifest as make_manifest

from encode_pipeline.artifacts import (
    Artifact,
    artifacts_by_id,
    artifacts_by_manifest_output_type,
    load_catalog,
    validate_artifact,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
INVENTORY_PATH = REPO_ROOT / "docs" / "architecture" / "artifact-inventory.yaml"
MANIFEST_SCRIPT = REPO_ROOT / "scripts" / "make_manifest.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_manifest_output_types(manifest_path):
    """AST-extract output_type strings from _add_row(...) calls."""
    with open(manifest_path) as fh:
        tree = ast.parse(fh.read(), filename=manifest_path)

    types = set()

    class AddRowVisitor(ast.NodeVisitor):
        def visit_Call(self, node):
            if (isinstance(node.func, ast.Name)
                    and node.func.id == "_add_row"):
                if len(node.args) >= 7:
                    val = self._extract_string(node.args[6])
                    if val:
                        types.add(val)
            self.generic_visit(node)

        @staticmethod
        def _extract_string(node):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                return node.value
            if isinstance(node, ast.JoinedStr):
                parts = []
                for v in node.values:
                    if isinstance(v, ast.Constant) and isinstance(v.value, str):
                        parts.append(v.value)
                    elif isinstance(v, ast.FormattedValue):
                        parts.append("<N>")
                return "".join(parts)
            return None

    AddRowVisitor().visit(tree)
    return types


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


def _run_manifest(config, samples_tsv, tmp_config, chrom_sizes=False):
    """Generate a manifest and return the set of output_type values."""
    cfg = _base_config()
    cfg.update(config)

    workdir, config_path, samples_path = tmp_config(
        config=cfg,
        samples=samples_tsv,
    )
    cfg["samples"] = str(samples_path)

    if chrom_sizes:
        cs_path = workdir / "chr.sizes"
        cs_path.write_text("", encoding="utf-8")
        cfg["genome_resources"]["hs"]["chrom_sizes"] = str(cs_path)
        with open(config_path, "w", encoding="utf-8") as fh:
            _write_yaml(fh, cfg)

    out = workdir / "manifest.tsv"

    old_argv = sys.argv
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    stdout_capture = io.StringIO()
    stderr_capture = io.StringIO()
    sys.stdout = stdout_capture
    sys.stderr = stderr_capture
    sys.argv = ["make_manifest.py", "--config", str(config_path), "--output", str(out)]
    try:
        make_manifest.main()
        returncode = 0
    except SystemExit as exc:
        returncode = exc.code if exc.code is not None else 0
    finally:
        sys.argv = old_argv
        sys.stdout = old_stdout
        sys.stderr = old_stderr

    assert returncode == 0, stderr_capture.getvalue()[-500:]

    rows = []
    if out.exists():
        with open(out, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["output_type"] for r in rows}


def _write_yaml(fh, data, indent=0):
    """Write a nested YAML mapping to a file handle."""
    prefix = "  " * indent
    for key, value in data.items():
        if isinstance(value, dict):
            fh.write(f"{prefix}{key}:\n")
            _write_yaml(fh, value, indent + 1)
        elif isinstance(value, bool):
            fh.write(f"{prefix}{key}: {str(value).lower()}\n")
        elif isinstance(value, str):
            if value:
                fh.write(f'{prefix}{key}: "{value}"\n')
            else:
                fh.write(f'{prefix}{key}: ""\n')
        elif isinstance(value, list):
            fh.write(f"{prefix}{key}:\n")
            for item in value:
                if isinstance(item, dict):
                    fh.write(f"{prefix}  -\n")
                    _write_yaml(fh, item, indent + 2)
                else:
                    fh.write(f"{prefix}  - {item}\n")
        else:
            fh.write(f"{prefix}{key}: {value}\n")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_inventory_entries_are_valid_artifacts():
    catalog = load_catalog(str(INVENTORY_PATH))
    assert isinstance(catalog, list)
    assert len(catalog) > 0
    assert all(isinstance(a, Artifact) for a in catalog)

    by_id = artifacts_by_id(catalog)
    assert len(by_id) == len(catalog)

    for entry in load_catalog(str(INVENTORY_PATH)):
        # validate_artifact operates on raw dicts, so reload as data and check
        pass

    # Re-load raw YAML to exercise validate_artifact on every raw entry
    import yaml
    with open(INVENTORY_PATH) as fh:
        data = yaml.safe_load(fh)
    for entry in data["artifacts"]:
        errors = validate_artifact(entry)
        assert errors == [], f"id={entry.get('id', '?')}: {errors}"


def test_manifest_output_types_exist_in_catalog():
    catalog = load_catalog(str(INVENTORY_PATH))
    by_mot = artifacts_by_manifest_output_type(catalog)
    catalog_mots = set(by_mot.keys())

    manifest_types = _extract_manifest_output_types(MANIFEST_SCRIPT)
    assert manifest_types, "No output_type strings extracted from make_manifest.py"

    missing = manifest_types - catalog_mots
    assert not missing, f"Manifest output_types missing from catalog: {sorted(missing)}"


def test_catalog_output_types_appear_in_generated_manifests(tmp_config):
    catalog = load_catalog(str(INVENTORY_PATH))
    catalog_mots = {a.manifest_output_type for a in catalog if a.manifest_output_type}

    columns_full = [
        "sample", "fastq_1", "fastq_2", "layout", "assay",
        "target", "peak_mode", "genome", "bowtie2_index",
        "experiment", "biological_replicate",
    ]

    all_types = set()

    # 1. chipseq narrow default (single sample)
    all_types |= _run_manifest(
        {},
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "chipseq",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n",
        tmp_config,
        chrom_sizes=True,
    )

    # 2. chipseq narrow IDR (2 bioreps)
    all_types |= _run_manifest(
        {
            "stage5": True,
            "idr": {"seed": 42, "threshold": 0.05, "rank": "p.value"},
            "qc": {"signal_tracks": True, "summary": True},
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "chipseq",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "chipseq",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
        chrom_sizes=True,
    )

    # 3. chipseq narrow consensus (2 bioreps, no IDR)
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "chipseq",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "chipseq",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 4. chipseq broad consensus (2 bioreps, no broad IDR)
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"chipseq_broad_experimental": False},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "chipseq",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "chipseq",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 4. chipseq broad IDR experimental
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"chipseq_broad_experimental": True},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "chipseq",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "chipseq",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 5. cuttag narrow consensus (2 bioreps, no narrow IDR)
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"cuttag_narrow": False},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "cuttag",
                     "H3K4me3", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "cuttag",
                     "H3K4me3", "narrow", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 6. cuttag narrow IDR
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"cuttag_narrow": True},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "cuttag",
                     "H3K4me3", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "cuttag",
                     "H3K4me3", "narrow", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 7. cuttag broad consensus (2 bioreps, no broad IDR)
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"cuttag_broad_experimental": False},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "cuttag",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "cuttag",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 8. cuttag broad IDR experimental
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"cuttag_broad_experimental": True},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "cuttag",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "cuttag",
                     "H3K27me3", "broad", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 9. cuttag SEACR consensus (2 bioreps)
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
            },
            "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "cuttag",
                     "H3K4me3", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "cuttag",
                     "H3K4me3", "narrow", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 10. atac narrow consensus + IDR (2 bioreps)
    all_types |= _run_manifest(
        {
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {"atac_narrow": True},
            },
        },
        "\t".join(columns_full) + "\n"
        + "\t".join(["S1", "R1.fq", "R2.fq", "PE", "atac",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["S2", "R3.fq", "R4.fq", "PE", "atac",
                     "CTCF", "narrow", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    # 11. MNase (2 bioreps)
    all_types |= _run_manifest(
        {},
        "\t".join(columns_full) + "\n"
        + "\t".join(["M1", "R1.fq", "R2.fq", "PE", "mnase",
                     "H3", "nucleosome", "hs", "/idx", "EXP1", "1"]) + "\n"
        + "\t".join(["M2", "R3.fq", "R4.fq", "PE", "mnase",
                     "H3", "nucleosome", "hs", "/idx", "EXP1", "2"]) + "\n",
        tmp_config,
    )

    missing = []
    for mot in sorted(catalog_mots):
        if "<N>" in mot:
            pattern = re.compile(re.escape(mot).replace(r"<N>", r"\d+"))
            if not any(pattern.fullmatch(t) for t in all_types):
                missing.append(mot)
        elif mot not in all_types:
            missing.append(mot)
    assert not missing, f"Catalog output_types not emitted by any manifest: {missing}"
