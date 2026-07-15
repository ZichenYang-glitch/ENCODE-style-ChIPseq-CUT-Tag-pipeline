"""Executable drift checks for the documented output contract."""

import re
from pathlib import Path

from encode_pipeline.artifacts import artifacts_by_id, load_catalog


REPO_ROOT = Path(__file__).resolve().parents[2]
INVENTORY_PATH = REPO_ROOT / "docs" / "architecture" / "artifact-inventory.yaml"
OUTPUT_CONTRACT_PATH = REPO_ROOT / "docs" / "output-contract.md"


def _parse_current_output_rows(text):
    start = text.find("## Current Output Types")
    end = text.find("## Manifest Field Schema")
    if start == -1 or end == -1 or end <= start:
        return {}, []

    rows = {}
    duplicates = []
    pattern = re.compile(r"\|\s*`([^`]+)`\s*\|[^|]*\|[^|]*\|\s*`([^`]+)`\s*\|")
    for output_type, path in pattern.findall(text[start:end]):
        if output_type in rows:
            duplicates.append(output_type)
        rows[output_type] = path
    return rows, duplicates


def _parse_current_output_producers(text):
    start = text.find("## Current Output Types")
    end = text.find("## Manifest Field Schema")
    if start == -1 or end == -1 or end <= start:
        return {}

    pattern = re.compile(
        r"\|\s*`([^`]+)`\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*`[^`]+`\s*\|"
    )
    return {
        output_type: {"tool": tool.strip(), "rule": rule.strip()}
        for output_type, tool, rule in pattern.findall(text[start:end])
    }


def _normalize_documented_path(path):
    return path.replace("<exp>", "<experiment>")


def test_output_contract_keeps_parseable_section_boundaries():
    text = OUTPUT_CONTRACT_PATH.read_text(encoding="utf-8")

    assert "## Current Output Types" in text
    assert "## Manifest Field Schema" in text
    assert text.index("## Current Output Types") < text.index(
        "## Manifest Field Schema"
    )


def test_output_contract_table_matches_artifact_inventory():
    text = OUTPUT_CONTRACT_PATH.read_text(encoding="utf-8")
    documented, duplicates = _parse_current_output_rows(text)
    inventory = artifacts_by_id(load_catalog(str(INVENTORY_PATH)))

    assert documented, "No output rows parsed from docs/output-contract.md"
    assert not duplicates, f"Duplicate documented output types: {sorted(duplicates)}"

    documented_ids = set(documented)
    inventory_ids = set(inventory)
    assert documented_ids == inventory_ids, (
        f"Only in output contract: {sorted(documented_ids - inventory_ids)}; "
        f"only in artifact inventory: {sorted(inventory_ids - documented_ids)}"
    )

    mismatches = {
        output_type: {
            "documented": _normalize_documented_path(documented[output_type]),
            "inventory": inventory[output_type].path_template,
        }
        for output_type in sorted(documented_ids)
        if _normalize_documented_path(documented[output_type])
        != inventory[output_type].path_template
    }
    assert not mismatches, f"Output path template drift: {mismatches}"


def test_output_contract_idr_producers_match_artifact_inventory():
    documented = _parse_current_output_producers(
        OUTPUT_CONTRACT_PATH.read_text(encoding="utf-8")
    )
    inventory = artifacts_by_id(load_catalog(str(INVENTORY_PATH)))
    idr_output_types = {
        "atac_macs3_narrow_idr_final_peak",
        "atac_macs3_narrow_idr_summary",
        "cuttag_macs3_narrow_idr_final_peak",
        "cuttag_macs3_narrow_idr_summary",
        "chipseq_macs3_broad_idr_final_peak",
        "chipseq_macs3_broad_idr_summary",
        "cuttag_macs3_broad_idr_final_peak",
        "cuttag_macs3_broad_idr_summary",
    }

    assert idr_output_types <= documented.keys()
    mismatches = {
        output_type: {
            "documented": documented[output_type],
            "inventory": {
                "tool": inventory[output_type].tool,
                "rule": inventory[output_type].producing_rule,
            },
        }
        for output_type in sorted(idr_output_types)
        if documented[output_type]
        != {
            "tool": inventory[output_type].tool,
            "rule": inventory[output_type].producing_rule,
        }
    }
    assert not mismatches, f"IDR producer metadata drift: {mismatches}"
