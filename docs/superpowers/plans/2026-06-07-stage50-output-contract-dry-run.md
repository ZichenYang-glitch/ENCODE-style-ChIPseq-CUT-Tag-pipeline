# Stage 50: Output Contract Generation Dry-Run — Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans inline.

**Goal:** Create `test/test_stage50_output_contract_dry_run.py` with 8 checks verifying output-contract.md ↔ artifact inventory bidirectional equivalence, including path_template matching.

---

### Task 1: Create `test/test_stage50_output_contract_dry_run.py`

```python
#!/usr/bin/env python3
"""Stage 50: Output-contract generation dry-run.

Verifies bidirectional equivalence between artifact inventory and
docs/output-contract.md output-type tables — including path_template
matching.  Dry-run only — does not write docs/output-contract.md.
"""

import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from lib.artifact import load_artifacts

_PASS = 0
_FAIL = 0


def _check(name, condition, detail=""):
    global _PASS, _FAIL
    if condition:
        print(f"PASS: {name}")
        _PASS += 1
    else:
        msg = f"FAIL: {name}"
        if detail:
            msg += f"  -- {detail}"
        print(msg)
        _FAIL += 1


def _normalize_oc_path(path):
    """Normalize OC path placeholders to match inventory conventions.

    output-contract.md uses <exp> for experiment-level paths while
    artifact-inventory.yaml uses <experiment>.  This is the only
    documented mismatch pattern (21 entries).  All sample-level paths
    match exactly.
    """
    return path.replace("<exp>", "<experiment>")


def _parse_oc_rows(oc_text):
    """Parse output-contract.md 'Current Output Types' section.

    Returns (rows, duplicates) where rows is dict[output_type, path]
    and duplicates is list of output_types that appeared more than once.
    """
    start = oc_text.find("## Current Output Types")
    end = oc_text.find("## Manifest Field Schema")
    if start == -1 or end == -1:
        return {}, []

    section = oc_text[start:end]
    rows = {}
    duplicates = []
    for m in re.finditer(
        r'\|\s*`([^`]+)`\s*\|[^|]*\|[^|]*\|\s*`([^`]+)`\s*\|',
        section
    ):
        oc_type = m.group(1)
        oc_path = m.group(2)
        if oc_type in rows:
            duplicates.append(oc_type)
        rows[oc_type] = oc_path
    return rows, duplicates


def main():
    global _PASS, _FAIL

    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")

    # 1. load_artifacts succeeds
    try:
        artifacts = load_artifacts(inv_path)
    except Exception as e:
        _check("1-load_artifacts_succeeds", False,
               f"load_artifacts raised {type(e).__name__}: {e}")
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1)
    _check("1-load_artifacts_succeeds", True)

    # Build inventory lookup
    inv_by_id = {a.id: a for a in artifacts}

    # 2. output-contract.md has section markers
    oc_path = os.path.join(REPO_ROOT, "docs", "output-contract.md")
    with open(oc_path) as fh:
        oc_text = fh.read()
    has_current = "## Current Output Types" in oc_text
    has_manifest = "## Manifest Field Schema" in oc_text
    _check("2-oc_section_markers", has_current and has_manifest,
           f"Current Output Types: {has_current}, "
           f"Manifest Field Schema: {has_manifest}")

    # 3. parse OC rows
    oc_rows, duplicate_oc_types = _parse_oc_rows(oc_text)
    _check("3-oc_rows_non_empty", len(oc_rows) > 0,
           f"parsed {len(oc_rows)} OC rows")

    # 4. no duplicate OC types (detected during parse)
    _check("4-no_duplicate_oc_types", len(duplicate_oc_types) == 0,
           f"duplicate OC types: {sorted(duplicate_oc_types)}")

    # 5. OC types in inventory
    oc_missing = set(oc_rows) - set(inv_by_id)
    _check("5-oc_types_in_inventory", len(oc_missing) == 0,
           f"OC types not in inventory: {sorted(oc_missing)}")

    # 6. Inventory IDs in OC
    inv_extra = set(inv_by_id) - set(oc_rows)
    _check("6-inventory_ids_in_oc", len(inv_extra) == 0,
           f"inventory IDs not in OC: {sorted(inv_extra)}")

    # 7. path_templates match (with normalization)
    path_mismatches = []
    for oc_type in sorted(oc_rows):
        if oc_type not in inv_by_id:
            continue
        normalized_oc = _normalize_oc_path(oc_rows[oc_type])
        inv_path = inv_by_id[oc_type].path_template
        if normalized_oc != inv_path:
            path_mismatches.append(
                f"{oc_type}: OC='{normalized_oc}' vs INV='{inv_path}'"
            )
    _check("7-path_templates_match", len(path_mismatches) == 0,
           f"mismatches: {path_mismatches}")

    # 8. dry-run output type set matches
    inv_ids = set(inv_by_id.keys())
    oc_set = set(oc_rows.keys())
    _check("8-dry_run_output_type_set_matches", inv_ids == oc_set,
           f"inventory only: {sorted(inv_ids - oc_set)}; "
           f"OC only: {sorted(oc_set - inv_ids)}")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
```

---

### Task 2: Update artifact roadmap

Add Stage 50 section. Target helpers → Stage 51+.

---

### Task 3: Run full test suite

```bash
python3 test/test_stage50_output_contract_dry_run.py        # 8/8
python3 test/test_stage49_manifest_artifact_contract.py     # 8/8
python3 test/test_stage43_artifact_inventory.py             # 10/10
python3 test/test_stage45_artifact_model.py                 # 29/29
python3 test/test_stage47_mnase_path_contract.py            # 17/17
python3 test/test_stage28_release_readiness.py              # 11/11
python3 test/test_no_hardcoded_paths.py                     # PASS
git diff --check
```

### Task 4: Commit

```bash
git add test/test_stage50_output_contract_dry_run.py \
  docs/architecture/artifact-roadmap.md \
  docs/superpowers/specs/2026-06-07-stage50-output-contract-dry-run-design.md \
  docs/superpowers/plans/2026-06-07-stage50-output-contract-dry-run.md

git commit -m "$(cat <<'EOF'
test: add output contract dry-run checks (Stage 50)

Dry-run verification of bidirectional equivalence between
docs/output-contract.md and artifact inventory, including
path_template matching with <exp>/<experiment> placeholder normalization.
No docs generation adoption — read-only contract validation.
EOF
)"
```
