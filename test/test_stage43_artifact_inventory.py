#!/usr/bin/env python3
"""Stage 43 artifact inventory consistency tests.

Uses Stage 45's Artifact model (load_artifacts, Artifact) as the
canonical inventory read path.

Checks:
  1. load_artifacts succeeds and returns entries
  2. All artifact ids are unique
  3. No path_template contains absolute workstation paths
  4. No entry mentions future-feature keywords as implemented
  5. Every manifest output_type from make_manifest.py has a matching
     inventory entry (using ast to inspect _add_row calls)
  6. No extra manifest types in inventory not emitted by manifest
  7. Output-contract types match inventory ids (both directions)
  8. Schema completeness owned by Artifact dataclass
"""

import ast
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from dataclasses import fields
from lib.artifact import Artifact, load_artifacts, artifacts_by_manifest_output_type

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


def _extract_manifest_output_types(manifest_path):
    """Parse make_manifest.py AST to extract output_type strings from
    _add_row(...) calls. Returns a set of normalized output_type values.

    For f-string JoinedStr values like f"biorep{br}_final_bam", normalizes
    to biorep<N>_final_bam.
    """
    with open(manifest_path) as fh:
        tree = ast.parse(fh.read(), filename=manifest_path)

    types = set()

    class AddRowVisitor(ast.NodeVisitor):
        def visit_Call(self, node):
            # Match _add_row(rows, ...) calls
            if (isinstance(node.func, ast.Name)
                    and node.func.id == "_add_row"):
                # The output_type is always the 2nd positional argument
                # after the rows list: _add_row(rows, sid, exp, assay, ...)
                # Actually looking at the code pattern, output_type is
                # positional arg index 4 (0-based): rows, sid, exp, assay,
                # target, genome, output_type, method, path, ...
                # Let's find the string literal or f-string
                if len(node.args) >= 7:
                    arg = node.args[6]  # output_type is 7th positional arg (0-indexed 6)
                    val = _extract_string(arg)
                    if val:
                        types.add(val)
            self.generic_visit(node)

    def _extract_string(node):
        """Extract a string value from an AST node, normalizing f-strings."""
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value
        if isinstance(node, ast.JoinedStr):
            # Reconstruct f-string pattern with <N> placeholder
            parts = []
            for v in node.values:
                if isinstance(v, ast.Constant) and isinstance(v.value, str):
                    parts.append(v.value)
                elif isinstance(v, ast.FormattedValue):
                    parts.append("<N>")
            return "".join(parts)
        return None

    AddRowVisitor().visit(tree)

    # Also catch _add_row calls inside the manifest that use the
    # _add_row(rows, sid, exp, assay, target, genome, output_type, ...)
    # pattern. The visitor above catches all.
    return types


def main():
    global _PASS, _FAIL

    # 1. Load via Stage 45 Artifact model
    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")

    try:
        artifacts = load_artifacts(inv_path)
    except Exception as e:
        _check("1-load_artifacts_succeeds", False,
               f"load_artifacts raised {type(e).__name__}: {e}")
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1)

    _check("1-load_artifacts_succeeds", True)
    _check("1b-has_entries", len(artifacts) >= 50,
           f"expected >=50, got {len(artifacts)}")

    # 2. Unique IDs
    seen = set()
    dupes = []
    for a in artifacts:
        if a.id in seen:
            dupes.append(a.id)
        seen.add(a.id)
    _check("2-unique_ids", len(dupes) == 0,
           f"duplicates: {dupes}")

    # 3. No absolute workstation paths
    BAD_PREFIXES = ("/home/", "/data/", "/mnt/")
    bad_paths = []
    for a in artifacts:
        if any(a.path_template.startswith(p) for p in BAD_PREFIXES):
            bad_paths.append(f"{a.id}: {a.path_template}")
    _check("3-no_workstation_paths", len(bad_paths) == 0,
           f"bad paths: {bad_paths}")

    # 4. No artifact-model keywords as implemented features
    BAD_KEYWORDS = ("Artifact dataclass", "artifact_path()",
                    "paths.smk", "AssayPolicy YAML",
                    "global target resolver")
    bad_notes = []
    for a in artifacts:
        for kw in BAD_KEYWORDS:
            if kw.lower() in a.notes.lower():
                bad_notes.append(f"{a.id}: mentions '{kw}'")
    _check("4-no_future_keywords_in_notes", len(bad_notes) == 0,
           f"entries: {bad_notes}")

    # 5. Every manifest output_type has an inventory entry
    manifest_path = os.path.join(REPO_ROOT, "scripts", "make_manifest.py")
    manifest_types = _extract_manifest_output_types(manifest_path)

    inventory_manifest_types = {a.manifest_output_type for a in artifacts
                                if a.manifest_output_type is not None}

    missing = manifest_types - inventory_manifest_types
    _check("5-manifest_coverage", len(missing) == 0,
           f"manifest types not in inventory: {sorted(missing)}")

    # 6. Every non-null inventory manifest_output_type is emitted by manifest
    by_mot = artifacts_by_manifest_output_type(artifacts)
    non_null_mt = set(by_mot.keys())
    extra_mt = non_null_mt - manifest_types
    _check("6-no_extra_manifest_types", len(extra_mt) == 0,
           f"inventory types not in manifest: {sorted(extra_mt)}")

    # 7. Output-contract cross-reference: inventory ids match exactly.
    #    Extract output_type entries from between "## Current Output Types"
    #    and "## Gating Conditions" — the only sections with output tables.
    oc_path = os.path.join(REPO_ROOT, "docs", "output-contract.md")
    with open(oc_path) as fh:
        oc_text = fh.read()
    start = oc_text.find("## Current Output Types")
    end = oc_text.find("## Manifest Field Schema")
    if start == -1 or end == -1:
        oc_types = set()
    else:
        section = oc_text[start:end]
        oc_types = set(re.findall(r'^\| `([^`]+)` \|', section, re.MULTILINE))
    inventory_ids_all = {a.id for a in artifacts}
    oc_missing = oc_types - inventory_ids_all
    _check("7-oc_types_in_inventory", len(oc_missing) == 0,
           f"OC types not in inventory: {sorted(oc_missing)}")
    inv_extra = inventory_ids_all - oc_types
    _check("8-inventory_ids_in_oc", len(inv_extra) == 0,
           f"inventory ids not in OC: {sorted(inv_extra)}")

    # 9. Schema completeness — owned by Artifact dataclass
    artifact_field_names = {f.name for f in fields(Artifact)}
    schema_errors = []
    for a in artifacts:
        if not isinstance(a, Artifact):
            schema_errors.append(f"not an Artifact: {type(a).__name__}")
            continue
        actual_fields = set(vars(a).keys())
        missing = artifact_field_names - actual_fields
        extra = actual_fields - artifact_field_names
        if missing:
            schema_errors.append(f"{a.id}: missing fields {sorted(missing)}")
        if extra:
            schema_errors.append(f"{a.id}: extra fields {sorted(extra)}")
    _check("9-schema_completeness",
           len(artifact_field_names) == 13 and len(schema_errors) == 0,
           f"field count={len(artifact_field_names)}; errors: {schema_errors[:5]}{'...' if len(schema_errors) > 5 else ''}")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
