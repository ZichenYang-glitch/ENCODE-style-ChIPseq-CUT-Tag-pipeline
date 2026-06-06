#!/usr/bin/env python3
"""Stage 43 artifact inventory consistency tests.

Checks:
  1. YAML parses successfully (fails if PyYAML unavailable)
  2. All artifact ids are unique
  3. No path_template contains absolute workstation paths
  4. No entry mentions future-feature keywords as implemented
  5. Every manifest output_type from make_manifest.py has a matching
     inventory entry (using ast to inspect _add_row calls)
"""

import ast
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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

    # 1. Load YAML (fail if unavailable)
    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")
    import yaml

    with open(inv_path) as fh:
        data = yaml.safe_load(fh)

    _check("1-yaml_parses", data is not None and "artifacts" in data)
    if data is None or "artifacts" not in data:
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1 if _FAIL else 0)

    artifacts = data["artifacts"]
    _check("1b-has_entries", len(artifacts) >= 50,
           f"expected >=50, got {len(artifacts)}")

    # 2. Unique IDs
    seen = set()
    dupes = []
    for a in artifacts:
        aid = a["id"]
        if aid in seen:
            dupes.append(aid)
        seen.add(aid)
    _check("2-unique_ids", len(dupes) == 0,
           f"duplicates: {dupes}")

    # 3. No absolute workstation paths
    BAD_PREFIXES = ("/home/", "/data/", "/mnt/")
    bad_paths = []
    for a in artifacts:
        pt = a.get("path_template", "")
        if any(pt.startswith(p) for p in BAD_PREFIXES):
            bad_paths.append(f"{a['id']}: {pt}")
    _check("3-no_workstation_paths", len(bad_paths) == 0,
           f"bad paths: {bad_paths}")

    # 4. No artifact-model keywords as implemented features
    BAD_KEYWORDS = ("Artifact dataclass", "artifact_path()",
                    "paths.smk", "AssayPolicy YAML",
                    "global target resolver")
    bad_notes = []
    for a in artifacts:
        notes = a.get("notes", "") or ""
        for kw in BAD_KEYWORDS:
            if kw.lower() in notes.lower():
                bad_notes.append(f"{a['id']}: mentions '{kw}'")
    _check("4-no_future_keywords_in_notes", len(bad_notes) == 0,
           f"entries: {bad_notes}")

    # 5. Every manifest output_type has an inventory entry
    manifest_path = os.path.join(REPO_ROOT, "scripts", "make_manifest.py")
    manifest_types = _extract_manifest_output_types(manifest_path)

    # Also add known dynamic f-string types that AST correctly normalizes
    # to <N> placeholders.  The AST visitor already does this.

    inventory_manifest_types = set()
    for a in artifacts:
        mt = a.get("manifest_output_type")
        if mt is not None:
            inventory_manifest_types.add(mt)

    missing = manifest_types - inventory_manifest_types
    _check("5-manifest_coverage", len(missing) == 0,
           f"manifest types not in inventory: {sorted(missing)}")

    # 6. Every non-null inventory manifest_output_type is emitted by manifest
    non_null_mt = {a["manifest_output_type"] for a in artifacts
                   if a.get("manifest_output_type") is not None}
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
    inventory_ids_all = {a["id"] for a in artifacts}
    oc_missing = oc_types - inventory_ids_all
    _check("7-oc_types_in_inventory", len(oc_missing) == 0,
           f"OC types not in inventory: {sorted(oc_missing)}")
    inv_extra = inventory_ids_all - oc_types
    _check("8-inventory_ids_in_oc", len(inv_extra) == 0,
           f"inventory ids not in OC: {sorted(inv_extra)}")

    # 9. Schema completeness — every entry has exactly the 13 required fields,
    #    and pipeline_done / rule_all are booleans.
    REQUIRED_FIELDS = {
        "id", "description", "scope", "level", "assay_gate",
        "path_template", "producing_rule", "tool", "manifest_output_type",
        "pipeline_done", "rule_all", "config_gate", "notes",
    }
    schema_errors = []
    for a in artifacts:
        keys = set(a.keys())
        missing_fields = REQUIRED_FIELDS - keys
        extra_fields = keys - REQUIRED_FIELDS
        if missing_fields:
            schema_errors.append(f"{a['id']}: missing {sorted(missing_fields)}")
        if extra_fields:
            schema_errors.append(f"{a['id']}: extra {sorted(extra_fields)}")
        if not isinstance(a.get("pipeline_done"), bool):
            schema_errors.append(
                f"{a['id']}: pipeline_done not bool ({type(a.get('pipeline_done')).__name__})")
        if not isinstance(a.get("rule_all"), bool):
            schema_errors.append(
                f"{a['id']}: rule_all not bool ({type(a.get('rule_all')).__name__})")
    _check("9-schema_completeness", len(schema_errors) == 0,
           f"schema errors: {schema_errors[:5]}{'...' if len(schema_errors) > 5 else ''}")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
