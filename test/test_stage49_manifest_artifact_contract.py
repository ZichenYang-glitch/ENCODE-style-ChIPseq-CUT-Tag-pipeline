#!/usr/bin/env python3
"""Stage 49: Manifest-artifact contract tests.

Verifies that manifest output_type vocabulary matches the artifact
inventory and that tracked Artifacts satisfy basic invariants.
"""

import ast
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from lib.artifact import load_artifacts, artifacts_by_manifest_output_type

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
    """AST-extract output_type strings from _add_row(...) calls."""
    with open(manifest_path) as fh:
        tree = ast.parse(fh.read(), filename=manifest_path)

    types = set()

    class AddRowVisitor(ast.NodeVisitor):
        def visit_Call(self, node):
            if (isinstance(node.func, ast.Name)
                    and node.func.id == "_add_row"):
                if len(node.args) >= 7:
                    arg = node.args[6]
                    val = _extract_string(arg)
                    if val:
                        types.add(val)
            self.generic_visit(node)

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

    # 2. AST extracts manifest output_type set
    manifest_path = os.path.join(REPO_ROOT, "scripts", "make_manifest.py")
    manifest_types = _extract_manifest_output_types(manifest_path)
    _check("2-manifest_types_from_ast",
           isinstance(manifest_types, set) and len(manifest_types) > 0,
           f"got {len(manifest_types)} types")

    # 3. no duplicate manifest_output_type in inventory
    try:
        by_mot = artifacts_by_manifest_output_type(artifacts)
    except ValueError as e:
        _check("3-no_duplicate_manifest_output_type", False,
               f"ValueError: {e}")
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1)
    _check("3-no_duplicate_manifest_output_type", True)

    # 4. bidirectional manifest ↔ inventory
    inventory_mot = set(by_mot.keys())
    missing_from_inv = manifest_types - inventory_mot
    extra_in_inv = inventory_mot - manifest_types
    _check("4-manifest_vs_inventory_bidirectional",
           len(missing_from_inv) == 0 and len(extra_in_inv) == 0,
           f"missing from inventory: {sorted(missing_from_inv)}; "
           f"extra in inventory: {sorted(extra_in_inv)}")

    # 5. tracked mapping values match keys
    _check("5-tracked_mapping_values_match_keys",
           all(a.manifest_output_type == key for key, a in by_mot.items()),
           "value.manifest_output_type != key")

    # 6. producing_rule non-empty for tracked artifacts
    bad_rules = [(k, a.producing_rule) for k, a in by_mot.items()
                 if not a.producing_rule]
    _check("6-producing_rule_non_empty", len(bad_rules) == 0,
           f"empty producing_rule: {bad_rules}")

    # 7. path_template starts with results/
    bad_paths = [(k, a.path_template) for k, a in by_mot.items()
                 if not a.path_template.startswith("results/")]
    _check("7-path_starts_results", len(bad_paths) == 0,
           f"bad paths: {bad_paths}")

    # 8. tracked artifacts non-empty
    _check("8-tracked_artifacts_non_empty", len(by_mot) > 0,
           f"by_mot has {len(by_mot)} entries")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
