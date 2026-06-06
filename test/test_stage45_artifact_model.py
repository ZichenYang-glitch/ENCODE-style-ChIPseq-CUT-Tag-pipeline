#!/usr/bin/env python3
"""Stage 45 artifact model tests.

Validates the Artifact dataclass, loader, and validate_artifact() helper.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from lib.artifact import (
    Artifact,
    load_artifacts,
    validate_artifact,
    VALID_SCOPES,
    VALID_LEVELS,
    VALID_ASSAY_GATES,
)

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


def main():
    global _PASS, _FAIL

    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")

    # 1. Load returns non-empty list
    artifacts = load_artifacts(inv_path)
    _check("1-load_artifacts_returns_list",
           isinstance(artifacts, list) and len(artifacts) > 0,
           f"got {type(artifacts).__name__} with {len(artifacts) if isinstance(artifacts, list) else 'N/A'} entries")

    # 2. Entry count
    _check("2-entry_count", len(artifacts) == 62,
           f"expected 62, got {len(artifacts)}")

    # 3. All entries are Artifact instances
    all_art = all(isinstance(a, Artifact) for a in artifacts)
    bad_types = [a for a in artifacts if not isinstance(a, Artifact)]
    _check("3-all_are_artifacts", all_art,
           f"{len(bad_types)} entries not Artifact instances")

    # 4. Unique IDs
    ids = [a.id for a in artifacts]
    dupes = sorted(set(i for i in ids if ids.count(i) > 1))
    _check("4-unique_ids", len(dupes) == 0,
           f"duplicates: {dupes}")

    # 5. Frozen -- mutation raises FrozenInstanceError
    try:
        artifacts[0].id = "mutated"
        _check("5-frozen_no_mutation", False,
               "expected FrozenInstanceError but no exception raised")
    except Exception as e:
        _check("5-frozen_no_mutation",
               type(e).__name__ == "FrozenInstanceError",
               f"got {type(e).__name__}: {e}")

    # 6. 13 fields via dataclasses.fields()
    from dataclasses import fields
    field_names = {f.name for f in fields(Artifact)}
    _check("6-schema_13_fields", field_names == {
        "id", "description", "scope", "level", "assay_gate",
        "path_template", "producing_rule", "tool",
        "manifest_output_type", "pipeline_done", "rule_all",
        "config_gate", "notes",
    })

    # 7. Valid scope
    bad_scope = [(a.id, a.scope) for a in artifacts
                 if a.scope not in VALID_SCOPES]
    _check("7-valid_scope", len(bad_scope) == 0,
           f"invalid scopes: {bad_scope}")

    # 8. Valid level
    bad_level = [(a.id, a.level) for a in artifacts
                 if a.level not in VALID_LEVELS]
    _check("8-valid_level", len(bad_level) == 0,
           f"invalid levels: {bad_level}")

    # 9. Valid assay_gate
    bad_gate = [(a.id, a.assay_gate) for a in artifacts
                if a.assay_gate not in VALID_ASSAY_GATES]
    _check("9-valid_assay_gate", len(bad_gate) == 0,
           f"invalid assay_gates: {bad_gate}")

    # 10. path_template starts with results/
    bad_path = [(a.id, a.path_template) for a in artifacts
                if not a.path_template.startswith("results/")]
    _check("10-path_starts_with_results", len(bad_path) == 0,
           f"bad paths: {bad_path}")

    # 11. manifest_output_type is str or None
    bad_mot = [(a.id, type(a.manifest_output_type).__name__)
               for a in artifacts
               if a.manifest_output_type is not None
               and not isinstance(a.manifest_output_type, str)]
    _check("11-manifest_output_type_str_or_none", len(bad_mot) == 0,
           f"bad types: {bad_mot}")

    # 12. pipeline_done is bool
    bad_pd = [(a.id, type(a.pipeline_done).__name__) for a in artifacts
              if not isinstance(a.pipeline_done, bool)]
    _check("12-pipeline_done_bool", len(bad_pd) == 0,
           f"not bool: {bad_pd}")

    # 13. rule_all is bool
    bad_ra = [(a.id, type(a.rule_all).__name__) for a in artifacts
              if not isinstance(a.rule_all, bool)]
    _check("13-rule_all_bool", len(bad_ra) == 0,
           f"not bool: {bad_ra}")

    # 14. Round-trip from dict
    d = {
        "id": "test_roundtrip",
        "description": "round-trip test",
        "scope": "project",
        "level": "project",
        "assay_gate": "all",
        "path_template": "results/test/test.txt",
        "producing_rule": "test_rule",
        "tool": "test_tool",
        "manifest_output_type": None,
        "pipeline_done": True,
        "rule_all": False,
        "config_gate": None,
        "notes": "",
    }
    a = Artifact(**d)
    _check("14-from_dict_roundtrip",
           a.id == "test_roundtrip" and a.pipeline_done is True
           and a.rule_all is False and a.manifest_output_type is None)

    # 15. validate_artifact rejects bad scope
    bad = dict(d)
    bad["scope"] = "invalid_scope"
    errors = validate_artifact(bad)
    _check("15-validate_rejects_bad_scope", len(errors) > 0,
           f"expected errors, got {errors}")

    # 16. validate_artifact rejects missing fields
    minimal = {"id": "minimal"}
    errors = validate_artifact(minimal)
    _check("16-validate_rejects_missing_fields", len(errors) > 0,
           f"expected errors, got {errors}")

    # 17. validate_artifact(None) returns errors, no exception
    try:
        errors = validate_artifact(None)
        _check("17-validate_none_no_crash",
               isinstance(errors, list) and len(errors) > 0,
               f"expected non-empty list, got {errors!r}")
    except Exception as e:
        _check("17-validate_none_no_crash", False,
               f"unexpected exception: {type(e).__name__}: {e}")

    # 18. path_template = None returns errors
    bad_pt = dict(d)
    bad_pt["path_template"] = None
    errors = validate_artifact(bad_pt)
    _check("18-path_template_none_rejected", len(errors) > 0,
           f"expected errors, got {errors}")

    # 19. config_gate = 7 returns errors
    bad_cg = dict(d)
    bad_cg["config_gate"] = 7
    errors = validate_artifact(bad_cg)
    _check("19-config_gate_int_rejected", len(errors) > 0,
           f"expected errors, got {errors}")

    # 20. non-string id rejected
    bad_id = dict(d)
    bad_id["id"] = 42
    errors = validate_artifact(bad_id)
    _check("20-non_string_id_rejected", len(errors) > 0,
           f"expected errors, got {errors}")

    # 21. non-string notes/tool rejected
    bad_notes = dict(d)
    bad_notes["notes"] = ["list", "not", "string"]
    errors = validate_artifact(bad_notes)
    _check("21-non_string_notes_rejected", len(errors) > 0,
           f"expected errors, got {errors}")

    bad_tool = dict(d)
    bad_tool["tool"] = 3.14
    errors = validate_artifact(bad_tool)
    _check("22-non_string_tool_rejected", len(errors) > 0,
           f"expected errors, got {errors}")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
