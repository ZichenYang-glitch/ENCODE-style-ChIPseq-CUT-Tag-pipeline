# Config Validator Decomposition Checkpoint 2

> State as of PR71 (`05bd74a refactor(samples): split validate_replicate_groups into helpers (#71)`).

This document audits `src/encode_pipeline/config/validator.py` after the
sample-sheet extraction phase. PR68–PR71 pinned and then extracted sample
loading, per-row validation, and replicate-group validation into
`encode_pipeline.samples.*`. This checkpoint records the remaining surface in
`validator.py`, the new `samples.*` modules, the current import shape, stale
aliases, and a risk-ranked roadmap for the next PRs. No production code changes
are included; this is architecture audit and planning only.

---

## 1. Size snapshot

| File | Lines | Notes |
|---|---|---|
| `src/encode_pipeline/config/validator.py` (current) | **662** | Down from ~2,070 at Phase 0 and 1,217 at PR66. |
| `src/encode_pipeline/config/validator.py` (PR66) | 1,217 | After config-block extraction, before sample-sheet extraction. |
| `src/encode_pipeline/config/validator.py` (Phase 0) | ~2,070 | Monolithic validator. |

Net reduction since Phase 0: **~1,408 lines** (~68%).

---

## 2. Extracted modules inventory

### Config-block modules (from PR66 and earlier)

| Module | Lines | Owns |
|---|---|---|
| `encode_pipeline.config.defaults` | 114 | Pure constants, allowed-value sets, static default ranges. |
| `encode_pipeline.config.coercion` | 33 | Stateless type coercion: `coerce_int`, `coerce_bool`. |
| `encode_pipeline.config.genome` | 120 | Effective genome size validation, genome resource mapping, Picard/TSS resource gates. |
| `encode_pipeline.config.tools` | 311 | Tool parameter validation for `tool_parameters` block. |
| `encode_pipeline.config.reproducibility` | 217 | Reproducibility/IDR/consensus policy validation and warnings. |
| `encode_pipeline.config.qc` | 45 | QC config block normalization. |
| `encode_pipeline.config.cuttag` | 115 | CUT&Tag config validation. |
| `encode_pipeline.config.mnase` | 173 | MNase config validation. |

### Sample-sheet modules (new since PR66 checkpoint)

| Module | Lines | Owns |
|---|---|---|
| `encode_pipeline.errors` | ~5 | Canonical home for `ValidationError`. |
| `encode_pipeline.config.errors` | ~3 | Backward-compatible re-export of `ValidationError`. |
| `encode_pipeline.samples.load` | ~14 | Public compatibility wrapper/re-export for `load_and_validate_samples` and `ValidationError`. |
| `encode_pipeline.samples.loader` | ~105 | `load_and_validate_samples` orchestration: TSV parsing, required-column checks, duplicate-ID tracking, Pass 1/2/3 calls, and strict-input gating. |
| `encode_pipeline.samples.rows` | ~222 | Per-row normalization and validation (`validate_and_build_sample`). |
| `encode_pipeline.samples.strict` | ~79 | Strict-input file-existence validation (`validate_strict_inputs`, Bowtie2 index checks). |
| `encode_pipeline.samples.replicates` | ~271 | Cross-replicate validation (`validate_replicate_groups`) and IDR eligibility gates. |
| `encode_pipeline.samples.validate` | ~54 | Convenience `validate_sample_row` wrapper for tests/CLI. |
| `encode_pipeline.samples.models` | (existing) | Typed `SampleRecord` dataclass. |

All extracted modules depend only on `encode_pipeline.errors.ValidationError`,
`encode_pipeline.config.defaults`, and (where necessary) stdlib/os. None import
`validator.py`.

---

## 3. Remaining responsibilities in `validator.py`

### Public API (intentionally retained for backward compatibility)

| Name | Kind | Notes |
|---|---|---|
| `ValidationError` | exception alias | Points to `encode_pipeline.errors.ValidationError` via `config.errors`. |
| `validate_config()` | orchestrator | ~210 lines. Top-level config normalization and stage-gate checks. |
| `load_and_validate_samples` | lazy re-export | Resolved on first attribute access to `encode_pipeline.samples.load.load_and_validate_samples`. |
| `validate_replicate_groups` | lazy re-export | Resolved on first attribute access to `encode_pipeline.samples.replicates.validate_replicate_groups`. |
| `validate_picard_reference_resources()` | wrapper | Delegates to `encode_pipeline.config.genome`. |
| `validate_tss_annotation_resources()` | wrapper | Delegates to `encode_pipeline.config.genome`. |
| `main()` | CLI entry | `scripts/validate_samples.py` entry point. |

### Private helpers

| Name | Kind | Notes |
|---|---|---|
| `_coerce_int()` | wrapper | Delegates to `encode_pipeline.config.coercion.coerce_int`. |
| `_validate_effective_genome_size()` | wrapper | Delegates to `encode_pipeline.config.genome`. |
| `_validate_genome_resources()` | wrapper | Delegates to `encode_pipeline.config.genome`. |
| `_validate_tool_params()` | wrapper | Delegates to `encode_pipeline.config.tools`. |
| `_validate_idr_settings()` | wrapper | Delegates to `encode_pipeline.config.reproducibility`. |
| `_validate_reproducibility()` | wrapper | Delegates to `encode_pipeline.config.reproducibility`. |
| `_validate_cuttag_config()` | wrapper | Delegates to `encode_pipeline.config.cuttag`. |
| `_validate_mnase_config()` | wrapper | Delegates to `encode_pipeline.config.mnase`. |
| `_validate_qc_config()` | wrapper | Delegates to `encode_pipeline.config.qc`. |
| `_load_yaml()` | CLI helper | Prefers PyYAML, falls back to `_parse_config_minimal`. |
| `_parse_config_minimal()` | CLI helper | ~148-line stdlib YAML parser for standalone CLI use. |

### Stale / deprecated aliases (cleanup candidates)

| Alias | Points to | Status |
|---|---|---|
| `_SAMPLE_ID_RE` | `defaults.SAMPLE_ID_RE` | Deprecated; kept for backward compatibility. |
| `_SANITIZE_RE` | `defaults.SANITIZE_RE` | Deprecated; kept for backward compatibility. |
| `_BT2_STANDARD` | `defaults.BT2_STANDARD` | Deprecated; kept for backward compatibility. |
| `_BT2_LARGE` | `defaults.BT2_LARGE` | Deprecated; kept for backward compatibility. |

These aliases are no longer referenced inside the codebase but are left in place
in case external code imports them. A future PR can remove them after a
deprecation search across the repo and any documented public-API audit.

---

## 4. Import and dependency shape

### Current import graph

```
encode_pipeline.config
    └── encode_pipeline.config.validate
            ├── encode_pipeline.config.errors (ValidationError)
            └── encode_pipeline.config.validator
                    ├── encode_pipeline.config.coercion
                    ├── encode_pipeline.config.cuttag
                    ├── encode_pipeline.config.defaults
                    ├── encode_pipeline.config.errors
                    ├── encode_pipeline.config.genome
                    ├── encode_pipeline.config.mnase
                    ├── encode_pipeline.config.qc
                    ├── encode_pipeline.config.reproducibility
                    ├── encode_pipeline.config.tools
                    └── (lazy) encode_pipeline.samples.load
                        └── encode_pipeline.samples.replicates
```

### Circular-import mitigation

PR69 moved the canonical `ValidationError` to `encode_pipeline.errors` so that
`samples.*` modules can import it without pulling in `config.validator`.
`config.validator` re-exports `load_and_validate_samples` and
`validate_replicate_groups` lazily via module-level `__getattr__` to avoid a
partially initialized `samples.load` during `config/__init__.py` import.

### `config/__init__.py` eager-import risk

`encode_pipeline.config.__init__.py` currently does:

```python
from encode_pipeline.config.validate import (
    validate_config,
    validate_picard_reference_resources,
    validate_tss_annotation_resources,
)
```

Because `config.validate` imports `config.validator`, importing
`encode_pipeline.config` eagerly pulls in the entire validator module (and,
indirectly, all config-block submodules). This is functional today but couples
the package import to a large surface. A future PR should evaluate and likely
lighten this: either import only `ValidationError` and selected thin helpers,
or make the top-level `validate_config` import lazy. Any such change must be
preceded by characterization tests for `import encode_pipeline.config` and the
public names it is expected to expose, because workflow/Snakefile code and
standalone scripts rely on `from encode_pipeline.config import validate_config`.

---

## 5. Test safety net

### Sample-sheet characterization tests

- `test/samples/test_boundary.py` — shared `ValidationError` identity,
  backward-compatible re-exports, forbidden-import checks, default
  `error_cls`, subprocess fresh-import tests.
- `test/samples/test_replicates.py` — 42+ direct tests for
  `validate_replicate_groups` covering grouping, consistency, duplicate
  replicate combos, control consistency, Stage 5 / Stage 55 / Stage 64 /
  Stage 65 IDR eligibility gates.

### Config-block characterization tests

- `test/config/test_coercion.py`
- `test/config/test_genome_resources.py`
- `test/config/test_tool_parameters.py`
- `test/config/test_reproducibility.py`
- `test/config/test_qc.py`
- `test/config/test_cuttag.py`
- `test/config/test_mnase.py`

### Helper unit tests

- `test/config/test_*_helpers.py`

### Integration / smoke / behavior guards

- `test/config/test_validation.py`
- `test/test_stage8_smoke_profiles.py`
- `scripts/validate_samples.py --config config/config.yaml`
- `test/test_no_hardcoded_paths.py`
- `test/check_snakemake_lint.py`

---

## 6. Risk-ranked next-PR roadmap

| PR | Goal | Risk | Safety net |
|---|---|---|---|
| **PR73** | Pin `_load_yaml` / `_parse_config_minimal` / `main()` CLI behavior. | Low-Medium | Add tests that exercise the standalone CLI with the stdlib YAML fallback and with PyYAML, snapshot the normalized output, and assert exit codes. |
| **PR74** | Extract YAML parsing / CLI entry helper to a lightweight module (e.g. `encode_pipeline.cli.validate` or `encode_pipeline.config.yaml_loader`). | Medium | Characterization tests from PR73 must pass unchanged. Verify `scripts/validate_samples.py` still works. Check import-time side effects. |
| **PR75** | Lighten `config/__init__.py` import behavior. | Medium | First add characterization tests for `import encode_pipeline.config` and `from encode_pipeline.config import validate_config`. Then make the change and rerun the import tests, smoke profiles, and the full test suite. |
| **PR76** | Reduce `validator.py` to `validate_config` orchestrator only. Remove thin wrappers if they are no longer used, and delete stale aliases after a deprecation search. | Medium-High | Full test suite, smoke profiles, DAG snapshots, and standalone CLI must pass. Search the repo for any imports of `_SAMPLE_ID_RE`, `_SANITIZE_RE`, `_BT2_STANDARD`, `_BT2_LARGE`, and the wrapper functions before deleting. |

### Recommended order

1. **PR73/PR74** — CLI/YAML extraction is mostly isolated and reduces
   `validator.py` by ~165 lines.
2. **PR75** — After CLI/YAML is out of `validator.py`, the import-shape risk is
   smaller and easier to characterize.
3. **PR76** — Final cleanup once the orchestrator is the only remaining
   substantial block and all consumers are covered by tests.

---

## 7. Explicit non-goals and out-of-scope items

- No production code changes in this checkpoint PR.
- No changes to validation behavior, error messages, warning text, or defaults.
- No schema, Snakefile, rule, CI, lockfile, container, profile, or scientific
  changes.
- No migration of legacy `test/test_stage53*.py` / `test/test_stage55*.py`
  scripts in this phase.
- No new `encode_pipeline.config.cli` module unless a future design explicitly
  approves it.
- **`docs/superpowers/` and `research/` are out of scope** for this PR and the
  immediate decomposition roadmap. Any planning artifacts in those directories
  are not part of the production-audit deliverable.

---

## 8. Summary

`validator.py` has shrunk from ~2,070 lines to **662 lines**. The sample-sheet
phase (PR68–PR71) successfully moved loading, row validation, strict-input
checks, and replicate-group validation into focused `encode_pipeline.samples.*`
modules while preserving all public APIs through lazy re-exports. The remaining
large blocks are `validate_config()` (~210 lines) and the CLI/YAML helpers
(~165 lines). The next safe extraction is the YAML/CLI surface, followed by an
import-shape cleanup and a final orchestrator-only consolidation.
