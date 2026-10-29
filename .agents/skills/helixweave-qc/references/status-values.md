# QC status values

Use the frozenset that belongs to the column being read. Do not add members,
and do not add pass/fail cutoffs. Checklist numbers in
`docs/qc-interpretation.md` ("Practical Checklists", "Known Caveats") are
guidelines. Manifest `qc_flag` is always `NA`.

## Distinctions

These are the distinctions in
[../../_shared/entries.md](../../_shared/entries.md). The set that actually
allows a number is the frozenset for that column.

- numeric `0`: a measured zero only when the accompanying status says the
  value is defined (`ok`, `reported`, or another status that contract code
  allows to carry a number)
- `NA`: the cell is absent or undefined. Do not coerce it to zero
- `not_applicable`: the metric does not apply (for example SE template
  metrics, or SE `max_insert`)
- `not_enabled`: the producing switch is off
- `missing`: the producer was expected and the input is absent
- `parse_failed`: the source was present but rejected by the strict reader
- `no_data`, `tool_failed`, `unresolved_contigs`, `no_mapped_reads`:
  distinct contract statuses; read the column's allowed set before renaming
  them

`not_collected` and `no_valid_templates` are members of policy and template
sets below. They are not aliases of `missing` or of numeric `0`.

## Where each family is defined

| Set | Defined in |
| --- | --- |
| `ANALYSIS_TEMPLATE_METRICS_STATUSES` | `src/encode_pipeline/qc/contract.py` |
| `POLICY_ANALYSIS_TEMPLATE_STATUSES` | `src/encode_pipeline/qc/contract.py` |
| `POLICY_ORGANELLE_STATUSES` | `src/encode_pipeline/qc/contract.py` |
| `POLICY_CROSS_CORRELATION_STATUSES` | `src/encode_pipeline/qc/contract.py` |
| `CROSS_CORRELATION_EXECUTION_STATUSES` | `src/encode_pipeline/qc/contract.py` |
| `CROSS_CORRELATION_INTERPRETATION_STATUSES` | `src/encode_pipeline/qc/contract.py` |
| `CROSS_CORRELATION_THRESHOLD_HINTS` | `src/encode_pipeline/qc/contract.py` |
| `CROSS_CORRELATION_TARGET_TYPES` | `src/encode_pipeline/qc/contract.py` |
| `POLICY_DUPLICATE_POLICIES` | `src/encode_pipeline/qc/contract.py` |
| `POLICY_BLACKLIST_MODES` | `src/encode_pipeline/qc/contract.py` |
| `POLICY_ORGANELLE_MODES` | `src/encode_pipeline/qc/contract.py` |
| `_SAMPLE_STATUSES` | `src/encode_pipeline/adapters/encode_qc.py` |
| `_CROSS_CORRELATION_STATUSES` | `src/encode_pipeline/adapters/encode_qc.py` |
| `_CC_VALUE_STATUSES` | `src/encode_pipeline/adapters/encode_qc.py` |
| `_MNASE_ANALYSIS_STATUSES` | `src/encode_pipeline/adapters/encode_qc.py` |
| `_ALLOWED_INPUT_STAGES` | `src/encode_pipeline/qc/policy.py` |

`_ALLOWED_INPUT_STAGES` names denominators, not metric statuses. Input-stage
keys are `ANALYSIS_POLICY_INPUT_STAGE_KEYS` and
`ANALYSIS_POLICY_V2_INPUT_STAGE_KEYS` in `contract.py`.

Sample-summary `cross_correlation_status` uses `_CROSS_CORRELATION_STATUSES`.
The detail TSV uses `CROSS_CORRELATION_EXECUTION_STATUSES` and
`CROSS_CORRELATION_INTERPRETATION_STATUSES`. Do not rename one family into
the other. `CROSS_CORRELATION_THRESHOLD_HINTS` is the closed set for
`threshold_hint`. `validate_cross_correlation_detail` in `contract.py`
decides when a hint may be other than `not_evaluated`. A hint is not
biological acceptance of the run. `CROSS_CORRELATION_TARGET_TYPES` is the
configured selector. On the detail row, `validate_cross_correlation_detail`
allows `target_type` `tf`, `histone`, or `unknown`.

`peak_mode_status` on pooled QC is not a named frozenset. `encode_qc.py`
accepts `ok`, `mismatch`, and `unknown` and rejects anything else as
`EncodeQcSummaryError`. A platform metric flag
(`_ALLOWED_FLAGS` in `src/encode_pipeline/services/qc_summary_indexing.py`)
is not manifest `qc_flag`.

Fragment-size `status` has no frozenset in `contract.py`. Keep the cell
text. An unrecognized status is not a zero.

## What NA means

`NA` means the cell is absent or undefined. It is not a measured zero.

On a version 2 sample row, `encode_qc.py` requires `organelle_reads`,
`organelle_total_mapped_reads`, and `organelle_fraction` to be `NA`
together unless `organelle_status` is `ok`, and all three defined when it
is `ok`. NSC and RSC stay `NA` when the cross-correlation status is outside
`_CC_VALUE_STATUSES`. When every `project_qc_summary` input is v1, the v1
header is preserved; if any input is v2, v1 rows are padded with `NA` for
fields that version did not define (`docs/output-contract.md`, "QC Contract
(v2)"). That padding is an undefined field, not a measurement.

JSON `null` in `stage_counts` is the same kind of absence. `policy.py`
rejects a non-null `organelle_reads` when `organelle_status` is not `ok`.

## What parse_failed means

`parse_failed` means a source was present and the strict reader rejected it.
`EncodeQcSummaryError` and `AnalysisPolicyError` fail the document. They do
not become zeroes.

Unknown `schema_version`, extra columns, duplicate columns, a wrong numeric
type, and a status outside the column's frozenset are parse failures. For
analysis-template and MNase analysis counts, `parse_failed` is carried with
`NA` counts. A provided but unreadable analysis-template source is
`parse_failed`, not a zero (`scripts/write_analysis_policy.py`).
`docs/output-contract.md` says a missing or malformed analysis BAM is
`parse_failed` rather than a fabricated zero. Keep `not_collected` when the
policy record uses it: the writer uses that status when no template
component and no legacy source were provided.

## Manifest statuses

`results/multiqc/result_manifest.tsv` `status` is `present`, `missing`, or
`not_applicable` (`docs/output-contract.md`, "Manifest Field Schema").
`qc_flag` is always `NA`. Automatic threshold conclusions are not
implemented.

Manifest `missing` means the path was absent when the manifest was written.
Metric status `missing` means the producer was expected and its input was
absent. Manifest `not_applicable` means that product does not apply to the
row. Do not treat either manifest word as a measured zero.

## Same word on two records

Read the status on the record being cited. `scripts/assemble_qc_summary.py`
records organelle status `missing` when the mito table is absent.
`scripts/write_analysis_policy.py` records `not_applicable` when that path
was not provided, and `parse_failed` when the file was present but invalid.
Do not rename one record's status into the other's.

## When zero is defined

A zero is a measurement only for a status the contract allows to carry a
number.

- Organelle fraction `0` is defined for status `ok` when mapped reads are
  positive and organelle reads are zero. `unresolved_contigs` and
  `no_mapped_reads` keep `organelle_reads` and `fraction` at `NA`
  (`scripts/calc_organelle_metrics.py`). On `mito_fraction.tsv`,
  `total_mapped_reads` is always the idxstats total the writer stored:
  `0` for `no_mapped_reads`, and often a positive count for
  `unresolved_contigs`. That cell is not an organelle fraction. The
  sample-summary reader still requires all three organelle value columns
  to be `NA` unless status is `ok`.
- `no_valid_templates` reports a legitimate zero valid-template count.
  `ok` requires a positive valid count. `not_applicable`, `not_collected`,
  and `parse_failed` keep those counts null (`policy.py`).
- SE template metrics and SE `resolved_max_insert` use `not_applicable`
  (`POLICY_MAX_INSERT_NOT_APPLICABLE`), with `NA` counts where the template
  contract says so.
