# PR142 Semantic ENCODE Config Design

## Objective and boundary

PR142 replaces historical engine-stage switches in the default ENCODE
authoring experience with two stable scientific concepts:

```yaml
replicate_analysis:
  enabled: true

chipseq_idr:
  enabled: false
```

The adapter schema advances from `1.0.0` to `1.1.0`. Existing YAML using the
legacy engine keys remains valid, and the Snakemake workflow, rule names,
targets, output paths, and scientific validator remain unchanged. This PR does
not add QC sources, another adapter, or platform-wide config normalization.

## Baseline findings

The public config schema currently exposes the two legacy keys directly. The
scientific validator and Snakefile consume those keys, while validated input
snapshots persist the original `WorkflowInputs`. Workspace planning calls the
adapter again and renders the adapter-private validated config. That existing
separation provides the required compatibility boundary without changing the
platform or persistence layers.

The advanced `reproducibility` block is a broader policy covering consensus
and multiple assay/peak-mode IDR paths. `chipseq_idr` is retained because it
names the narrower ChIP-seq narrow-peak IDR gate that the legacy engine switch
controls; it neither aliases nor overrides the broader policy block. Its
description makes the eligible two-biological-replicate scope explicit.

## Considered approaches

### Teach the scientific validator the semantic fields

This would make the mature validator and CLI own authoring aliases, broaden the
scientific contract, and risk changing non-platform callers. It is rejected.

### Normalize in the API or snapshot service

This would make workflow-neutral platform code understand ENCODE keys and
would replace the user-authored semantic document in the immutable snapshot.
It is rejected.

### Translate an adapter-owned validation copy

The ENCODE adapter deep-copies the submitted config, validates the two semantic
objects, resolves legacy/semantic compatibility, removes semantic fields from
that private copy, and supplies strict boolean legacy keys to the existing
scientific validator. This is selected. Snapshot creation still canonicalizes
the untouched submitted `WorkflowInputs`; workspace planning renders only the
validated engine config.

## Public schema 1.1.0

Both semantic properties are strict object schemas with a required `enabled`
boolean and `additionalProperties: false`. Defaults are complete objects:
replicate analysis is enabled and ChIP-seq IDR is disabled. Titles and
descriptions are scientific and contain no engine-stage numbering. The former
legacy properties are absent from `properties`; top-level
`additionalProperties: true` remains necessary for partial-schema advanced
YAML and legacy compatibility.

All three document IDs end in `/1.1.0`. Coverage, modes, and uniform input
ceilings are unchanged. The generic workbench explicitly supports the already
published `1.0.0` contract and the new `1.1.0` contract while continuing to
reject unknown versions. This keeps an independent adapter's compatible
`1.0.0` authoring schema usable instead of coupling every adapter revision to
the ENCODE schema revision. Its default form, deterministic Review payload,
controlled browser fixture, screenshots, and request assertions use only
semantic fields. A user who deliberately imports legacy YAML may still see the
exact keys they supplied in the transparent request preview; the product does
not generate or label those fields.

## Adapter translation and issues

Translation handles each semantic/legacy pair independently:

1. Neither supplied: apply the existing engine default.
2. Semantic only: require exactly `{enabled: <bool>}` and set the private
   engine key to that value.
3. Legacy only: pass it unchanged to the existing scientific validator.
4. Both supplied and equivalent after the existing `coerce_bool`
   normalization used by the scientific config layer: accept and emit one
   `ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED` warning. The compatible legacy set
   is boolean or case-insensitive, non-trimmed `"true"`/`"false"`; null,
   numbers, `yes`, and whitespace-padded strings remain invalid.
5. Both supplied and different: fail closed with
   `ENCODE_CONFIG_SEMANTIC_CONFLICT`.

Malformed semantic objects fail with `ENCODE_CONFIG_SEMANTIC_INVALID` before
scientific validation. Non-mappings, missing `enabled`, null, string booleans,
unknown nested keys, and non-boolean values are rejected. Public issues use
controlled messages and semantic paths only; they contain no submitted value,
legacy key, exception text, path, or environment data. An invalid legacy value
continues through the existing `ENCODE_CONFIG_INVALID` path instead of being
misclassified as a semantic conflict. If both pairs use a deprecated alias,
the warning is deduplicated into one stable issue.

The translated private config is the only value passed to
`validate_config()`. Existing dependencies still apply, including the rule
that ChIP-seq IDR cannot be enabled when replicate analysis is disabled.
Scientific failures retain the established `ENCODE_CONFIG_INVALID` envelope.
Validation warnings flow through `ValidationService`, the validation API, and
snapshot validation evidence. Workspace planning preserves the warning and
adds its existing planning info issue. API tests require status 200, one
warning, a usable snapshot, semantic-only path data, and no legacy key or value
in `technical_message` or context.

The input workbench retains Issue severity. Successful validation warnings use
an amber `role=status` treatment and do not disable the validated snapshot or
Create action; errors retain the red alert treatment. This is a local fix in
`ValidatedSubmission`, not a global issue-display redesign.

RJSF renders the canonical draft with `emptyObjectFields: skipDefaults`.
Initial drafts still receive adapter defaults once through
`createDefaultObject`, but importing a legacy-only or partial advanced YAML
document does not cause a later Form edit to inject fields the user never
submitted. A regression uses legacy values opposite the semantic defaults,
edits a known Form field, and requires YAML, Review, and validation to retain
the exact legacy switches without adding semantic aliases.

## Snapshot, workspace, and build identity

`ValidatedInputService` already snapshots the caller's original inputs after
adapter validation. Tests will require a semantic-only snapshot payload and no
injected legacy engine fields. A run reconstructed from that snapshot is
validated and translated again during workspace planning. The materialized
`config/config.yaml` contains strict boolean engine keys and no semantic
authoring objects, so the current scientific validator and Snakemake DAG see
the same configuration as before.

The build identity already fingerprints all Python files under
`src/encode_pipeline`, including the schema and translator. Therefore schema
and translation changes are bound to validation, preflight, and worker
execution without a new identity scheme or migration.

## Browser gate and compatibility evidence

The controlled results fixture converts its copied legacy tiny-profile
switches into semantic authoring objects before exposing config to the browser.
The underlying scientific tiny profile remains unchanged. The real Playwright
path proves authoring, validation, immutable snapshot creation, run creation,
preflight, explicit start, RUNNING, and SUCCEEDED. It asserts that both the
validate request and the visible default workbench contain no legacy stage
keys. The default Config form is captured at 1440x900, 1024x768, 390x844, and
360x800 and must show accessible `Replicate analysis` and `ChIP-seq IDR`
labels with the correct defaults and no horizontal overflow. The validate
request must contain both strict semantic objects and no legacy tokens.

Advanced-YAML compatibility gets a separate state test: user-supplied legacy
fields survive YAML -> Form -> YAML and Review unchanged. That transparent
preview is explicitly distinguished from the default semantic Review so the
workbench never hides submitted request content.

Unit and integration tests cover schema structure, every boolean mapping pair,
legacy-only and semantic-only inputs, equal aliases, conflicts, malformed
objects, warning propagation, original snapshot payload, workspace engine
config, and a scientific-validator round trip. The compatibility matrix also
covers legacy bool/string/case spellings, invalid legacy values, and the
existing orthogonality of `chipseq_idr.enabled` and advanced
`reproducibility.idr.chipseq_narrow`. Existing legacy scientific DAG tests
remain unchanged and supply compatibility evidence.

## Non-goals and residual risks

The partial schema still allows advanced YAML properties by design, so it
cannot prohibit a knowledgeable user from supplying legacy keys. Hiding those
user-supplied values in Review would make the request preview dishonest and is
not attempted. The advanced `reproducibility` policy remains a distinct legacy
scientific surface and is not redesigned. Immutable workflow bundles, QC
expansion, another adapter, authentication, multi-user isolation, HPC, object
storage, and Agent writes remain out of scope.
