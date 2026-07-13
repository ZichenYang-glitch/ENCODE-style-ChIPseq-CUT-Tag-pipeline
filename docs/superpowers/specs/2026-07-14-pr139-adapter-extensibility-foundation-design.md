# PR139 Adapter Extensibility Foundation Design

## Scope

PR139 establishes a reusable, workflow-neutral conformance gate for workflow
adapters. It does not add a second scientific adapter, dynamically discover
packages, change the worker lifecycle, or claim that an arbitrary external
workflow can already use the bundled local execution command composition.

The deliverable is deliberately split between two trust boundaries:

1. `WorkflowRegistry` rejects adapter declarations that can be checked without
   workflow inputs or filesystem access.
2. A packaged, pytest-independent conformance runner exercises adapter behavior
   with adapter-owned fixtures and rejects declarations that do not match the
   returned `Result` and platform-neutral value types.

The ENCODE adapter and a test-only minimal adapter run through the same suite.
The minimal adapter also passes through existing registry, workflow-info,
validation, workspace-planning, artifact-indexing, and QC-indexing services
without adding workflow-specific branches to those services.

## Current architecture assessment

The current structural `WorkflowAdapter` protocol and constructor-injected
`WorkflowRegistry` are sufficient for deterministic registration, metadata,
schema, validation, workspace planning, and adapter-owned artifact/QC mapping.
Platform services already resolve by `workflow_id`; they do not import ENCODE
config keys, sample columns, output paths, or QC metric names.

Two existing deployment boundaries remain intentionally explicit:

- `CommandBuilder` composes the bundled ENCODE Snakemake entrypoint rather than
  delegating the adapter `build_command` method.
- `WorkflowBuildIdentityProvider` fingerprints this repository's controlled
  source tree.

Those boundaries mean PR139 is an adapter contract and service-conformance
foundation, not a promise that an independently packaged workflow can enter the
current durable worker without later command/source-bundle composition work.
The conformance suite tests the adapter command contract honestly: a declared
`command` capability must return a `CommandSpec`; an adapter that does not
declare it must return a structured failure. It does not route the test-only
adapter through the bundled ENCODE `CommandBuilder`.

## Canonical capability vocabulary

The platform publishes one bounded vocabulary for the current contract:

- `validation`
- `dag_preview`
- `workspace_plan`
- `command`
- `input_authoring`
- `artifact_extract`
- `qc_summary_extract`

`WorkflowCapabilities` preserves declaration order but rejects duplicates,
invalid tokens, and unknown names. Extending the vocabulary therefore requires
an intentional platform contract change instead of silently exposing an inert
or misspelled capability through the API.

Registration additionally requires `WorkflowCapabilities` as the concrete
declaration type. The optional `qc_summary_extract` label and the
`QcSummaryExtractingAdapter` protocol must agree in both directions. Other
methods are part of the base structural protocol and need fixtures to test, so
their dynamic consistency belongs to the conformance runner rather than
registry construction.

`WorkspacePlanner` will check `workspace_plan` before delegation, matching the
existing validation, artifact, and QC service behavior. Unsupported capability
paths return stable platform issues and do not call adapter code.

## Reusable conformance runner

`encode_pipeline.testing.adapter_conformance` is installed with the package but
does not import pytest. A future adapter package can construct an
`AdapterConformanceCase` in its own test suite and call
`verify_adapter_conformance(case)`.

Each case supplies:

- an adapter instance;
- valid and invalid `WorkflowInputs`;
- a fresh planning workspace path that must remain untouched by pure planning;
- an adapter-prepared artifact workspace;
- platform-vetted QC source documents.

The runner verifies:

- concrete metadata/capability declarations and registry round-trip;
- stable, fresh, JSON-safe `WorkflowSchema` responses;
- validation success and failure behavior;
- declared versus unsupported DAG, workspace, command, and artifact behavior;
- platform-neutral result value types and tuple collections;
- exact optional QC protocol/capability agreement, unique source types, and
  platform-neutral QC candidate types;
- structured failure results for undeclared callable capabilities.

The runner raises `AdapterConformanceError` with a bounded contract-coordinate
message. It never includes input contents, paths, exception text, or adapter
private return values. Adapter-controlled property, method, schema-serialization,
and JSON Schema failures pass through one safe invocation boundary; the public
error suppresses the original exception chain so ordinary tracebacks cannot
recover private exception details.

## Test-only minimal adapter

The minimal adapter lives only under `test/`. It implements every current
capability using deterministic in-memory values and a tiny workspace fixture.
It contains no ENCODE fields or imports. Its tests register it alongside normal
platform services and prove that metadata, schema, validation, workspace,
command, artifact, and QC contracts are exercised without adding a workflow ID
branch to production business logic.

The ENCODE case uses the existing tiny deterministic input profile, a
materialized adapter workspace, and bounded QC fixtures. Existing scientific
adapter tests remain authoritative for detailed manifest vocabulary and QC
semantics; the shared suite checks only the cross-adapter contract.

## Discovery decision

PR139 does not add entry-point discovery. Explicit registry construction is the
current deployment allowlist shared by API and worker processes. Loading entry
points now would add package selection, duplicate/version policy, import-failure
isolation, and build-identity questions without a real external adapter.

When the first external RNA-seq package is implemented, discovery should use
Python standard `importlib.metadata.entry_points` with a named entry-point group
and explicit allowlisting. No custom module scanner or plugin loader is needed.
Hi-TrAC/TracPre2 follows RNA-seq so the discovery contract is first exercised by
the simpler workflow family.

## Failure and security semantics

- Static declaration corruption fails registry construction.
- Dynamic declaration/behavior mismatch fails conformance before release.
- Unsupported workspace planning fails before adapter invocation.
- Platform artifact and QC services continue to validate lifecycle, build
  identity, paths, file types, sizes, candidate types, and persistence.
- The conformance runner does not weaken any runtime safety boundary and does
  not treat an adapter as trusted filesystem or lifecycle code.

## Non-goals and residual limits

PR139 does not implement RNA-seq, Hi-TrAC/TracPre2, entry-point loading,
immutable workflow bundles, generic runtime command dispatch, authentication,
multi-user isolation, HPC/Kubernetes, object storage, or Agent write actions.
The bundled local command and source fingerprint remain ENCODE deployment
composition that a later real external adapter must address explicitly.
