# Workflow Adapter Conformance

The adapter conformance suite is the release-time contract gate for workflow
adapters. It complements runtime platform validation; it does not make adapter
filesystem output trusted.

## Use the suite

Install the project development dependencies so JSON Schema 2020-12 validation
is available, then create one adapter-owned test fixture:

```python
from encode_pipeline.testing import (
    AdapterConformanceCase,
    verify_adapter_conformance,
)


def test_adapter_conformance(adapter_case: AdapterConformanceCase) -> None:
    verify_adapter_conformance(adapter_case)
```

`AdapterConformanceCase` requires:

- the adapter instance;
- valid and invalid `WorkflowInputs`;
- a fresh absolute path for pure workspace planning;
- an adapter-prepared absolute artifact workspace;
- any platform-vetted `QcSourceDocument` fixtures.

The fixture owns scientific input and output meaning. The shared runner owns
only workflow-neutral contract checks.

Adapter-controlled descriptors, methods, and schema serialization are invoked
through a bounded diagnostic boundary. Failures expose only the contract
coordinate and do not chain the adapter's original exception into a traceback.

## Capability rules

The current platform vocabulary is:

- `validation`
- `dag_preview`
- `workspace_plan`
- `command`
- `input_authoring`
- `artifact_extract`
- `qc_summary_extract`

Names are exact, unique lowercase tokens. Unknown names fail adapter
construction rather than appearing as inert API promises.

For a declared callable capability, the adapter must return a successful
`Result` containing the corresponding platform value type. For an undeclared
base-protocol capability, the method must return a structured failure. The QC
capability is optional at the protocol level, so its declaration and
`QcSummaryExtractingAdapter` implementation must agree.

All adapters exposed through the current authoring platform declare
`input_authoring` and return a fresh, stable, valid JSON Schema 2020-12 contract.

## What the suite proves

The suite checks metadata, schema, validation, DAG preview, workspace planning,
command construction, artifact candidates, QC candidates, and capability
truthfulness. A test-only minimal adapter proves that the registry,
workflow-info, validation, workspace-planning, artifact-indexing, and
QC-indexing services do not require an ENCODE workflow-ID branch.

Artifact and QC platform services still enforce lifecycle state, build
identity, path containment, regular-file and size rules, candidate types, and
atomic persistence. Conformance never bypasses those runtime checks.

## Current deployment limit

The local `CommandBuilder` and `WorkflowBuildIdentityProvider` remain composed
for this repository's bundled ENCODE Snakemake source. The tested conformance
seam does not imply that an arbitrary external package can enter the durable
worker without a later, explicit source/command composition design.

Dynamic discovery remains deferred until a second adapter and its deployment
boundary are approved. Bulk RNA-seq is the current research candidate, not a
committed adapter or delivery date. If external discovery is later selected,
use Python standard `importlib.metadata.entry_points` with an explicit
deployment allowlist and deterministic duplicate/version and load-failure
policy; do not add a custom module scanner.
