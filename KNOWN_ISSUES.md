# Known Issues and Follow-ups

This file records current, non-blocking limitations. Completed implementation
history belongs in `CHANGELOG.md`, Git history, and `docs/release-checks/`.
The maintained delivery sequence is the
[HelixWeave product roadmap](docs/development/workflow-platform-agent-roadmap.md).

## Scientific scope

HelixWeave ships two scientific adapters. Their supported paths and evidence
must not be generalized beyond the contracts below.

### ENCODE-style epigenomics

The ENCODE-inspired workflow does not claim complete ENCODE pipeline parity
for every assay, mark, replicate design, or QC metric.

- Broad-peak IDR is experimental and opt-in. It is not the default final-output
  policy. See [reproducibility policy](docs/reproducibility-policy.md).
- The opt-in broad-peak IDR warning says consensus remains primary, while the
  current DAG and result manifest select the broad-IDR artifact as the final
  reproducibility output. This pre-existing policy/implementation mismatch
  needs an explicit scientific decision before either behavior is changed.
- SEACR CUT&Tag outputs use consensus behavior; SEACR IDR has no supported rank
  scheme and is outside the current policy.
- IDR eligibility remains bounded by assay, peak mode, and replicate policy.
  Larger or irregular replicate designs must use the documented consensus
  behavior rather than inferred IDR semantics.
- ATAC broad peaks, footprinting, and ATAC-specific nucleosome positioning are
  not supported.
- MNase fragment classes, dyad/occupancy tracks, and QC are supported, but a
  dedicated nucleosome caller such as DANPOS3, iNPS, or SEM is not implemented.
- CUT&Tag spike-in normalization and a multi-caller comparison framework are
  not implemented.
- GC-bias metrics and optional `plotFingerprint` QC are not implemented.

### Bulk RNA-seq

- The adapter is pinned to nf-core/rnaseq 3.26.0. Authoring and validation are
  available without execution assets, but create and start remain fail-closed
  until the complete operator-owned runtime and transcriptome binding pass
  live admission.
- Nextflow, the JDK/plugin/source closure, OCI images, references, STAR/Salmon
  indexes, and SortMeRNA database/index assets are not bundled or downloaded
  automatically.
- The controlled synthetic acceptance fixture proves runtime, lifecycle,
  cleanup, artifact, QC, and product contracts. It does not establish
  biological validity or production-scale performance.
- RiboDetector/GPU, Tower/Wave, cloud/HPC executors, and cross-attempt resume
  are outside the current product boundary.

Assay-specific truth is maintained in [assay policy](docs/assay-policy.md),
[configuration](docs/configuration.md), and
[QC interpretation](docs/qc-interpretation.md).

## Reference data and external validation

- No real public FASTQ, BAM, BigWig, reference genome, or index is bundled.
- Public-data reports describe external executions and must not be interpreted
  as reproducible downloads unless their external data identifiers and
  environment are still available.
- Site-specific reference preparation, storage throughput, and scheduler
  behavior remain operator responsibilities.
- Existing real-data, container, and release evidence is retained under
  [`docs/release-checks/`](docs/release-checks/).

## Workflow operations

- Large FE/ppois bedGraph conversion can create temporary disk pressure.
  Limit concurrent conversions with the documented Snakemake resource and
  provide adequate scratch space.
- MultiQC sample labels can drift when upstream FASTQ names differ across data
  sources. A project-specific replacement-name mapping remains the recommended
  mitigation.
- Some optional QC rules tolerate unavailable optional outputs. Operators who
  require strict completeness should verify the result manifest and QC index,
  not only the terminal workflow status.
- Network filesystems may require a larger Snakemake latency wait and careful
  workspace placement.

See [container usage](docs/container-usage.md),
[environment guidance](docs/environments.md), and the
[local platform runtime](docs/development/local-platform-runtime.md).

## Platform scope

- The supported product boundary is local or small trusted teams with SQLite,
  Redis/RQ, and filesystem workspaces.
- The default registry contains the ENCODE-style epigenomics and Bulk RNA-seq
  adapters. An unavailable optional Bulk RNA-seq runtime does not prevent safe
  authoring or validation and does not make the base platform unavailable.
- Authentication, multi-tenant isolation, PostgreSQL, object storage, HPC
  schedulers, Kubernetes, and remote workspace semantics are not implemented.
- SQLite is canonical lifecycle and result-metadata state. Queue state must not
  be used as an alternative source of truth.
- Artifact downloads are limited to indexed workspace artifacts; arbitrary
  paths and symlink escapes are intentionally rejected.
- The Agent surface is read-only. It cannot submit, start, cancel, modify, or
  delete runs, and generated explanations are not provenance.
- Omics Intake integration is limited to a pinned, read-only Bundle 0.2
  inspection boundary for ENCODE authoring inputs. It does not import producer
  code, persist Bundle provenance, create a snapshot, or authorize execution.
- v0.3.0 does not publish a HelixWeave application, ENCODE runner, or
  Bulk RNA-seq container image.

See the [architecture overview](docs/architecture/platform-overview.md) for
the durable ownership and safety boundaries.

## Quality and delivery scope

The maintenance and quality baseline is complete. Deterministic tests now have
distinct PR-fast and full-main selections, coverage consumes the single pytest
producer's evidence, and platform, scientific, and container real-execution
jobs run independently on manual dispatch, nightly schedules, and releases.
The measured test inventory and enforced coverage floors are maintained in the
[quality baseline](docs/development/coverage-policy.md).

Real-execution jobs are intentionally not required pull-request checks. A
high-risk branch needs an explicit workflow dispatch when those environments
are relevant. CI timing budgets are review signals rather than pass/fail
thresholds, and configuring branch-protection contexts remains a maintainer
operation outside repository content. See the
[development harness](docs/development/harness.md) for the current tier model.

### Frontend dependency and bundle maintenance

The v0.3.0 candidate lock reports 12 npm package-level findings: three are in
production dependency paths (`fast-uri`, `react-router`, and
`react-router-dom`), while nine are confined to the Vite, Vitest, Orval, and
OpenAPI-generation toolchain. The audit did not classify any finding as a
v0.3.0 release blocker within the supported single-workstation or trusted-team
boundary: schema URI parsing is not used as a network allow/deny decision,
navigation targets are internal and constrained, there is no server-side
rendering, Vitest UI/API/browser mode is not enabled, and generated OpenAPI
input is repository-controlled. Vite is nevertheless part of the supported
local trial launcher and retains a browser-to-loopback attack surface. It must
remain bound to loopback and must not receive untrusted requests; exposing it
to a LAN or other untrusted network is outside the supported trial boundary.

Maintenance should update `fast-uri` to the compatible fixed `3.1.4` lock
version, then test coordinated React Router 7.18.1+, Vite 6.4.3+/Vitest 3.2.6+,
and Orval/JS-YAML upgrades as their upstream ranges permit. Do not use
`npm audit fix --force` or an unreviewed major upgrade merely to reduce the
aggregate count.

The lazy `/workflows/:workflowId/new-run` authoring chunk is approximately
1.20 MB minified (385 KB gzip). It is not loaded by the workflow catalog or
detail routes; most of its size comes from CodeMirror/Lezer and RJSF/AJV.
Future bounded work can align the duplicate `@codemirror/view` versions and
load the YAML editor only when YAML mode is selected. Raising Vite's warning
limit would only hide the maintenance signal and is not a fix.

## Reporting a new issue

A useful issue should identify the assay or platform layer, minimal inputs,
expected and observed behavior, exact command or API operation, environment,
and whether the problem reproduces with a committed test profile. Do not attach
secrets, private sample payloads, absolute workspace paths, or large data files.
