# Known Issues and Follow-ups

This file records current, non-blocking limitations. Completed implementation
history belongs in `CHANGELOG.md`, Git history, and `docs/release-checks/`.
The maintained delivery sequence is the
[HelixWeave product roadmap](docs/development/workflow-platform-agent-roadmap.md).

## Scientific scope

HelixWeave ships an ENCODE-inspired workflow; it does not claim complete ENCODE
pipeline parity for every assay, mark, replicate design, or QC metric.

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
- Authentication, multi-tenant isolation, PostgreSQL, object storage, HPC
  schedulers, Kubernetes, and remote workspace semantics are not implemented.
- SQLite is canonical lifecycle and result-metadata state. Queue state must not
  be used as an alternative source of truth.
- Artifact downloads are limited to indexed workspace artifacts; arbitrary
  paths and symlink escapes are intentionally rejected.
- The Agent surface is read-only. It cannot submit, start, cancel, modify, or
  delete runs, and generated explanations are not provenance.

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

## Reporting a new issue

A useful issue should identify the assay or platform layer, minimal inputs,
expected and observed behavior, exact command or API operation, environment,
and whether the problem reproduces with a committed test profile. Do not attach
secrets, private sample payloads, absolute workspace paths, or large data files.
