# Bulk RNA-seq Offline Execution Boundary

Status: Accepted for phased implementation. PR #151 enables only an explicitly
runtime-composed adapter; product exposure remains deferred.

## Context and decision

HelixWeave implements `bulk-rnaseq` as a Python `WorkflowAdapter` around the
released nf-core/rnaseq 3.26.0 pipeline. The scientific workflow remains
upstream DSL2/Nextflow; workflow-neutral platform code does not depend on
nf-core parameters, Nextflow process names, or RNA-seq output paths. No custom
pipeline Groovy is maintained by HelixWeave.

This controlled black-box boundary is preferred to a Snakemake rewrite because
reimplementing the mature upstream workflow would create substantial scientific
equivalence and upgrade risk. A thin wrapper that downloads a release at launch
is also rejected because network resolution cannot provide an auditable,
reproducible, offline build. nf-core/rnaseq remains MIT licensed; reference and
rRNA database licensing, and the licenses of third-party tool/container
contents, remain the operator's responsibility.

PR #150 fixed the versioned 1.0.0 authoring contract and its three ownership
layers:

- **Standard** contains stable typed product semantics, including explicit
  reference identity and SortMeRNA/Bowtie2 rRNA-removal configuration.
- **Advanced** is a reviewed allowlist validated against the immutable 3.26.0
  upstream schema. Unknown, conflicting, and output-shaping parameters fail
  closed.
- **Platform-owned** includes `input`, `outdir`, launch/work/cache directories,
  executor, profile, resume, reports, plugins, containers, networking, and all
  CLI or Nextflow configuration. Callers cannot supply raw argv, shell text, or
  Groovy configuration.

PR #151 does not expand those parameters. It adds the immutable, offline-first
execution boundary described below.

## Immutable runtime identity

An execution-capable adapter exists only when an operator supplies a complete
local asset binding. Every coordinate is verified before workspace planning or
build-identity capture; a missing, extra, replaced, symlinked, non-regular, or
digest-mismatched asset fails closed with a path-redacted doctor result.

- The pipeline is nf-core/rnaseq 3.26.0 at commit
  `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`. The committed 829-file manifest
  has SHA-256
  `dc75d105ad26b381197268ef67c44da3107b694e788a3c693237924a86ead774` and
  identifies source tree
  `4f779bd8934d41896fbf137ff31158b02daacd7efbb40ed8cee55b9c8f757722`.
  Execution uses that verified local tree, never a tag, Git URL, or nf-core
  launcher lookup.
- Nextflow is fixed to the official 25.04.3 distribution, build 5949, SHA-256
  `53c232cdd8a9419d2c205dc7c6c4dd2646182c997300e6439a453099e28aa21a`.
  The JDK remains a server-owned deployment prerequisite rather than a
  vendored scientific asset; a real version/configuration canary and every
  run's configuration preflight prove compatibility. Its patch identity is
  not currently part of the scientific build digest and remains an operator
  reproducibility constraint.
- The only pipeline plugin is `nf-schema` 2.5.1. Its official archive, metadata,
  and extracted local tree are independently hashed. `NXF_PLUGINS_DIR` points
  only to that verified local installation and `NXF_OFFLINE=true`; plugin
  resolution may not use the network.
- The 78-process pinned source universe is audited by source path and file
  digest: 56 processes are in the supported route and 22 have an explicit
  exclusion reason. A reserved negative-label selector assigns every upstream
  process a non-existent image ID before the 56 verified `withName` selectors
  override it. An unclassified or accidentally reachable process therefore
  fails locally instead of falling back to its source-level image tag.
- Container assignments in pinned Nextflow config are a separate precedence
  surface from process-module declarations. All 78 such assignments in the
  source tree are classified and digest-bound; 77 belong only to the unselected
  ARM profile or upstream tests. The sole default-loaded assignment targets
  the `PARABRICKS_STARGENOMEGENERATE` alias. Nextflow 25.04.3 gives that alias
  selector priority over both an original-process selector and the
  negative-label default. The platform therefore repeats the exact audited
  alias selector after its verified selectors and assigns the non-existent
  image ID. A pinned-engine priority canary proves that the later same-selector
  override wins. New or changed default config assignments fail the audit and
  require an explicit verified-image or deny decision.
- Each supported image retains its upstream coordinate and raw distribution
  manifest/OCI digest as provenance. Large images are not committed to Git. An
  operator-owned availability lock must match the inventory exactly and bind
  every process to a bounded Docker archive by size and SHA-256. The raw
  manifest binds the image config; the config and every rootfs diff ID bind the
  archive. Runtime execution uses the verified local config image ID because a
  Docker save/load round trip does not preserve `RepoDigests` reliably.
- Runtime policy disables pipeline, plugin, and container fetching, Wave,
  Tower, Fusion, and alternative package/container engines. Docker is launched
  against one verified local Unix socket with `--pull=never`, `--network=none`,
  and automatic removal; `docker.registry` is empty so Nextflow cannot qualify
  a config image ID as a registry reference. A lock records availability but
  never authorizes a load or pull.

The adapter build identity frames the adapter version, source manifest/tree,
Nextflow distribution, plugin archive/tree, container inventory and
availability lock, and the fixed upstream schemas. Changing any coordinate
produces a different durable build identity.

## Workspace and command boundary

Workspace planning first verifies FASTQ, FASTA, GTF, optional reference index,
and rRNA resources without modifying them. It then deterministically emits:

- a sorted nf-core samplesheet and canonical params JSON;
- a platform-owned Nextflow config with the local executor, exact container
  selectors, a default-deny container selector, explicit audited alias-level
  overrides, report paths, and no-network Docker policy;
- execution and cache identity documents; and
- isolated `engine/launch`, `engine/work`, `engine/cache`, `NXF_HOME`, temporary,
  log, report, and result directories.

The command uses fixed argv with `shell=False`: the pinned local Nextflow binary
runs the pinned local source path with the generated config and params file,
the Docker profile, `-offline`, the server-owned work directory, and a bounded
identity-derived run name. The launch directory is inside the run workspace.
An explicit command-owned `nextflow config` invocation is the preflight; the
platform no longer invents engine-specific dry-run flags. Users cannot add
tokens or override paths, profile, executor, reports, plugins, or containers.
The runtime asset doctor and the process runner are bound to the same canonical
Docker executable and local socket identity; a mismatched server composition
fails before launch. Source and runtime assets are operator-owned; submitted
input and reference paths must be absolute, server-local, regular resources.
`process.stageInMode='copy'` prevents tasks from modifying source inputs, at the
cost of additional local disk and staging time.

Prebuilt STAR and Salmon indexes require adapter-private exact manifests, not
only directory digests. Each manifest binds the FASTA/GTF identities, complete
file tree, pinned producer process, tool version, immutable runtime image, and
the relevant 3.26.0 build semantics. STAR records its SA-index strategy and
splice-junction defaults; Salmon records gentrome/k-mer/GENCODE semantics.
Both are checked against the run's index-shaping Advanced parameters
(`skip_gtf_filter`, `skip_gtf_transcript_filter`, and, for Salmon,
`gffread_transcript_fasta`). A self-declared index from another producer or a
manifest incompatible with the normalized run fails before command creation.

Cache coordinates bind the complete runtime build, normalized inputs, verified
FASTQ/reference/rRNA closure, and Nextflow version. Cache ownership is therefore
server-side and incompatible identities cannot share it. PR #151 deliberately
sets `resume=false` and never emits `-resume`: safe retry/resume requires a
durable attempt/session lifecycle and is a later decision, not a caller option.

Hashing proves the resources observed during planning, but Nextflow opens their
paths later. Deployment must therefore keep runtime assets, FASTQ, references,
and database bundles read-only and stable for the run. Atomic content-store
staging that removes this remaining path-replacement window is outside PR #151;
the execution boundary fails closed on races it can observe and records the
complete closure identity used for the run.

## Lifecycle and logging

No RNA-seq-specific lifecycle is introduced. SQLite remains canonical for run
state and durable build identity; Redis/RQ remains the queue and worker boundary.
The existing worker claim, callback, cancellation acknowledgement, and terminal
state rules continue to apply.

Nextflow is executed through the existing `ProcessRunner` and RQ horse-process
group boundary. A temporary Linux subreaper tracks descendants by PID/start-time
identity, including detached and late-forked children. Cancellation, timeout,
callback failure, and exceptional exits apply bounded TERM/KILL cleanup. The
platform then repeatedly stops, kills, removes, and observes only containers
carrying the workspace-derived scope label on the same local daemon; cleanup
uses a server-owned bounded quiescence contract. By default it requires two
continuous seconds with no labelled container, resets that interval after
activity, caps observation at twenty seconds, and reserves five seconds for a
final cleanup and audit. During ordinary execution a cleanup or audit failure
maps to the existing `FAILED` terminal status as an infrastructure failure; it
does not add a workflow-specific terminal status. During cancellation and RQ
stop callbacks the same failure never acknowledges cancellation, leaving the
canonical state diagnosable for operational reconciliation.
Docker exposes no causal fence for an already in-flight create request; a
daemon operation that materializes only after the complete bounded budget is
an unavoidable external-runtime residual rather than a claim of unbounded
cleanup.

Stdout/stderr are bounded and redacted before durable callbacks. The adapter
also declares the Nextflow preflight and execution logs as managed logs; the
platform reads them descriptor-relatively with no-follow/type/size/race checks,
applies the same size- and candidate-work-bounded literal redaction as captured
stdout/stderr, rewrites the safe bytes, and persists only the sanitized content.
Workspace, runtime, Docker, reference, and submitted input paths are private
redaction values and are not copied into public issues.

## rRNA-removal boundary

An rRNA manifest is not a trusted list of arbitrary paths. PR #151 reads it
descriptor-relatively with `O_NOFOLLOW`, applies byte/entry/file bounds, rejects
URLs, absolute or traversing entries, duplicates, symlinks, FIFO/devices, and
replacement races, and hashes the manifest plus every referenced regular file.
Path components are restricted to a conservative ASCII shell-safe set because
the fixed upstream modules interpolate those paths into task scripts. Effective
staged names, names after `.gz` decompression, and Bowtie2 logical basenames
must also be unique, so separate manifest directories cannot create a late
Nextflow staging or header-prefix collision.
The resulting database closure participates in input/reference, execution, and
cache identity evidence.

For a prebuilt SortMeRNA index, an exact manifest fixes SortMeRNA 4.3.7, every
index file, and the database-closure identity. Without a prebuilt index, PR #151
permits only one verified database file and composes one fixed upstream route:
the duplicate preparation index is disabled and the read-filtering subworkflow
builds exactly one index with fixed resources, retry policy, arguments, and the
verified SortMeRNA image. The strategy token and database closure enter build,
execution, and cache identity. This establishes deterministic composition, not
scientific acceptance; a real pinned-container gate remains required before
product exposure. Multi-file no-index construction continues to fail closed.

For Bowtie2 paired-end input, the fixed upstream route retains pairs only when
both mates are unmapped. It does not publish the final paired filtered FASTQ for
`save_non_ribo_reads`. The product contract therefore permits saved filtered
reads for Bowtie2 single-end runs only and rejects that request for paired or
mixed layouts. RiboDetector and GPU RiboDetector remain outside schema 1.0.0.

## Delivery boundary and consequences

The runtime-composed adapter may truthfully declare workspace and command
capabilities after its assets are verified. The default adapter instance remains
contract-only, and `bulk-rnaseq` is not added to the default registry, API, or
frontend in PR #151. No real STAR+Salmon scientific acceptance is claimed.

Artifact discovery and machine-readable QC extraction remain PR #152; HTML
scraping will not define their core contract. Tiny `rapid_quant` qualification
and the full local STAR+Salmon acceptance gate remain PR #153; `rapid_quant` is
not the default product analysis. Default registry and product UI exposure
remain PR #154.

Official upstream coordinates:
[nf-core/rnaseq 3.26.0](https://github.com/nf-core/rnaseq/releases/tag/3.26.0),
[immutable source](https://github.com/nf-core/rnaseq/tree/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4),
[parameter schema](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/nextflow_schema.json),
[license](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/LICENSE), and
[Nextflow 25.04.3](https://github.com/nextflow-io/nextflow/releases/tag/v25.04.3).
