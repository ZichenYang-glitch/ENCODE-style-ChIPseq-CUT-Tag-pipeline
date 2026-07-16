# Bulk RNA-seq Offline Execution Boundary

Status: Accepted for phased implementation. PR #151 enables only an explicitly
runtime-composed adapter and PR #152 adds its closed result contract; product
exposure remains deferred.

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
  Its JVM is the operator-staged Amazon Corretto 21.0.7.6.1 Linux x64 release
  (runtime `21.0.7+6-LTS`, GPL-2.0 with Classpath Exception). The release
  matches the Corretto 21 family used by the official Nextflow 25.04.3
  container baseline while turning that moving family into an exact asset.
  The release
  archive is fixed by size and SHA-256
  `8bb627728d147e7507b2e38a5ef872549e895da50c2685d435c0d4c15ba95eb4`;
  the 457-file expanded tree is
  `a910f21c80d8d1fcc24ded6a5759b9f5ddd9c64643f0d683df3ef243a35a8ae8`
  and `bin/java` is
  `d4c2eea1b5fd0e16bcdf080dc2d162dcd58f3d1d87faa17fc66e0fc473272b9c`.
  No JDK binary is committed to Git. Admission verifies the
  archive, closed no-symlink tree, executable, exact `java -version`, exact
  Nextflow version/build, and configuration canary. JDK archive/tree/java
  identities are part of runtime, build, and cache identity.
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
Nextflow distribution, JDK archive/tree/executable, plugin archive/tree,
container inventory and availability lock, and the fixed upstream schemas.
Changing any coordinate produces a different durable build identity.

Heavyweight admission is separated from run-time liveness. A process-local
`RuntimeAssetAdmission` performs the full source, Nextflow, JDK, plugin, OCI
archive/layer, and canary verification once, records a descriptor-relative
metadata witness (`dev`, inode, mode, link count, size, mtime, and ctime), and
reuses one immutable verified object across build capture, planning, command
construction, and doctor calls. A cache hit reads no source/plugin/JDK/archive
payload bytes: it compares bounded metadata before and after an exact Docker
config/rootfs availability check. Any witness or Docker endpoint change drops
the evidence and triggers a complete new admission; a Docker liveness failure
never falls back to stale evidence. PID changes also discard the cache, and a
worker constructs a fresh admission from the raw operator binding rather than
receiving or deserializing a verified object.

The production contract contains 829 source files, 457 JDK files, 106 plugin
tree files, and 34 unique OCI images for 56 supported processes. One admission
therefore has at most 1,499 top-level content passes, including at most 34
whole-archive hashes and 34 archive/config/layer closure parses. Every process
that shares an upstream image coordinate must name the same OCI digest,
canonical archive, and distribution-manifest closure; a divergent operator
lock fails before the second closure is read. Before this split, the four
preflight calls and three worker calls would each repeat those image passes:
136 plus 102 archive hashes and the same number of closure parses. The bounded
API/preflight process and a fresh worker now each perform one admission (34 +
34), while subsequent calls perform only metadata observations and Docker
liveness checks. Tiny archive tests validate behavior; these production counts
are asserted independently from the committed manifests.

## Workspace and command boundary

Workspace planning first verifies FASTQ, FASTA, GTF, optional reference index,
and rRNA resources without modifying them. It then deterministically emits:

- a sorted official five-column nf-core samplesheet
  (`sample,fastq_1,fastq_2,strandedness,seq_platform`) with normalized
  per-row `ILLUMINA`, and canonical params JSON with no conflicting global
  `seq_platform`;
- a platform-owned Nextflow config with the local executor, exact container
  selectors, a default-deny container selector, explicit audited alias-level
  overrides, report paths, and no-network Docker policy;
- execution and cache identity documents; and
- isolated `engine/launch`, `engine/work`, `engine/cache`, `NXF_HOME`, temporary,
  log, report, and result directories.

The command uses fixed argv with `shell=False`: the pinned local Nextflow binary
runs the pinned local source path with the generated config and params file,
the Docker profile, `-offline`, the server-owned work directory, and a bounded
identity-derived run name. Main execution and preflight each have exactly one
hard `-C` configuration entry and no soft `-c`. The generated file first
`includeConfig`s the verified source tree's manifest-bound `nextflow.config`,
then applies all platform overrides. Nextflow therefore ignores
`launchDir/nextflow.config`, `$NXF_HOME/config`, and every other implicit config
source. A real pinned canary places conflicting marker/executor settings in
both implicit locations and proves they are absent while the Docker profile,
nf-schema plugin, container selector, and fixed nf-core manifest still parse.
The launch directory is inside the run workspace. An explicit command-owned
`nextflow config` invocation is the preflight; the platform no longer invents
engine-specific dry-run flags. Users cannot add tokens or override paths,
profile, executor, reports, plugins, or containers.

`JAVA_HOME`, `NXF_JAVA_HOME`, and `JAVA_CMD` point to the verified Corretto
tree; `PATH` contains only its `bin` plus fixed system directories before the
process runner adds its fixed Docker CLI directory. `LD_LIBRARY_PATH` and every
inherited Conda/Mamba coordinate are cleared. Java option variables are not in
the process runner's inheritance allowlist. Runtime/JDK paths are redaction
values. The host Conda environment, `JAVA_HOME`, or Java option variables
cannot select or modify the JVM.
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

## Result boundary

PR #152 adds a separate, explicitly runtime-composed results adapter. It
declares `artifact_extract` and `qc_summary_extract` only in addition to the
verified workspace and command capabilities; the contract-only adapter and the
default registry remain unchanged. Its versioned results contract is fixed to
nf-core/rnaseq 3.26.0 STAR+Salmon and has SHA-256
`2e2a64389b66a6fcfb59d83281562b8973bbf3846542a4f20ddaa3bda3e0fe9f`.

Artifact discovery derives a finite list of exact paths from normalized sample
identities and output-shaping parameters. It does not recursively glob the
results tree. It checks the fixed STAR+Salmon namespaces for foreign sample
identities, reads the two bounded MultiQC sample-threshold tables when present,
and treats every other required-path absence as failure. `save_reference` and
StringTie remain fail-closed until their complete trees have manifests.
Qualimap, dupRadar, DESeq2 QC plots, preseq, generated scripts, opaque R
objects, raw tool logs, and files known to contain commands or private paths
are explicit public exclusions rather than dynamically discovered artifacts.
MultiQC 1.33 HTML is also excluded because its default report embeds analysis
file paths; publication requires a future pinned privacy-override canary.
Nextflow execution reports and logs remain platform-private outside the
results root; moving them would not make their commands, work paths, and engine
state suitable public artifacts.

The workflow-neutral artifact and QC services enforce descriptor-relative
no-follow reads, regular-file and component checks, stable identity reopens, a
1 TiB per-artifact ceiling, and finite candidate/source/metric ceilings derived
to cover the platform's 1,000-row authoring bound. The bulk result contract
records their exact values, including 16 MiB per QC source and 256 MiB total.
Persisted artifact identity, size, produced time, and metric produced time
remain platform concerns. The bulk adapter supplies stable output types, media
types, sample/run scope, assay, and deterministic ordering. The authoring
contract has no experiment identity, so it is not invented in result metadata.

MultiQC 1.33 filename cleaning is not injective over the original 1.0.0 sample
identifier alphabet. The adapter therefore reserves every literal from the
pinned MultiQC defaults and nf-core config that could rewrite an authored
sample ID. Validation rejects those IDs before execution rather than guessing a
cleaned name or permitting two biological samples to collapse into one report
identity. The exact upstream config digests and 159-literal policy are part of
the versioned result contract. In addition, canonical validation builds a
run-global owner map for canonical IDs, the pinned authored-PE `_1`/`_2`
`table_sample_merge` aliases, actual post-UMI downstream aliases, and exact
FASTQ `simpleName` replacements. Conflicting owners fail with one stable,
value-free issue before workspace planning. Ownership includes both canonical
sample and mate role, so a basename cannot silently reassign R1 to R2. Repeated
lanes claim the same owner and are valid. The graph mirrors MultiQC 1.33's
fixed order: prepended nf-core `extra_fn_clean_exts`, default `fn_clean_exts`,
`fn_clean_trim`, then global exact name replacement. It resolves every observed
FASTQ identity through the same global replacement table before checking final
sample-and-mate ownership. It therefore rejects canonical IDs changed by trim
tokens, replacement keys invalidated by cleaning, duplicate PE source keys,
cross-sample replacement/no-op collisions, and the R2 role collapse created
when upstream suppresses both replacements after an R1 no-op. This
route-specific graph is not applied when MultiQC is disabled.
The machine-data extractor repeats the Cutadapt row-owner check so a bypassed
or stale validation snapshot cannot cause dictionary overwrite.

This reservation is a deliberate pre-release correctness amendment to the
authoring schema labeled 1.0.0 in PR #150. The workflow has never been in the
default registry or product API, so no user-visible 1.0.0 input contract has
been released. The amendment is made before exposure; after product release,
an accepted-set change requires a new schema version and compatibility policy.

QC uses a closed source and metric catalog: FastQC 0.12.1 ZIP data, fastp 1.0.1
JSON, sanitized MultiQC 1.33 Cutadapt and Picard tables, STAR final summaries,
Salmon metadata, featureCounts summaries, and selected RSeQC 5.0.4 native
tables. MultiQC HTML is neither published nor parsed. Percentages are converted directly to
`Decimal` fractions, unknown or duplicate semantic coordinates fail, and
missing optional tools never produce synthetic zeroes. Trimmed FastQC `Total
Sequences` provides an exact final retained-read count; paired mates must
agree. MultiQC Cutadapt integer fields provide exact input and
post-adapter/quality base counts before later length/pair filtering. FastQC
abbreviates `Total Bases`, while the pinned Trim Galore module does not publish
its exact JSON and its raw report contains command-line fields, so PR #152 does
not mislabel either source as an exact final retained-base metric. A future
runtime-produced sanitized derivative is required for that value.
RSeQC transcript integrity numbers retain their upstream 0–100 score semantics;
the workflow-neutral QC unit vocabulary therefore adds `score` rather than
mislabeling TIN as a fraction or ratio. The RNA-seq adapter, not the platform,
owns the TIN-specific range check. RSeQC 5.0.4 `tin.py` also applies a
case-sensitive global `basename.replace("bam", "")`, while the pinned nf-core
wrapper expects the unmodified BAM basename. When the effective module set
contains TIN, sample IDs containing lowercase `bam` are therefore rejected
before workspace planning; IDs are never rewritten. CSI indexing removes TIN
upstream and does not activate this restriction.

The RSeQC inner-distance contract follows the fixed module's output labels:
`<sample>.inner_distance.txt` is per-read-pair distance detail and is published
as `bulk_rnaseq.rseqc.inner_distance.distance`; it is not a summary.
`<sample>.inner_distance_mean.txt` contains the mean, median, and standard
deviation summary and remains the `.mean` artifact. Salmon 1.10.3 writes
`num_processed` and `num_mapped` from its observed-fragment and mapped-fragment
counters. Public metric keys therefore use `processed_fragments` and
`mapped_fragments`, and `mapping_fraction` is explicitly fragment-based. One PE
fragment is one read pair; one SE fragment is one single read.

Discovery and indexing each reopen path components and reject observable
replacement, but they are not one atomic directory snapshot. The existing
input-path TOCTOU and the result replacement interval between adapter discovery
and platform persistence remain explicit deployment risks.

Execution composition identity includes the adapter variant: `runtime-v1` for
the execution-only adapter and `results-v1` for the result-capable adapter.
Consequently their build, workspace, cache, and resume identities cannot be
interchanged even when they share the same pinned runtime assets. The controlled
implementation manifest also binds the production SQLite composition,
SQLAlchemy repository/models, Alembic environment, and the complete current
upgrade-to-head revision set. An added unlisted revision or a change to atomic
artifact replacement/QC invalidation makes the old manifest fail closed and
changes the results build identity; binding only the in-memory repository is
insufficient.

## Delivery boundary and consequences

The runtime-composed adapter may truthfully declare workspace and command
capabilities after its assets are verified. The default adapter instance remains
contract-only, and `bulk-rnaseq` is not added to the default registry, API, or
frontend in PR #151. No real STAR+Salmon scientific acceptance is claimed.

PR #152 supplies deterministic artifact discovery and machine-readable QC
extraction without HTML scraping, but does not claim scientific execution
acceptance. Tiny `rapid_quant` qualification and the full local STAR+Salmon
acceptance gate remain PR #153; `rapid_quant` is not the default product
analysis. Default registry and product UI exposure remain PR #154.

Official upstream coordinates:
[nf-core/rnaseq 3.26.0](https://github.com/nf-core/rnaseq/releases/tag/3.26.0),
[immutable source](https://github.com/nf-core/rnaseq/tree/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4),
[parameter schema](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/nextflow_schema.json),
[license](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/LICENSE), and
[Nextflow 25.04.3](https://github.com/nextflow-io/nextflow/releases/tag/v25.04.3).
The pinned JVM is the official
[Amazon Corretto 21.0.7.6.1 release](https://github.com/corretto/corretto-21/releases/tag/21.0.7.6.1),
including its release-bound
[GPLv2 license and Classpath Exception](https://github.com/corretto/corretto-21/blob/21.0.7.6.1/LICENSE).
