# Bulk RNA-seq Adapter Selection

Status: Accepted for phased implementation; PR #150 remains contract-only.

## Decision

HelixWeave will add `bulk-rnaseq` as a Python `WorkflowAdapter` around a
locally available, immutable copy of nf-core/rnaseq 3.26.0. The upstream
release tag is a lightweight tag at commit
`e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`. PR #150 ships only the adapter's
authoring and validation contract plus the exact upstream parameter and sample
schemas needed to prove that identity. It does not make the adapter runnable or
register it as a default product workflow.

The fixed upstream files are recorded in
`src/encode_pipeline/contracts/nfcore_rnaseq/provenance.json`. In particular,
the parameter schema SHA-256 is
`8f2f84a25c0aec65a18234cf01acdd74f2385e8dfac8417e4bad23a70bfb4388`.
The schema's own `$id` points at the movable `master` branch and is therefore
not an identity coordinate. nf-core/rnaseq is distributed under the MIT
license; the exact license bytes are retained beside the schemas.

This controlled black-box design is preferred to a Snakemake rewrite because
nf-core/rnaseq already carries a mature scientific workflow, documented
parameters, tests, and release practice. Reimplementing that behavior would
create avoidable scientific equivalence and upgrade risk. A thin wrapper that
downloads a tag at run time is also rejected: a movable network dependency
cannot satisfy offline execution, reproducibility, source review, or durable
build identity. Runtime source, Nextflow, plugins, containers, and references
will each receive immutable identities before execution is enabled.

The pipeline and execution engine remain upstream DSL2/Nextflow. The
HelixWeave adapter remains Python and translates platform contracts; no custom
Groovy is required. Workflow-neutral platform, API, persistence, worker, and
frontend code must not depend on nf-core parameter names, Nextflow profiles,
process names, output paths, or validation payloads.

## Parameter ownership

The versioned adapter schema has three closed layers:

- **Standard** exposes stable HelixWeave semantics. Schema 1.0.0 fixes the
  default analysis to STAR alignment plus Salmon quantification and provides
  typed trimming, UMI, ribosomal RNA removal, explicit reference and index
  identity, QC, output, and per-row strandedness controls. The initial
  contract intentionally does not advertise HISAT2 as an
  expression-quantification route.
- **Advanced** is a small allowlist of native scientific parameters from the
  exact 3.26.0 schema. Each value is checked against its unmodified upstream
  property schema and stricter adapter semantics. Unknown parameters and
  Standard aliases fail closed; schema-version upgrades require an explicit
  review of the allowlist and output-impact classification.
- **Platform-owned** includes samplesheet and output paths, launch/work
  directories, executor/profile/resume/cache, trace/report locations,
  Nextflow/Groovy configuration, plugins, Tower/Wave and networking, remote
  configuration, notification, and hardware/container selection. Standard
  resource objects may name immutable reference assets, but callers cannot
  directly set the corresponding native path parameters. Callers also cannot
  provide CLI token lists, shell strings, `extra_*_args`, or raw configuration.

Schema 1.0.0 models ribosomal RNA removal as a closed Standard object with an
`enabled` switch, a `sortmerna` or `bowtie2` tool choice, a
`save_filtered_reads` selection, and adapter-owned manifest/index resource
objects carrying SHA-256 identities. The adapter alone maps these semantics to
`remove_ribo_rna`, `ribo_removal_tool`, `save_non_ribo_reads`,
`ribo_database_manifest`, and `sortmerna_index`. Disabled removal accepts no
tool or resource fields. Both supported tools require an explicit database
manifest, and a prebuilt index is accepted only for SortMeRNA. This prevents
the upstream default manifest, which contains a remote database URL, from
silently becoming a runtime dependency. A manifest digest identifies only the
manifest bytes; PR #151 must reject URL entries inside it, verify the local
database assets it names, and prove that any SortMeRNA index identity is bound
to that exact database reference set. Database licensing remains the operator's
responsibility and is not established by a digest.

The fixed upstream Bowtie2 paired-end route does not publish the filtered
FASTQ generated after its BAM conversion, so schema 1.0.0 permits
`save_filtered_reads` for Bowtie2 only when every submitted sample is
single-end. SortMeRNA supports that selection for both layouts. Supplying a
SortMeRNA index is route/namespace affecting, not content-only: in 3.26.0 the
no-index path wires index building in both genome preparation and the filtering
subworkflow, while prebuilt directories and archives take reuse or unpack
routes. PR #151 must qualify this behavior and either require a verified
prebuilt index or prove that the no-index route is safe and deterministic.

The adapter requires explicit Illumina SE/PE layout and strandedness; it never
infers scientific metadata from file names. Repeated rows with one sample ID
represent upstream-compatible technical libraries or lanes and must agree on
layout and strandedness. Biological replicates use different sample IDs.
Validation is lexical and semantic only and performs no filesystem probe.
The controlled MVP deliberately accepts only absolute, whitespace-free
`.fq.gz` or `.fastq.gz` read paths, a stricter subset of the upstream schema;
uncompressed FASTQ can be added only through a later schema-version decision.
A supplied Salmon index is consumed only for `auto` strandedness inference in
the fixed STAR+Salmon route; otherwise it is rejected rather than silently
ignored. UMI R2 fields likewise require an entirely paired-end submission.

Every accepted native target and every adapter-fixed route invariant is
classified as content-only, sample-set, artifact-set, route/namespace,
additive, subtractive, or filename/format affecting. That map is a required
input to the later artifact contract; accepting or emitting an output-shaping
parameter without recording its effect is not allowed.

## Scientific scope

The MVP ends at trimming, STAR alignment, Salmon quantification, gene and
transcript expression matrices, alignment/QC outputs, FastQC, and MultiQC.
Optional reference-based rRNA removal is supported through SortMeRNA or
Bowtie2 before alignment and quantification.
DESeq2 transformations may be used for exploratory QC/PCA only. Differential
expression significance, UMI variants outside the typed contract, fusion,
single-cell, HPC/cloud, and upload systems are out of scope.

nf-core/rnaseq 3.26.0 marks RiboDetector as experimental and warns that known
ONNX multiprocessing behavior can hang in containerized runs. It also makes
`use_gpu_ribodetector` a hardware/container concern. Neither RiboDetector nor
its GPU switch is part of the formal 1.0.0 Standard or Advanced surface. Any
future exposure requires a deliberate schema revision and real containerized
CPU execution gate; GPU exposure additionally requires a real NVIDIA runtime
gate covering cancellation, timeout, determinism, and resource reporting.

The upstream `rapid_quant` profile skips alignment and merged quantification
and therefore changes product semantics. It may be qualified later as a tiny,
explicit CI gate, but it is not the Standard default or an interchangeable
production analysis profile.

## Phased delivery boundary

- **PR #150:** this decision, pinned source/schema/license evidence, adapter
  skeleton, schema 1.0.0, strict validation, and conformance tests only.
- **PR #151:** immutable offline runtime assets, Nextflow version and plugin
  policy, workspace/command construction, logs, cancellation, timeout, resume,
  cache, source/build identity, and verified local rRNA database/index
  materialization. It must not fetch a pipeline or reference at run time.
- **PR #152:** deterministic artifact inventory and machine-readable QC
  extraction; HTML scraping cannot define a core result contract.
- **PR #153:** explicit tiny `rapid_quant` CI gate with its reduced semantics,
  followed by a local STAR+Salmon acceptance gate using controlled data.
- **PR #154:** default registry, API/frontend authoring and results exposure,
  documentation, and release qualification only after the preceding gates are
  green.

Official upstream evidence: the
[3.26.0 release](https://github.com/nf-core/rnaseq/releases/tag/3.26.0),
[immutable parameter schema](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/nextflow_schema.json),
[immutable sample schema](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/assets/schema_input.json),
[versioned usage contract](https://nf-co.re/rnaseq/3.26.0/docs/usage/),
[license](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/LICENSE), and
[`rapid_quant` profile](https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/conf/rapid_quant.config).
