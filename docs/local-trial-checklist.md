# HelixWeave v0.3.0 Local Trial Checklist

Use this checklist to evaluate the exact v0.3.0 candidate on one workstation.
It covers the supported local/small trusted-team product boundary, both bundled
workflows, the browser journey, persistence, and cleanup. It is not a production
deployment or scientific benchmarking guide.

## Trial asset and scope

- [ ] Obtain the wheel, sdist, source-trial archive, and `SHA256SUMS` from the
      same candidate commit.
- [ ] Run `sha256sum --check SHA256SUMS` before extracting or installing.
- [ ] Confirm the source archive expands under one
      `helixweave-v0.3.0/` directory and does not require a Git checkout.
- [ ] Confirm the candidate identity is HelixWeave v0.3.0 while the Python
      distribution/import/CLI compatibility identities remain
      `encode-pipeline`, `encode_pipeline`, and `encode-*`.
- [ ] Confirm no OCI image, JDK, reference, index, FASTQ, result database,
      environment, cache, or credential is bundled.

The wheel and sdist are the compatibility Python distribution. The source-trial
archive is the complete no-Git-checkout product tree containing the frontend,
workflows, scripts, documentation, and lock files. v0.3.0 does not publish a
container image.

Supported product scope:

- file-backed SQLite as canonical lifecycle and result-metadata state;
- Redis/RQ as the local queue and worker boundary;
- local filesystem workspaces;
- ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase; and
- `bulk-rnaseq`, fixed to nf-core/rnaseq 3.26.0.

Authentication, multi-user isolation, HPC/cloud schedulers, object storage,
Tower/Wave, production support SLAs, and cross-attempt resume are not provided.

## Compatibility distribution checks

The wheel and sdist checks are deliberately separate from the source-trial
browser journey. From a directory outside the extracted archive and any Git
checkout, create two new virtual environments and install one verified
compatibility artifact into each:

```bash
TRIAL_ROOT=$(mktemp -d)
ASSET_ROOT=/absolute/path/to/verified/assets

python3.12 -m venv "$TRIAL_ROOT/wheel"
env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  "$TRIAL_ROOT/wheel/bin/python" -m pip install \
  "$ASSET_ROOT/encode_pipeline-0.3.0-py3-none-any.whl[api]"
env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  "$TRIAL_ROOT/wheel/bin/python" -m pip check

python3.12 -m venv "$TRIAL_ROOT/sdist"
env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  "$TRIAL_ROOT/sdist/bin/python" -m pip install \
  "$ASSET_ROOT/encode_pipeline-0.3.0.tar.gz[api]"
env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  "$TRIAL_ROOT/sdist/bin/python" -m pip check
```

Set `ASSET_ROOT` to the flat directory whose `SHA256SUMS` was verified. Run all
checks below through each environment's Python from the same outside
directory. Do not use an editable install, `PYTHONPATH`, user site packages, or
files from the source-trial tree for these two distribution checks.

- [ ] `python -c 'import importlib.metadata; print(importlib.metadata.version("encode-pipeline"))'`
      prints `0.3.0`.
- [ ] `python -c 'import encode_pipeline; print(encode_pipeline.__version__)'`
      prints `0.3.0`.
- [ ] `encode-validate --help`, `encode-manifest --help`,
      `encode-dag --help`, and `encode-worker --help` start from the installed
      distribution.
- [ ] The installed API can export OpenAPI and apply its packaged Alembic
      migrations without reading the source tree.

## Source-trial environment

The complete browser/scientific product is intentionally exercised from the
separately checksummed source-trial archive, not from the compatibility wheel.
Inside the extracted archive, create its locked environment and install its
exact source package without dependency resolution:

```bash
unset PYTHONPATH
export PYTHONNOUSERSITE=1
export HOME="$PWD/.local/trial-home"
export MAMBA_ROOT_PREFIX="$PWD/.local/micromamba"
export PLAYWRIGHT_BROWSERS_PATH="$PWD/.local/playwright-browsers"
mkdir -p "$HOME"

micromamba create -p .local/envs/ci-fast --file workflow/envs/ci-fast.lock
./.local/envs/ci-fast/bin/python -m pip install --no-index --no-deps \
  --no-build-isolation -e ".[api]"
./.local/envs/ci-fast/bin/python -m pip check
npm --prefix frontend ci
npm --prefix frontend exec -- playwright install chromium
export PATH="$PWD/.local/envs/ci-fast/bin:$PATH"
```

The launcher intentionally runs the exact extracted `src/`, workflow, and
frontend tree. This proves the no-Git-checkout source-trial asset; it is not
additional evidence that the wheel contains the complete product. A clean
Linux host may also need Playwright's separately approved
`install-deps chromium` host operation.

## Doctor and availability

Run the maintained doctor before opening ports or creating runtime data:

```bash
python scripts/run_local_platform.py --doctor
```

- [ ] Required Python, Redis, Snakemake, Node.js, npm, frontend, and platform
      checks are successful.
- [ ] The output lists exactly
      `encode-style-chipseq-cuttag-atac-mnase` and `bulk-rnaseq`.
- [ ] ENCODE authoring/execution status is explicit.
- [ ] With no Bulk runtime binding, Bulk reports
      `authoring=available`, `execution=not_configured`, and
      `WORKFLOW_EXECUTION_NOT_CONFIGURED`.
- [ ] Missing optional Bulk execution assets do not make the core platform
      doctor fail.
- [ ] A partial, malformed, or drifted Bulk binding reports `unavailable`
      rather than runnable.
- [ ] Output contains stable reason codes and actionable documentation, but no
      private runtime paths, environment values, Docker socket details,
      credentials, or exception text.

`not_configured` means the operator has not supplied a complete server-owned
execution binding. `unavailable` means supplied execution coordinates could not
be safely admitted. Neither state disables Bulk input authoring or validation;
both must keep create/start fail-closed.

## Start the local product

Choose a task-owned runtime directory and start the deterministic trial:

```bash
python scripts/run_local_platform.py \
  --input-authoring-demo \
  --runtime-root "$PWD/.local/v0.3.0-trial"
```

- [ ] The launcher reports one browser URL and the generated demo-input path.
- [ ] Redis, FastAPI, one RQ worker, and the frontend become ready.
- [ ] `/api/v1/workflows/` returns both workflows.
- [ ] Pressing Ctrl-C once stops every process started by the launcher.

The launcher may reuse an explicitly configured compatible Redis server; it
does not flush or stop a server it did not start.

## ENCODE workflow journey

Use the generated JSON config and TSV path printed by
`--input-authoring-demo`.

- [ ] Open **Workflows** and select the ENCODE-style workflow.
- [ ] Open **Author inputs / New analysis**.
- [ ] Import the generated TSV and paste `resultsConfig` in YAML mode.
- [ ] Review and validate the inputs.
- [ ] Confirm validation creates an immutable snapshot.
- [ ] Create the run, wait for successful preflight, and select **Start run**.
- [ ] Confirm SQLite/RQ/worker lifecycle reaches `SUCCEEDED`.
- [ ] Confirm **Runs** and run detail show stable activity and persisted logs.
- [ ] Confirm **Artifacts** and **QC** show the deterministic result manifest,
      QC summary, and indexed metrics.
- [ ] Open an artifact detail and download it using the displayed artifact
      generation and revision; confirm the downloaded bytes match the indexed
      source.

This deterministic tiny route proves authoring, validation, snapshot,
submission, execution, lifecycle, artifact/QC, and revision-bound download
contracts. It does not prove biological validity or production-scale
performance.

## Bulk RNA-seq unavailable journey

Perform this section without configuring a Bulk runtime.

- [ ] Select **Bulk RNA-seq** from the same workflow list.
- [ ] Confirm its detail page identifies nf-core/rnaseq 3.26.0 and shows
      execution as not configured.
- [ ] Exercise schema-driven Form and YAML authoring for representative SE, PE,
      and repeated-lane samples.
- [ ] Review strandedness, UMI, trimming, and ribosomal RNA removal fields;
      confirm executor, profile, work directory, container, plugin, network,
      raw CLI, and Groovy config are not user-editable.
- [ ] Validate safe inputs and confirm structured validation remains available.
- [ ] Confirm the submission area clearly says execution is unavailable.
- [ ] Confirm both API create and start boundaries fail closed with stable,
      redacted public errors; a disabled browser button alone is not evidence.
- [ ] Confirm no Nextflow process, container, network pull, or scientific
      workspace is started.

The runtime-available STAR+Salmon/SortMeRNA, SE/PE/repeated-lane, artifact/QC,
offline, lifecycle, timeout, and cancellation evidence is the protected exact
execution Gate recorded for the reviewed closure. Do not recreate a rootful
runner merely for this ordinary local trial. A changed execution closure
requires a new protected Gate under the release checklist.

## Browser and responsive review

At each viewport below, inspect the workflow list, both workflow details,
authoring, Review, run history, run detail, Artifacts, QC, and the
revision-bound download action:

- [ ] 1440×900
- [ ] 1024×768
- [ ] 390×844
- [ ] 360×800

For every viewport:

- [ ] no page-level horizontal overflow;
- [ ] no overlapping or clipped controls;
- [ ] navigation, tabs, forms, tables/cards, dialogs, and download actions
      remain keyboard- and pointer-operable;
- [ ] unavailable execution messaging remains visible without hiding
      authoring/validation; and
- [ ] no ENCODE-specific wording is hard-coded into generic Runs, Artifacts,
      or QC surfaces.

Save screenshots with the candidate commit and viewport in their evidence
metadata. Screenshots must not expose local paths, environment values, private
or non-synthetic sample identifiers, or credentials.

## Optional Omics Intake Bundle boundary

Omics Intake is an optional external producer, not a bundled third workflow or
runtime dependency.

- [ ] HelixWeave only consumes a reviewed, versioned Bundle through its
      read-only import boundary.
- [ ] The import does not acquire data, mutate the producer project, trust
      arbitrary paths, or bypass adapter validation.
- [ ] A user can complete both bundled workflow trials without installing or
      running Omics Intake.

See the
[Bundle consumption boundary](architecture/omics-intake-bundle-consumption.md)
for its exact contract.

## Persistence and upgrade

- [ ] A fresh start creates a file-backed SQLite database at the selected
      runtime root.
- [ ] API and worker use the same absolute database/workspace coordinates.
- [ ] Normal API restart preserves runs, events, logs, artifacts, QC, and
      lifecycle state.
- [ ] The maintained `20260714_07` fixture upgrades to `20260717_08` with the
      release test.
- [ ] Before upgrading an existing candidate database, stop services and
      follow the
      [v0.3.0 database upgrade procedure](development/local-platform-runtime.md#v030-database-upgrade).

## Cleanup and evidence

Press Ctrl-C and wait for launcher shutdown before cleanup.

- [ ] No launcher-owned API, frontend, RQ worker/work horse, Redis, Snakemake,
      Nextflow, or container process remains.
- [ ] No task-owned port or Unix socket remains bound.
- [ ] Preserve the SQLite database and screenshots if they are release
      evidence; otherwise remove only the explicitly selected trial runtime.
- [ ] Do not delete an unowned runtime directory, shared Redis data, Docker
      resources, user caches, references, or another checkout.
- [ ] Do not run `docker system prune`.

Record:

- exact candidate commit and tree;
- asset SHA-256 values;
- installed package/API version;
- doctor output with both workflow states;
- database revision before/after any upgrade;
- commands and test outcomes;
- four viewport screenshots;
- inherited or rerun protected Bulk Gate identity; and
- anything not run plus the reason.
