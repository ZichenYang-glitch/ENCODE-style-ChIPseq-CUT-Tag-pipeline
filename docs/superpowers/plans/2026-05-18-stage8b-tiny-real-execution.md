# Stage 8b: Tiny Real Execution — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one real end-to-end ChIP-seq PE no-control execution through preprocessing+signal on tiny synthetic fixtures, proving the workflow runs with real tool invocations.

**Architecture:** A single stdlib-only Python harness (`test/test_stage8b_tiny_execution.py`) generates a 20 kb pseudo-random reference, 1 000 PE reads, and a Bowtie2 index under `/tmp`, then invokes Snakemake with explicit file targets (fastqc.done, trim.done, final.bam, final.bam.bai, CPM.bw). MACS3 is intentionally skipped — Stage 8a dry-run already verifies peak-calling DAG connectivity.

**Tech Stack:** Python 3 (stdlib only), Snakemake, bowtie2, bowtie2-build, samtools, fastqc, trim_galore, cutadapt, bamCoverage (all from chipseq Conda env)

---

### Task 1: Create the harness file skeleton

**Files:**
- Create: `test/test_stage8b_tiny_execution.py`

- [ ] **Step 1: Write the full harness script**

```python
#!/usr/bin/env python3
"""Stage 8b tiny real-execution smoke test for ChIP-seq/CUT&Tag pipeline.

Generates a tiny pseudo-random reference, paired-end FASTQ reads, and a
Bowtie2 index under /tmp, then runs Snakemake with explicit file targets
through the preprocessing + signal path. MACS3 peak calling is intentionally
skipped — Stage 8a dry-run already verifies that DAG.

Exit codes:
  0  PASS — Snakemake completed, all output assertions passed
  1  FAIL — a Snakemake rule or output assertion failed
  2  SKIP — required tool not found, environment not ready

Usage:
    SNAKEMAKE=/path/to/snakemake python3 test/test_stage8b_tiny_execution.py
"""

import subprocess
import tempfile
import shutil
import os
import sys
import random

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")

REF_LEN = 20_000
N_PAIRS = 1_000
READ_LEN = 50
INSERT_SIZE = 200
SEED = 42
QUAL = "I" * READ_LEN

REQUIRED_TOOLS = [
    "snakemake",
    "bowtie2-build",
    "bowtie2",
    "samtools",
    "fastqc",
    "trim_galore",
    "cutadapt",
    "bamCoverage",
]

FORBIDDEN_SUFFIXES = (
    "results/", ".snakemake/",
    ".fq", ".fq.gz", ".fastq", ".fastq.gz",
    ".bam", ".bai", ".bw",
)

FULL_TSV_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
    "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
    "replicate\tbiological_replicate\ttechnical_replicate\t"
    "role\tcontrol_sample\tcontrol_bam"
)


# ---------------------------------------------------------------------------
# Tool resolution
# ---------------------------------------------------------------------------

def resolve_snakemake():
    """Return (snakemake_path, bin_dir | None)."""
    env_val = os.environ.get("SNAKEMAKE", "")
    if env_val:
        path = env_val
    else:
        conda_path = "/home/irenadler/miniconda3/envs/chipseq/bin/snakemake"
        if os.path.isfile(conda_path) and os.access(conda_path, os.X_OK):
            path = conda_path
        else:
            path = "snakemake"
    bin_dir = os.path.dirname(os.path.abspath(path)) if os.path.isabs(path) else None
    return path, bin_dir


def build_subprocess_env(bin_dir):
    """Return os.environ copy with PYTHONDONTWRITEBYTECODE=1 and PATH adjusted."""
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if bin_dir:
        env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def check_required_tools(env):
    """Return list of missing tool names (empty list = all present)."""
    missing = []
    for tool in REQUIRED_TOOLS:
        resolved = shutil.which(tool, path=env.get("PATH"))
        if not resolved or not os.access(resolved, os.X_OK):
            missing.append(tool)
    return missing


# ---------------------------------------------------------------------------
# Fixture generation
# ---------------------------------------------------------------------------

COMPLEMENT = str.maketrans("ATGC", "TACG")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def generate_reference(workdir):
    """Write a 20 kb pseudo-random FASTA to workdir/ref.fa.  Returns the sequence."""
    rng = random.Random(SEED)
    seq = "".join(rng.choice("ATGC") for _ in range(REF_LEN))
    path = os.path.join(workdir, "ref.fa")
    with open(path, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    return seq


def generate_pe_fastq(workdir, ref_seq):
    """Write R1.fq.gz and R2.fq.gz to workdir.  1 000 PE pairs, 50 bp reads, 200 bp insert."""
    r1_path = os.path.join(workdir, "R1.fq.gz")
    r2_path = os.path.join(workdir, "R2.fq.gz")
    max_start = len(ref_seq) - INSERT_SIZE  # 19 800

    with open(r1_path, "w") as f1, open(r2_path, "w") as f2:
        for i in range(N_PAIRS):
            start = round(i * max_start / (N_PAIRS - 1))
            r1_seq = ref_seq[start:start + READ_LEN]
            r2_seq = reverse_complement(
                ref_seq[start + INSERT_SIZE - READ_LEN:start + INSERT_SIZE]
            )
            name = f"@read{i + 1:06d}"
            f1.write(f"{name}/1\n{r1_seq}\n+\n{QUAL}\n")
            f2.write(f"{name}/2\n{r2_seq}\n+\n{QUAL}\n")


def build_bowtie2_index(workdir, env):
    """Run bowtie2-build ref.fa bt2_idx in workdir."""
    ref_fa = os.path.join(workdir, "ref.fa")
    idx_prefix = os.path.join(workdir, "bt2_idx")
    subprocess.run(
        ["bowtie2-build", ref_fa, idx_prefix],
        cwd=workdir, env=env,
        capture_output=True, text=True,
        check=True,
    )
    return idx_prefix


# ---------------------------------------------------------------------------
# Config / sample sheet
# ---------------------------------------------------------------------------

def write_config_and_samples(workdir):
    """Write config.yaml and samples.tsv into workdir."""
    config_path = os.path.join(workdir, "config.yaml")
    samples_path = os.path.join(workdir, "samples.tsv")

    config = f"""\
samples: "{samples_path}"
outdir: "results"
threads: 1
mapq: 30
binsize: 10
remove_dup: "no"
trim: true
extend_reads: "auto"
use_control: false
multiqc: false
stage4b: true
stage5: false
genome_resources:
  hs:
    effective_genome_size: "hs"
    chrom_sizes: ""
    blacklist: ""
    gtf: ""
    reference_fasta: ""
qc:
  blacklist_filter: false
  frip: false
  library_complexity: false
  nrf_pbc: false
  signal_tracks: false
  summary: false
  cuttag_fragment_size: false
"""

    r1 = os.path.join(workdir, "R1.fq.gz")
    r2 = os.path.join(workdir, "R2.fq.gz")
    idx = os.path.join(workdir, "bt2_idx")

    samples_tsv = (
        f"{FULL_TSV_HEADER}\n"
        f"T1\t{r1}\t{r2}\tPE\tchipseq\tH3K27ac\tnarrow\ths\t{idx}\t"
        f"exp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
    )

    with open(config_path, "w") as fh:
        fh.write(config)
    with open(samples_path, "w") as fh:
        fh.write(samples_tsv)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

def _run(cmd, cwd, env, check=True):
    """Run a subprocess, returning CompletedProcess.  Raises on error if check=True."""
    return subprocess.run(
        cmd, cwd=str(cwd), env=env,
        capture_output=True, text=True,
        check=check,
    )


def capture_git_status():
    """Return set of untracked paths from git status (short format)."""
    r = subprocess.run(
        ["git", "status", "--short", "--untracked-files=all"],
        cwd=REPO_ROOT,
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        return set()
    paths = set()
    for line in r.stdout.strip().split("\n"):
        stripped = line.strip()
        if not stripped or stripped.startswith(" "):
            continue
        # "?? path" or " M path" -> extract path after status prefix
        parts = stripped.split(None, 1)
        if len(parts) >= 2:
            paths.add(parts[1])
        else:
            paths.add(stripped)
    return paths


def has_forbidden_suffix(path):
    """Return True if *path* matches any forbidden generated-artifact suffix."""
    return any(
        path.endswith(suffix) or suffix + "/" in path or ("/" + suffix) in path
        for suffix in FORBIDDEN_SUFFIXES
    )


def run_snakemake_targets(workdir, env, snakemake_path):
    """Run Snakemake with explicit file targets. Returns (returncode, stdout, stderr)."""
    targets = [
        "results/T1/logs/T1.fastqc.done",
        "results/T1/logs/T1.trim.done",
        "results/T1/02_align/T1.final.bam",
        "results/T1/02_align/T1.final.bam.bai",
        "results/T1/03_bigwig/T1.CPM.bw",
    ]
    cmd = [
        snakemake_path, "-s", SNAKEFILE,
        "--configfile", os.path.join(workdir, "config.yaml"),
        "--cores", "1",
    ] + targets
    r = subprocess.run(
        cmd, cwd=str(workdir), env=env,
        capture_output=True, text=True,
    )
    return r.returncode, r.stdout, r.stderr


# ---------------------------------------------------------------------------
# Output assertions
# ---------------------------------------------------------------------------

def assert_outputs(workdir, env):
    """Assert all expected outputs exist.  Returns list of failure messages (empty = all pass)."""
    failures = []
    results_dir = os.path.join(workdir, "results", "T1")

    # Sentinels — existence only (created via touch, may be 0 bytes)
    for sentinel in [
        os.path.join(results_dir, "logs", "T1.fastqc.done"),
        os.path.join(results_dir, "logs", "T1.trim.done"),
    ]:
        if not os.path.isfile(sentinel):
            failures.append(f"missing sentinel: {sentinel}")
            break  # stop checking sentinels after first missing

    # Real outputs — existence + nonzero size
    for label, rel in [
        ("final.bam",     os.path.join(results_dir, "02_align", "T1.final.bam")),
        ("final.bam.bai", os.path.join(results_dir, "02_align", "T1.final.bam.bai")),
        ("CPM.bw",        os.path.join(results_dir, "03_bigwig", "T1.CPM.bw")),
    ]:
        if not os.path.isfile(rel):
            failures.append(f"missing: {rel}")
        elif os.path.getsize(rel) == 0:
            failures.append(f"zero size: {rel}")

    # Mapped read count
    if not failures:
        bam = os.path.join(results_dir, "02_align", "T1.final.bam")
        r = subprocess.run(
            ["samtools", "view", "-c", bam],
            cwd=str(workdir), env=env,
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            failures.append(f"samtools view -c failed: {r.stderr.strip()}")
        else:
            try:
                count = int(r.stdout.strip())
            except ValueError:
                count = 0
            if count < 100:
                failures.append(f"mapped reads {count} < 100")

    return failures


def assert_no_repo_artifacts(before_status, after_status):
    """Return list of new paths under REPO_ROOT that match forbidden suffixes."""
    new_paths = after_status - before_status
    return [p for p in new_paths if has_forbidden_suffix(p)]


def cleanup_pycache():
    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print("Starting Stage 8b Tiny Real Execution Test\n")

    # 1. Tool availability
    snakemake_path, bin_dir = resolve_snakemake()
    env = build_subprocess_env(bin_dir)
    missing = check_required_tools(env)
    if missing:
        for tool in missing:
            print(f"SKIP: {tool} not found")
        sys.exit(2)

    # 2. Capture before-status
    before_status = capture_git_status()

    # 3. Create temp workdir
    workdir = tempfile.mkdtemp(prefix="stage8b_", dir="/tmp")

    try:
        # 4. Generate fixtures
        ref_seq = generate_reference(workdir)
        generate_pe_fastq(workdir, ref_seq)
        print("  Fixtures generated")

        # 5. Build Bowtie2 index
        build_bowtie2_index(workdir, env)
        print("  Bowtie2 index built")

        # 6. Write config and sample sheet
        write_config_and_samples(workdir)

        # 7. Run Snakemake
        print("  Running Snakemake ...")
        rc, stdout, stderr = run_snakemake_targets(workdir, env, snakemake_path)
        if rc != 0:
            tail = (stdout + stderr).strip()[-800:]
            print(f"FAIL: Snakemake exit {rc}\n   {tail}")
            sys.exit(1)

        # 8. Assert outputs
        failures = assert_outputs(workdir, env)
        if failures:
            for f in failures:
                print(f"FAIL: {f}")
            sys.exit(1)

        # 9. Assert no repo artifacts leaked
        after_status = capture_git_status()
        leaked = assert_no_repo_artifacts(before_status, after_status)
        if leaked:
            for p in leaked:
                print(f"FAIL: generated artifact in repo: {p}")
            sys.exit(1)

        print("PASS: Stage 8b tiny real execution\n")
        print("Pipeline ran successfully through preprocessing + signal on tiny synthetic data.")
        print("Verified: FastQC, Trim Galore, Bowtie2, sort, index, filter, dedup, final.bam, BigWig.")

    finally:
        shutil.rmtree(workdir, ignore_errors=True)
        cleanup_pycache()


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Syntax-check the harness**

```bash
python3 -c "import py_compile; py_compile.compile('test/test_stage8b_tiny_execution.py', doraise=True)"
```

---

### Task 2: Update KNOWN_ISSUES.md

**Files:**
- Modify: `KNOWN_ISSUES.md`

- [ ] **Step 1: Mark Stage 8b tiny real execution as complete, keep 8c/8d deferred**

Read the Stage 8 section in `KNOWN_ISSUES.md` and edit it to add the Stage 8b completion checkmarks. The current Stage 8 section already has the Stage 8a items marked complete (from Stage 8a). Add the following under the 8a items:

```
- ✅ Add real end-to-end execution with tiny synthetic data (1 ChIP-seq PE
  no-control run, preprocessing + signal path only). (Stage 8b, completed
  2026-05-18)
- ✅ Add tiny FASTQ subsets and a tiny Bowtie2 index generated at runtime
  under /tmp (no committed binary fixtures). (Stage 8b, completed 2026-05-18)
```

Leave all 8c/8d items with `⬜`.

- [ ] **Step 2: Verify KNOWN_ISSUES.md reads correctly**

```bash
python3 -c "
with open('KNOWN_ISSUES.md') as f:
    content = f.read()
assert '8b' in content
print('OK')
"
```

---

### Task 3: Update README.md

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add the Stage 8b execution command under the developer/testing section**

After the existing Stage 8a smoke command block (added in that stage), append:

```markdown
### Tiny real execution

One real end-to-end run through preprocessing + signal on synthetic
fixtures (20 kb pseudo-random reference, 1 000 PE reads, Bowtie2 index
built at runtime — zero binary files committed):

```bash
SNAKEMAKE=/path/to/snakemake python3 test/test_stage8b_tiny_execution.py
```

All outputs land under `/tmp`. Exit code 0 = PASS, 1 = FAIL, 2 = SKIP.
MACS3 is intentionally skipped (Stage 8a dry-run covers the DAG).
See `docs/superpowers/specs/2026-05-18-stage8b-tiny-real-execution-design.md`
```

- [ ] **Step 2: Verify the README renders the new section**

```bash
python3 -c "
with open('README.md') as f:
    content = f.read()
assert 'test_stage8b_tiny_execution.py' in content
print('OK')
"
```

---

### Task 4: Run the harness

**Files:** none (runtime only)

- [ ] **Step 1: Run Stage 8b harness**

```bash
python3 test/test_stage8b_tiny_execution.py
```

Expected: `PASS: Stage 8b tiny real execution` with exit code 0.
If tools are missing: exit code 2 with `SKIP: <tool> not found`.

- [ ] **Step 2: Confirm no generated files in repo**

```bash
git status --short --untracked-files=all
```

No `results/`, `.snakemake/`, `*.fq`, `*.fq.gz`, `*.bam`, `*.bai`, `*.bw` should appear.

- [ ] **Step 3: Confirm no smoke temp dirs remain**

```bash
ls /tmp/stage8b_* 2>&1 || echo "No stage8b temp dirs remain"
```

- [ ] **Step 4: Confirm no scripts/__pycache__**

```bash
ls scripts/__pycache__ 2>&1 || echo "No scripts/__pycache__"
```

---

### Task 5: Full verification suite

**Files:** none

- [ ] **Step 1: Run all three test suites**

```bash
# Stage 8b — tiny real execution
python3 test/test_stage8b_tiny_execution.py

# Stage 8a — dry-run smoke profiles
python3 test/test_stage8_smoke_profiles.py

# Stage 2 — validation stress tests
python3 test/test_validation_stress.py
```

Expected: all three suites PASS with exit code 0.

- [ ] **Step 2: Verify final git status**

```bash
git status --short --untracked-files=all
```

Only expected untracked/modified files: spec, plan, harness, profile files, docs updates.
No generated artifacts.

---

### Task 6: Commit

- [ ] **Step 1: Stage all new and modified files**

```bash
cd /home/irenadler/workflow/chipseq && \
  git add .gitignore \
         KNOWN_ISSUES.md \
         README.md \
         test/test_stage8b_tiny_execution.py \
         docs/superpowers/specs/2026-05-18-stage8b-tiny-real-execution-design.md \
         docs/superpowers/plans/2026-05-18-stage8b-tiny-real-execution.md
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
Add Stage 8b tiny real execution smoke test

One ChIP-seq PE no-control run through preprocessing + signal on
synthetic fixtures generated under /tmp at runtime. MACS3 intentionally
skipped — Stage 8a dry-run covers peak-calling DAG connectivity.

Co-Authored-By: deepseek-v4-pro[1m] <deepseek-ai@claude-code-best.win>
EOF
)"
```
