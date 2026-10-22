"""Run-identity guard counterexamples with real downstream tasks (Codex R2/R3).

A valid marker with a missing, corrupted, or mismatched effective-config
record must stop real production rules (here: real FastQC), and a symlinked
.helixweave directory must not receive record writes. These tests use tiny
controlled fixtures only; they are guard-behavior evidence, not scientific
library validation.

Usage:
    SNAKEMAKE=/path/to/snakemake python3 -m pytest -m real_execution \
        test/real_execution/test_run_identity_guard.py -v
"""

import gzip
import os
import shutil
import subprocess

import pytest
from _tool_resolver import require_external_tools, resolve_tool

TEST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TEST_ROOT)
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")

pytestmark = pytest.mark.real_execution

REQUIRED_TOOLS = ["snakemake", "fastqc"]

FULL_TSV_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
    "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
    "replicate\tbiological_replicate\ttechnical_replicate\t"
    "role\tcontrol_sample\tcontrol_bam"
)

FASTQC_HTML = "results/T1/01_qc/fastqc/raw/T1.raw.R1_fastqc.html"
MARKER = "results/.helixweave/analysis_identity.json"
RECORD = "results/.helixweave/effective_config.json"


def _env():
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    missing = [
        tool
        for tool in REQUIRED_TOOLS
        if not shutil.which(resolve_tool(tool, tool.upper()), path=env.get("PATH"))
    ]
    require_external_tools(missing, "run-identity guard real execution")
    return env


def _snakemake(env):
    return resolve_tool("snakemake", "SNAKEMAKE")


def _write_project(workdir):
    """Write a tiny PE project whose only real work is FastQC."""
    r1 = os.path.join(workdir, "R1.fq.gz")
    r2 = os.path.join(workdir, "R2.fq.gz")
    for path, mate in ((r1, 1), (r2, 2)):
        with gzip.open(path, "wt") as handle:
            for i in range(8):
                handle.write(f"@read{i:04d}/{mate}\nACGTACGT\n+\nIIIIIIII\n")
    config_path = os.path.join(workdir, "config.yaml")
    samples_path = os.path.join(workdir, "samples.tsv")
    config = (
        f'samples: "{samples_path}"\n'
        'outdir: "results"\n'
        "threads: 1\n"
        "trim: false\n"
        "use_control: false\n"
        "multiqc: false\n"
        "replicate_analysis: false\n"
        "chipseq_idr: false\n"
        'genome_resources: {hs: {effective_genome_size: "hs"}}\n'
        "qc:\n"
        "  blacklist_filter: false\n"
        "  frip: false\n"
        "  library_complexity: false\n"
        "  nrf_pbc: false\n"
        "  signal_tracks: false\n"
        "  summary: false\n"
        "  cuttag_fragment_size: false\n"
    )
    with open(config_path, "w") as handle:
        handle.write(config)
    samples = (
        f"{FULL_TSV_HEADER}\n"
        f"T1\t{r1}\t{r2}\tPE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
        f"exp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
    )
    with open(samples_path, "w") as handle:
        handle.write(samples)
    return config_path


def _run(env, workdir, config_path, target):
    return subprocess.run(
        [
            _snakemake(env),
            "-s",
            SNAKEFILE,
            "--configfile",
            config_path,
            "--cores",
            "1",
            target,
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=env,
    )


def test_real_fastqc_stops_when_record_missing_or_corrupt(tmp_path):
    env = _env()
    workdir = str(tmp_path)
    config_path = _write_project(workdir)

    # Positive control: a consistent run reaches the real FastQC rule.
    result = _run(env, workdir, config_path, FASTQC_HTML)
    assert result.returncode == 0, result.stdout + result.stderr
    assert os.path.exists(os.path.join(workdir, FASTQC_HTML))
    assert os.path.exists(os.path.join(workdir, MARKER))
    assert os.path.exists(os.path.join(workdir, RECORD))

    # Missing record: the real downstream task must not run.
    os.remove(os.path.join(workdir, RECORD))
    result = _run(env, workdir, config_path, FASTQC_HTML)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "record is missing" in output, output
    assert not os.path.exists(os.path.join(workdir, RECORD))

    # Corrupted record: also refused, and never auto-repaired in place.
    shutil.rmtree(os.path.join(workdir, "results"))
    result = _run(env, workdir, config_path, FASTQC_HTML)
    assert result.returncode == 0, result.stdout + result.stderr
    record_path = os.path.join(workdir, RECORD)
    with open(record_path, "w") as handle:
        handle.write("{invalid\n")
    result = _run(env, workdir, config_path, FASTQC_HTML)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "record is unreadable" in output, output
    with open(record_path) as handle:
        assert handle.read() == "{invalid\n"


def test_marker_rule_refuses_symlinked_metadata_dir(tmp_path):
    env = _env()
    workdir = str(tmp_path)
    config_path = _write_project(workdir)
    escape = tmp_path / "escape"
    escape.mkdir()
    os.makedirs(os.path.join(workdir, "results"))
    os.symlink(str(escape), os.path.join(workdir, "results", ".helixweave"))

    result = _run(env, workdir, config_path, MARKER)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "must not be a symlink" in output, output
    assert list(escape.iterdir()) == []
