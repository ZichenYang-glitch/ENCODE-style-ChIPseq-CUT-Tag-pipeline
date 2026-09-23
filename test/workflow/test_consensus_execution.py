"""Original consensus rule execution on tiny precomputed narrowPeak inputs."""

import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest


@pytest.fixture
def consensus_project(tmp_path, tmp_config, monkeypatch):
    """Use the current test interpreter for the original stdlib-only script."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "python3").symlink_to(sys.executable)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / ".cache"))

    def prepare(assay, *, outdir="results", parent=True, replicate_analysis=True):
        read1, read2 = (tmp_path / f"reads_{mate}.fq" for mate in (1, 2))
        for mate, path in enumerate((read1, read2), 1):
            path.write_text(f"@r/{mate}\nACGTACGT\n+\nIIIIIIII\n")
        columns = (
            "sample fastq_1 fastq_2 layout assay target peak_mode genome "
            "bowtie2_index experiment biological_replicate role"
        ).split()
        rows = ["\t".join(columns)]
        for biorep in (1, 2):
            rows.append(
                "\t".join(
                    map(
                        str,
                        (
                            f"S{biorep}",
                            read1,
                            read2,
                            "PE",
                            assay,
                            "CTCF",
                            "narrow",
                            "tiny",
                            tmp_path / "index",
                            "EXP",
                            biorep,
                            "treatment",
                        ),
                    )
                )
            )
        config = {
            "outdir": outdir,
            "threads": 1,
            "use_control": False,
            "trim": False,
            "multiqc": False,
            "replicate_analysis": replicate_analysis,
            "chipseq_idr": False,
            "qc": dict.fromkeys(
                (
                    "blacklist_filter",
                    "frip",
                    "library_complexity",
                    "nrf_pbc",
                    "signal_tracks",
                    "summary",
                    "cuttag_fragment_size",
                ),
                False,
            ),
            "genome_resources": {"tiny": {"effective_genome_size": 1000000}},
        }
        if parent is not None:
            # Leave the nested consensus switch omitted to exercise its default.
            config["reproducibility"] = {"enabled": parent}
        _, config_path, _ = tmp_config(config=config, samples="\n".join(rows) + "\n")
        base = Path(outdir) / "experiments/EXP/06_reproducibility"
        name = f"EXP.{assay}.macs3.narrow"
        inputs = []
        for biorep, start, signal, summit, lone in (
            (1, 100, 5, 20, 500),
            (2, 120, 10, 30, 800),
        ):
            path = (
                tmp_path
                / base
                / "consensus/biorep_peaks"
                / (f"EXP.biorep{biorep}.{assay}.macs3.narrow_peaks.narrowPeak")
            )
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                f"chr1\t{start}\t{start + 100}\tshared\t100\t.\t{signal}\t3\t2\t{summit}\n"
                f"chr1\t{lone}\t{lone + 100}\talone\t100\t.\t2\t3\t2\t20\n"
            )
            inputs.append(path)
        return {
            "config": config_path,
            "inputs": inputs,
            "peak": base / "consensus" / f"{name}.consensus.narrowPeak",
            "summary": base / "consensus" / f"{name}.consensus.summary.tsv",
            "final": base
            / "final"
            / f"{name}.replicate_validated.consensus.narrowPeak",
        }

    return prepare


def invoke(workdir, executable, snakefile, project, *, default=False, dry_run=False):
    argv = [
        executable,
        "-s",
        snakefile,
        "--workflow-profile",
        "none",
        "--configfile",
        project["config"],
        "--directory",
        str(workdir),
        "--cores",
        "1",
        "--printshellcmds",
    ]
    if dry_run:
        argv.append("--dry-run")
    if not default:
        # No conda activation: the consumer needs only Python stdlib. Avoid
        # environment metadata hashing; this does not test conda provisioning.
        argv.extend(["--drop-metadata", "--allowed-rules", "consensus_compute_narrow"])
        argv.extend(["--", str(project["peak"])])
    result = subprocess.run(argv, cwd=workdir, capture_output=True, text=True)
    record = workdir / "commands.jsonl"
    index = len(record.read_text().splitlines()) if record.exists() else 0
    stdout, stderr = f"command-{index}.stdout", f"command-{index}.stderr"
    (workdir / stdout).write_text(result.stdout)
    (workdir / stderr).write_text(result.stderr)
    with record.open("a") as handle:
        handle.write(
            json.dumps(
                dict(
                    argv=argv,
                    cwd=str(workdir),
                    exit_code=result.returncode,
                    stdout=stdout,
                    stderr=stderr,
                )
            )
            + "\n"
        )
    return result


@pytest.mark.parametrize(
    "assay,outdir",
    [
        ("chipseq", "results"),
        ("atac", "results"),
        ("cuttag", "results"),
        ("cuttag", "results with spaces"),
    ],
)
def test_original_narrow_rule_preserves_peak_and_summary(
    tmp_path,
    consensus_project,
    snakemake_executable,
    snakefile,
    run_validator,
    assay,
    outdir,
):
    project = consensus_project(assay, outdir=outdir)
    validation = run_validator(project["config"])
    assert validation.returncode == 0, validation.stderr
    dag = invoke(
        tmp_path, snakemake_executable, snakefile, project, default=True, dry_run=True
    )
    assert dag.returncode == 0, dag.stdout + dag.stderr
    assert "rule consensus_compute_narrow:" in dag.stdout + dag.stderr
    result = invoke(tmp_path, snakemake_executable, snakefile, project)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / project["peak"]).read_text() == (
        "chr1\t100\t220\tconsensus_peak_1\t1000\t.\t10\t-1\t-1\t50\n"
    )
    with (tmp_path / project["summary"]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    final = str(project["final"]) if assay == "cuttag" else ""
    assert len(rows) == 1
    assert rows[0] == {
        "experiment": "EXP",
        "assay": assay,
        "peak_mode": "narrow",
        "caller": "macs3",
        "n_bioreps": "2",
        "min_replicates": "2",
        "reciprocal_overlap": "0.5",
        "consensus_peak_count": "1",
        "support_distribution": '{"2": 1}',
        "biorep_labels": "1,2",
        "source_peak_files": json.dumps([str(p) for p in project["inputs"]]),
        "final_method": "consensus" if assay == "cuttag" else "none",
        "final_output": final,
    }
    if " " in outdir:
        assert " " in rows[0]["final_output"]
    # The final path is summary metadata, not another output of this consumer.
    assert not (tmp_path / project["final"]).exists()


@pytest.mark.parametrize("assay", ["chipseq", "atac"])
@pytest.mark.parametrize("parent", [None, False])
def test_default_parent_gate_does_not_schedule_consensus(
    tmp_path,
    consensus_project,
    snakemake_executable,
    snakefile,
    assay,
    parent,
):
    project = consensus_project(assay, parent=parent)
    result = invoke(
        tmp_path, snakemake_executable, snakefile, project, default=True, dry_run=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (
        re.search(r"(?m)^(?:local)?rule consensus_", result.stdout + result.stderr)
        is None
    )
    # Default target gating is not a ban on explicitly requesting this output.
    explicit = invoke(tmp_path, snakemake_executable, snakefile, project, dry_run=True)
    assert explicit.returncode == 0, explicit.stdout + explicit.stderr
    assert "rule consensus_compute_narrow:" in explicit.stdout + explicit.stderr


def test_replicate_analysis_gate_omits_consensus(
    tmp_path,
    consensus_project,
    snakemake_executable,
    snakefile,
):
    project = consensus_project("chipseq", replicate_analysis=False)
    result = invoke(
        tmp_path, snakemake_executable, snakefile, project, default=True, dry_run=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (
        re.search(r"(?m)^(?:local)?rule consensus_", result.stdout + result.stderr)
        is None
    )


@pytest.mark.real_execution
def test_original_narrow_rule_rejects_malformed_peaks(
    tmp_path,
    consensus_project,
    snakemake_executable,
    snakefile,
):
    # Snakemake's failed-job diagnostics query conda even without --use-conda.
    # Keep this dependency in the existing scientific execution CI tier.
    assert shutil.which("conda"), "Failed-rule diagnostics require conda on PATH"
    project = consensus_project("cuttag")
    project["inputs"][0].write_text("chr1\t100\t200\n")
    result = invoke(tmp_path, snakemake_executable, snakefile, project)
    assert result.returncode != 0
    assert "expected 10 narrowPeak columns, got 3" in result.stdout + result.stderr
    assert not (tmp_path / project["peak"]).exists()
    assert not (tmp_path / project["summary"]).exists()
