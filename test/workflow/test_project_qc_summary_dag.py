"""Project peak QC must not schedule peak callers for MNase treatments."""

import csv
import gzip
import json
import os
import re
import subprocess

import pytest


def _prepare_project(workdir, assays, summary, missing_fastq):
    """Write standalone DAG inputs; no alignment or scientific execution."""
    config_dir = workdir / "config"
    config_dir.mkdir()
    samples = []
    for sample, assay in assays:
        reads = [workdir / f"{sample}.R{mate}.fastq.gz" for mate in (1, 2)]
        if not missing_fastq:
            for mate, path in enumerate(reads, start=1):
                with gzip.open(path, "wt") as handle:
                    handle.write(f"@{sample}/{mate}\n")
                    handle.write("ACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")
        samples.append(
            {
                "sample": sample,
                "fastq_1": str(reads[0]),
                "fastq_2": str(reads[1]),
                "layout": "PE",
                "assay": assay,
                "target": "H3" if assay == "mnase" else "H3K27ac",
                "peak_mode": "nucleosome" if assay == "mnase" else "narrow",
                "genome": "tiny",
                # The DAG uses the prefix as a parameter, not an input file.
                "bowtie2_index": str(workdir / "tiny_index"),
                "experiment": f"EXP_{sample}",
                "biological_replicate": "1",
                "role": "treatment",
            }
        )
    samples_path = config_dir / "samples.tsv"
    with samples_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=samples[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(samples)
    (workdir / "chrom.sizes").write_text("chr1\t1000\n")
    (workdir / "blacklist.bed").write_text("chr1\t900\t910\n")
    config = {
        "samples": str(samples_path),
        "outdir": "results",
        "threads": 1,
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "genome_resources": {
            "tiny": {
                "effective_genome_size": 1000,
                "chrom_sizes": str(workdir / "chrom.sizes"),
                "blacklist": str(workdir / "blacklist.bed"),
            }
        },
    }
    if summary is not None:
        config["qc"] = {"summary": summary}
    # Also satisfy Snakefile's relative default configfile in the isolated cwd.
    path = config_dir / "config.yaml"
    path.write_text(json.dumps(config))
    return path


def _jobs(output, rule):
    return re.findall(rf"(?m)^(?:local)?rule {rule}:\n((?:[ \t]+[^\n]*\n)+)", output)


def _sample_jobs(output, rule):
    return {
        re.search(r"\bwildcards:.*\bsample=([^,\s]+)", job).group(1)
        for job in _jobs(output, rule)
    }


MIXED = (("CS1", "chipseq"), ("MN1", "mnase"))


@pytest.mark.parametrize(
    "assays,summary,missing_fastq",
    [
        pytest.param(MIXED, None, False, id="mixed-default-summary"),
        pytest.param(MIXED, True, False, id="mixed-summary-on"),
        pytest.param(MIXED, False, False, id="mixed-summary-off"),
        pytest.param((("MN1", "mnase"),), True, False, id="mnase-only"),
        pytest.param(
            (("CS1", "chipseq"), ("CS2", "chipseq")),
            True,
            False,
            id="peak-only",
        ),
        pytest.param(MIXED, None, True, id="missing-fastq"),
    ],
)
def test_project_summary_assay_dependencies(
    tmp_path,
    snakefile,
    snakemake_executable,
    run_validator,
    assays,
    summary,
    missing_fastq,
):
    config = _prepare_project(tmp_path, assays, summary, missing_fastq)
    validation = run_validator(config)
    assert validation.returncode == 0, validation.stderr

    env = os.environ.copy()
    env.pop("SNAKEMAKE_PROFILE", None)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    result = subprocess.run(
        [
            snakemake_executable,
            "-s",
            snakefile,
            "--workflow-profile",
            "none",
            "--directory",
            str(tmp_path),
            "--cores",
            "1",
            "--configfile",
            str(config),
            "--dry-run",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    (tmp_path / "dry-run.stdout").write_text(result.stdout)
    (tmp_path / "dry-run.stderr").write_text(result.stderr)
    output = result.stdout + result.stderr
    if missing_fastq:
        assert result.returncode != 0, output
        assert "MissingInputException in rule trim_galore" in output
        assert "does not use MACS3" not in output
        return
    assert result.returncode == 0, output

    peak_samples = {sample for sample, assay in assays if assay != "mnase"}
    mnase_samples = {sample for sample, assay in assays if assay == "mnase"}
    summary_enabled = summary is not False
    summaries = peak_samples if summary_enabled else set()
    assert _sample_jobs(output, "macs3_callpeak") == peak_samples
    assert _sample_jobs(output, "peak_counts") == summaries
    assert _sample_jobs(output, "qc_summary") == summaries
    assert _sample_jobs(output, "frip") == peak_samples
    for rule in (
        "mnase_qc_summary",
        "mnase_split_sub",
        "mnase_split_mono",
        "mnase_split_di",
        "mnase_dyad_bigwig",
        "mnase_mono_bigwig",
    ):
        assert _sample_jobs(output, rule) == mnase_samples

    projects = _jobs(output, "project_qc_summary")
    assert len(projects) == int(bool(summaries))
    if summaries:
        inputs = re.search(r"^    input: (.+)$", projects[0], re.MULTILINE).group(1)
        assert set(inputs.split(", ")) == {
            f"results/{sample}/01_qc/{sample}.qc_summary.tsv" for sample in summaries
        }
    # Both assay families retain their completion and result-manifest targets.
    assert _sample_jobs(output, "pipeline_done") == {sample for sample, _ in assays}
    manifests = _jobs(output, "result_manifest")
    assert len(manifests) == int(summary_enabled)
    if summary_enabled:
        inputs = re.search(r"^    input: (.+)$", manifests[0], re.MULTILINE).group(1)
        assert {
            f"results/{sample}/01_qc/{sample}.mnase_qc_summary.tsv"
            for sample in mnase_samples
        } <= set(inputs.split(", "))
