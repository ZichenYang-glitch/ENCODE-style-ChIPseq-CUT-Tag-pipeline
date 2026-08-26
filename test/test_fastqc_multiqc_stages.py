"""Focused DAG and real-report checks for staged FastQC / MultiQC."""

import gzip
import hashlib
from pathlib import Path
import subprocess

from _tool_resolver import resolve_tool


REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

TSV_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
    "technical_replicate\trole\tcontrol_sample\tcontrol_bam"
)


def _write_config(workdir, samples_path, *, trim):
    config_path = workdir / "config.yaml"
    config_path.write_text(
        f'samples: "{samples_path}"\n'
        'outdir: "results"\n'
        "threads: 1\n"
        "mapq: 30\n"
        "binsize: 10\n"
        'remove_dup: "no"\n'
        f"trim: {str(trim).lower()}\n"
        'extend_reads: "auto"\n'
        "use_control: false\n"
        "multiqc: false\n"
        "stage4b: false\n"
        "stage5: false\n"
        "genome_resources:\n"
        "  hs:\n"
        '    effective_genome_size: "hs"\n'
        '    chrom_sizes: ""\n'
        '    blacklist: ""\n'
        '    gtf: ""\n'
        '    reference_fasta: ""\n'
        "qc:\n"
        "  blacklist_filter: false\n"
        "  frip: false\n"
        "  library_complexity: false\n"
        "  nrf_pbc: false\n"
        "  signal_tracks: false\n"
        "  summary: false\n"
        "  cuttag_fragment_size: false\n",
        encoding="utf-8",
    )
    return config_path


def _write_fastq(path, *, mate):
    path.parent.mkdir(parents=True, exist_ok=True)
    sequence = "ACGT" * 12 + ("AC" if mate == 1 else "TG")
    quality = "I" * len(sequence)
    with gzip.open(path, "wt") as handle:
        for index in range(20):
            handle.write(f"@read{index:03d}/{mate}\n{sequence}\n+\n{quality}\n")


def _fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_size


def _run_snakemake(workdir, config_path, targets, *, dry_run=False, env=None):
    command = [
        SNAKEMAKE,
        "-s",
        str(SNAKEFILE),
        "--configfile",
        str(config_path),
        "--cores",
        "2",
    ]
    if dry_run:
        command.extend(["-n", "-p"])
    command.extend(targets)
    return subprocess.run(
        command,
        cwd=workdir,
        env=env,
        capture_output=True,
        text=True,
    )


def test_se_trim_false_dag_only_requests_raw_r1(tmp_path):
    reads = tmp_path / "inputs" / "reads.fastq.gz"
    _write_fastq(reads, mate=1)
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        TSV_HEADER
        + "\n"
        + f"SE1\t{reads}\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tunused\t"
        "exp1\tcondition\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    config = _write_config(tmp_path, samples, trim=False)

    result = _run_snakemake(
        tmp_path,
        config,
        ["results/SE1/logs/SE1.fastqc.done"],
        dry_run=True,
    )

    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "SE1.raw.R1_fastqc.html" in output
    assert "SE1.raw.R1_fastqc.zip" in output
    assert ".raw.R2_fastqc" not in output
    assert "/fastqc/trimmed/" not in output
    assert "stage=trimmed" not in output
    assert "rule trim_galore:" not in output
