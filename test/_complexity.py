"""Isolated inputs shared by complexity DAG and real-tool tests."""

import csv
import json
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def project(
    workdir, layout, *, remove_dup="yes", peak_mode="narrow", mapq=17, preseq=True
):
    config_dir = workdir / "config"
    config_dir.mkdir()
    reads = [workdir / f"reads_{mate}.fastq" for mate in (1, 2)]
    for mate, path in enumerate(reads, 1):
        path.write_text(f"@read/{mate}\nACGTACGT\n+\nIIIIIIII\n")
    row = dict(
        sample="S1",
        fastq_1=str(reads[0]),
        fastq_2=str(reads[1]) if layout == "PE" else "",
        layout=layout,
        assay="chipseq",
        target="H3K27ac",
        peak_mode=peak_mode,
        genome="tiny",
        bowtie2_index=str(workdir / "index"),
        experiment="EXP1",
        biological_replicate=1,
        role="treatment",
    )
    with (config_dir / "samples.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    qc = dict(
        blacklist_filter=False,
        frip=False,
        library_complexity=False,
        signal_tracks=False,
        summary=False,
        cuttag_fragment_size=False,
    )
    if preseq is not None:
        qc["preseq_complexity"] = preseq
    config = dict(
        samples=str(config_dir / "samples.tsv"),
        outdir="results",
        threads=1,
        mapq=mapq,
        remove_dup=remove_dup,
        use_control=False,
        multiqc=False,
        replicate_analysis=False,
        qc=qc,
        genome_resources={"tiny": {"effective_genome_size": 1000000}},
    )
    config_path = config_dir / "config.yaml"
    config_path.write_text(json.dumps(config))
    return config_path


def run(workdir, argv, *, check=True):
    """Keep each real invocation and its raw output in its isolated workdir."""
    argv = [str(value) for value in argv]
    path = workdir / "commands.jsonl"
    index = len(path.read_text().splitlines()) if path.exists() else 0
    result = subprocess.run(
        argv, cwd=workdir, capture_output=True, text=True, check=False
    )
    stdout = f"command-{index}.stdout"
    stderr = f"command-{index}.stderr"
    (workdir / stdout).write_text(result.stdout)
    (workdir / stderr).write_text(result.stderr)
    with path.open("a") as handle:
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
    if check:
        assert result.returncode == 0, result.stdout + result.stderr
    return result


def snakemake(
    workdir, executable, config, targets=(), *, dry_run=True, restricted=True
):
    command = [
        executable,
        "-s",
        REPO / "workflow/Snakefile",
        "--workflow-profile",
        "none",
        "--directory",
        workdir,
        "--cores",
        "1",
        "--configfile",
        config,
        "--printshellcmds",
    ]
    if dry_run:
        command.append("--dry-run")
    if restricted:
        command.extend(
            [
                "--allowed-rules",
                "samtools_filter",
                "samtools_index_filt",
                "duplicate_handling",
                "nrf_pbc",
                "preseq_complexity",
            ]
        )
    command.extend(["--", *targets])
    return run(workdir, command, check=False)
