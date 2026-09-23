"""Genuine phantompeakqualtools output and original rules on small synthetic BAMs.

These are unpaired, coordinate-sorted 50M records with correlated strand starts,
not a simulation of alignment or a claim about biological quality.
"""

import csv
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import sys

import pytest

from _complexity import REPO, run


pytestmark = pytest.mark.real_execution


@pytest.fixture
def cc_tools(tmp_path, monkeypatch):
    spp = shutil.which(os.environ.get("RUN_SPP", "run_spp.R"))
    samtools = shutil.which(os.environ.get("SAMTOOLS", "samtools"))
    assert spp and samtools, (
        "Real tests require locked phantompeakqualtools and samtools"
    )
    # run_spp.R uses /usr/bin/env Rscript. Select the same tool environment.
    rbin = Path(spp).resolve().parent
    assert (rbin / "Rscript").is_file()
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "python3").symlink_to(sys.executable)
    monkeypatch.setenv(
        "PATH", os.pathsep.join((str(bindir), str(rbin), os.environ.get("PATH", "")))
    )
    return spp, samtools


@pytest.fixture(params=["single", "multi"])
def cc_input(tmp_path, cc_tools, request):
    spp, samtools = cc_tools
    rng = random.Random(881)
    shifts = [150] if request.param == "single" else [100, 200, 300]
    records = []
    for shift in shifts:
        for i in range(10000):
            pos = rng.randint(1000, 900000)
            records.extend(
                [(pos, 0, f"m{shift}-{i}a"), (pos + shift - 49, 16, f"m{shift}-{i}b")]
            )
    sam = tmp_path / "input.sam"
    with sam.open("w") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n")
        for pos, flag, name in sorted(records):
            assert 1 <= pos <= 1000000 - 49
            handle.write(
                f"{name}\t{flag}\tchr1\t{pos}\t60\t50M\t*\t0\t0\t"
                + "A" * 50
                + "\t"
                + "I" * 50
                + "\n"
            )
    bam = tmp_path / "input.bam"
    run(tmp_path, [samtools, "view", "-b", "--no-PG", "-o", bam, sam])
    run(tmp_path, [samtools, "quickcheck", "-v", bam])
    (tmp_path / "input-identity.json").write_text(
        json.dumps(
            dict(
                seed=881,
                shifts=shifts,
                records=len(records),
                sha256=hashlib.sha256(bam.read_bytes()).hexdigest(),
            )
        )
    )
    qc = tmp_path / "automatic.cc.qc"
    run(
        tmp_path,
        [
            spp,
            f"-c={bam}",
            f"-savp={tmp_path / 'automatic.pdf'}",
            f"-out={qc}",
            "-x=-500:15",
            "-rf",
        ],
    )
    raw = qc.read_text().strip().split("\t")
    candidates = raw[2].split(",")
    assert len(candidates) == 1 if request.param == "single" else len(candidates) > 1
    assert all(float(value) > 0 for value in candidates)
    assert (tmp_path / "automatic.pdf").read_bytes().startswith(b"%PDF-")
    return bam, qc, raw


def assert_summary(path, raw, sample):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == [
            "sample",
            "cc_qc_file",
            "estimated_fragment_length",
            "phantom_peak",
            "nsc",
            "rsc",
            "quality_flag",
        ]
        rows = list(reader)
    nsc, rsc = float(raw[8]), float(raw[9])
    low = [
        name for name, bad in (("low_nsc", nsc < 1.05), ("low_rsc", rsc < 0.8)) if bad
    ]
    assert rows == [
        dict(
            sample=sample,
            cc_qc_file=f"{sample}.cc.qc",
            estimated_fragment_length=str(float(raw[2].split(",")[0])),
            phantom_peak=str(float(raw[4])),
            nsc=str(nsc),
            rsc=str(rsc),
            quality_flag="_".join(low) or "ok",
        )
    ]


def test_real_candidates_and_scalar_cli(tmp_path, cc_tools, cc_input):
    bam, qc, raw = cc_input
    if "," not in raw[2]:
        before = hashlib.sha256(bam.read_bytes()).hexdigest()
        forced = tmp_path / "forced.cc.qc"
        run(
            tmp_path,
            [
                cc_tools[0],
                f"-c={bam}",
                f"-savp={tmp_path / 'forced.pdf'}",
                f"-out={forced}",
                "-x=-500:15",
                "-rf",
                "-speak=0",
            ],
        )
        zero = forced.read_text().strip().split("\t")
        assert float(zero[2]) == 0
        assert zero[:2] == raw[:2]
        assert zero[4:8] == raw[4:8]  # Same curve's phantom peak and background.
        assert hashlib.sha256(bam.read_bytes()).hexdigest() == before
        assert (tmp_path / "forced.pdf").read_bytes().startswith(b"%PDF-")
    output = tmp_path / "cli.tsv"
    original = qc.read_bytes()
    run(
        tmp_path,
        [
            sys.executable,
            REPO / "scripts/parse_cross_correlation.py",
            "--input",
            qc,
            "--output",
            output,
        ],
    )
    assert_summary(output, raw, "automatic")
    assert qc.read_bytes() == original


def test_original_cross_correlation_and_summary(
    tmp_path, cc_input, tmp_config, snakemake_executable
):
    bam, _, automatic = cc_input
    columns = "sample fastq_1 fastq_2 layout assay target peak_mode genome bowtie2_index role".split()
    rows = [
        "\t".join(columns),
        "\t".join(
            map(
                str,
                (
                    "S1",
                    tmp_path / "unused.fq",
                    "",
                    "SE",
                    "chipseq",
                    "CTCF",
                    "narrow",
                    "tiny",
                    tmp_path / "index",
                    "treatment",
                ),
            )
        ),
    ]
    _, config, _ = tmp_config(
        config={
            "outdir": "results",
            "threads": 1,
            "use_control": False,
            "multiqc": False,
            "replicate_analysis": False,
            "qc": {"cross_correlation": True},
            "genome_resources": {"tiny": {"effective_genome_size": 1000000}},
        },
        samples="\n".join(rows) + "\n",
    )
    final = tmp_path / "results/S1/02_align/S1.final.bam"
    final.parent.mkdir(parents=True)
    shutil.copyfile(bam, final)
    summary = Path("results/multiqc/cross_correlation_summary.tsv")
    # Existing external toolchain via PATH, not a conda provisioning test.
    result = run(
        tmp_path,
        [
            snakemake_executable,
            "-s",
            REPO / "workflow/Snakefile",
            "--workflow-profile",
            "none",
            "--configfile",
            config,
            "--directory",
            tmp_path,
            "--cores",
            "1",
            "--printshellcmds",
            "--drop-metadata",
            "--allowed-rules",
            "cross_correlation",
            "cross_correlation_summary",
            "--",
            summary,
        ],
    )
    output = result.stdout + result.stderr
    assert "-p=1" in output and "-x=-500:15" in output and "-rf" in output
    raw_path = tmp_path / "results/S1/05_qc/cross_correlation/S1.cc.qc"
    raw = raw_path.read_text().strip().split("\t")
    assert (
        raw[1:] == automatic[1:]
    )  # Includes ordered candidates, NSC/RSC, quality tag.
    assert "-speak=" not in output
    assert_summary(tmp_path / summary, raw, "S1")
    assert raw_path.with_name("S1.cc.plot.pdf").read_bytes().startswith(b"%PDF-")
    assert final.read_bytes() == bam.read_bytes()
