"""Original rules/dispatch with fake tools: scheduling and argv, not scientific bigWigs."""

import base64
import csv
import json
import os
from pathlib import Path
import re
import sys

import pytest

from _complexity import REPO, run


# Coordinate-sorted, samtools-validated two-record BAM/BAI fixtures. SE has
# independent 50M reads at 101/301 (flags 0/16). PE is a complete 99/147 pair
# at 101/201, TLEN +/-150, 50M, with matching mates and 50-base SEQ/QUAL.
# Generated with: samtools view --no-PG -b -o input.bam input.sam; samtools index.
_BAMS = {
    "SE": [
        "H4sIBAAAAAAA/wYAQkMCAFsAc3L0ZTRgYGBw8HDhDPOzMtQz4wz2t0rOzy9KycxLLEnlcggO5Az2s0rOKDLk9AEqMAADLkagHlYgBokzODjxMwAArakSwUkAAAAfiwgEAAAAAAD/BgBCQwIASQArYoCAFCBmtvEUYgTSRkD8HwpAckWGDArMDAyCuIAGyaAIaqsOI8xWAUxbjahtKwAgbCZr7AAAAB+LCAQAAAAAAP8GAEJDAgAbAAMAAAAAAAAAAAA=",
        "QkFJAQEAAAACAAAASRIAAAEAAAAAAFwAAAAAAAAApgAAAAAASpIAAAIAAAAAAFwAAAAAAAAApgAAAAAAAgAAAAAAAAAAAAAAAAAAAAEAAAAAAFwAAAAAAAAAAAAAAAAA",
    ],
    "PE": [
        "H4sIBAAAAAAA/wYAQkMCAFsAc3L0ZTRgYGBw8HDhDPOzMtQz4wz2t0rOzy9KycxLLEnlcggO5Az2s0rOKDLk9AEqMAADLkagHlYgBokzODjxMwAArakSwUkAAAAfiwgEAAAAAAD/BgBCQwIAUQArYYCAFCBmtfEUYmRIZjCCip0A4mlAXJCYWcSgwMzAIIgLaJAMSpDsgNg7GW4vyC1Z////p4W9ANvbSwHwAAAAH4sIBAAAAAAA/wYAQkMCABsAAwAAAAAAAAAAAA==",
        "QkFJAQEAAAACAAAASRIAAAEAAAAAAFwAAAAAAAAArgAAAAAASpIAAAIAAAAAAFwAAAAAAAAArgAAAAAAAgAAAAAAAAAAAAAAAAAAAAEAAAAAAFwAAAAAAAAAAAAAAAAA",
    ],
}


@pytest.fixture
def project(tmp_path, monkeypatch, run_validator):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    # Preserve both original producer and consumer rules, including params
    # dispatch. Only executable boundaries are replaced. The .bw is argv JSON.
    tool = (
        f"#!{sys.executable}\n"
        + """import json, sys
from pathlib import Path
name = Path(sys.argv[0]).name
args = sys.argv[1:]
with Path("tool-calls.jsonl").open("a") as handle:
    handle.write(json.dumps({"tool": name, "argv": args}) + "\\n")
def value(flag):
    return args[args.index(flag) + 1]
if name == "macs3":
    assert args[0] == "callpeak"
    output = Path(value("--outdir"))
    output.mkdir(parents=True, exist_ok=True)
    suffix = "broadPeak" if "--broad" in args else "narrowPeak"
    row = "chr1\\t100\\t200\\tpeak\\t100\\t.\\t5\\t3\\t2"
    if suffix == "narrowPeak":
        row += "\\t50"
    (output / (value("-n") + "_peaks." + suffix)).write_text(row + "\\n")
    prediction = json.loads(Path("model.json").read_text())["prediction"]
    if "--nomodel" not in args and prediction is not None:
        print("predicted fragment length is", prediction)
    else:
        print("No predicted fragment length")
else:
    assert name == "bamCoverage"
    output = Path(value("-o"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(args))
"""
    )
    for name in ("macs3", "bamCoverage"):
        path = bindir / name
        path.write_text(tool)
        path.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))

    def prepare(
        *,
        assay="cuttag",
        layout="SE",
        mode="broad",
        ext="auto",
        role="treatment",
        prediction=150,
    ):
        config_dir = tmp_path / "config"
        config_dir.mkdir()
        reads = []
        for mate in (1, 2):
            path = tmp_path / f"R{mate}.fastq"
            path.write_text(f"@pair/{mate}\n" + "A" * 50 + "\n+\n" + "I" * 50 + "\n")
            reads.append(str(path))
        sample = dict(
            sample="S1",
            fastq_1=reads[0],
            fastq_2=reads[1] if layout == "PE" else "",
            layout=layout,
            assay=assay,
            target="H3K27me3",
            peak_mode=mode,
            genome="tiny",
            bowtie2_index=str(tmp_path / "index"),
            role=role,
            control_sample="",
        )
        samples = [sample]
        if role == "control":
            samples.append(
                dict(sample, sample="T1", role="treatment", control_sample="S1")
            )
        with (config_dir / "samples.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=sample, delimiter="\t")
            writer.writeheader()
            writer.writerows(samples)
        config = dict(
            samples=str(config_dir / "samples.tsv"),
            outdir="results",
            threads=1,
            use_control=role == "control",
            multiqc=False,
            replicate_analysis=False,
            qc={"signal_tracks": False},
            genome_resources={"tiny": {"effective_genome_size": 1000000}},
        )
        if ext is not None:
            config["extend_reads"] = ext
        config_path = config_dir / "config.yaml"
        config_path.write_text(json.dumps(config))
        validation = run_validator(config_path)
        assert validation.returncode == 0, validation.stderr
        bam = tmp_path / "results/S1/02_align/S1.final.bam"
        bam.parent.mkdir(parents=True)
        for suffix, encoded in zip(("", ".bai"), _BAMS[layout]):
            Path(str(bam) + suffix).write_bytes(base64.b64decode(encoded))
        (tmp_path / "model.json").write_text(json.dumps({"prediction": prediction}))
        return config_path

    return prepare


BW = "results/S1/03_bigwig/S1.CPM.bw"
PEAKS = "results/S1/04_peaks/S1"
LOG = "results/S1/logs/S1.macs3.log"


def invoke(workdir, executable, config, *, peaks=False, prioritize=BW, dry_run=False):
    argv = [
        executable,
        "-s",
        REPO / "workflow/Snakefile",
        "--workflow-profile",
        "none",
        "--directory",
        workdir,
        "--configfile",
        config,
        "--cores",
        "1",
        "--scheduler",
        "greedy",
        "--printshellcmds",
    ]
    if dry_run:
        argv.append("--dry-run")
    # Keep normal metadata: the incremental tests must exercise params tracking.
    argv.extend(
        [
            "--allowed-rules",
            "bamcoverage",
            "macs3_callpeak",
            "--prioritize",
            prioritize,
            "--",
            BW,
            *([PEAKS] if peaks else []),
        ]
    )
    return run(workdir, argv)


def calls(workdir):
    return [
        json.loads(line)
        for line in (workdir / "tool-calls.jsonl").read_text().splitlines()
    ]


def bam_inputs(result):
    match = re.search(
        r"(?:local)?rule bamcoverage:\n((?:[ \t]+[^\n]*\n)+)",
        result.stdout + result.stderr,
    )
    assert match, result.stdout + result.stderr
    return re.search(r"(?m)^    input: (.+)$", match[1])[1].split(", ")


def assert_coverage(workdir, extension):
    argv = json.loads((workdir / BW).read_text())
    assert argv == [
        "-b",
        "results/S1/02_align/S1.final.bam",
        "-o",
        BW,
        "--normalizeUsing",
        "CPM",
        "--binSize",
        "10",
        *extension,
        "--numberOfProcessors",
        "1",
    ]


@pytest.mark.parametrize("assay", ["cuttag", "chipseq"])
@pytest.mark.parametrize("ext", ["auto", "yes"])
@pytest.mark.parametrize(
    "prioritize", [BW, PEAKS], ids=["consumer-first", "producer-first"]
)
def test_first_run_waits_for_model(
    tmp_path, project, snakemake_executable, assay, ext, prioritize
):
    config = project(assay=assay, ext=ext)
    assert not (tmp_path / LOG).exists()
    result = invoke(
        tmp_path, snakemake_executable, config, peaks=True, prioritize=prioritize
    )
    assert_coverage(tmp_path, ["--extendReads", "150"])
    assert PEAKS in bam_inputs(result)
    assert [c["tool"] for c in calls(tmp_path)] == ["macs3", "bamCoverage"]


@pytest.mark.parametrize("ext", [None, "auto", "yes"], ids=["default", "auto", "yes"])
def test_explicit_bigwig_pulls_model_producer(
    tmp_path, project, snakemake_executable, ext
):
    config = project(ext=ext)
    result = invoke(tmp_path, snakemake_executable, config)
    assert_coverage(tmp_path, ["--extendReads", "150"])
    assert PEAKS in bam_inputs(result)
    assert [c["tool"] for c in calls(tmp_path)] == ["macs3", "bamCoverage"]


@pytest.mark.parametrize(
    "options,extension",
    [
        ({"ext": "no"}, []),
        ({"ext": 175}, ["--extendReads", "175"]),
        ({"layout": "PE", "ext": "auto"}, []),
        ({"layout": "PE", "ext": "yes"}, ["--extendReads"]),
        ({"role": "control", "ext": "auto"}, ["--extendReads", "200"]),
        ({"role": "control", "ext": "yes"}, ["--extendReads", "200"]),
        ({"mode": "narrow", "ext": "auto"}, ["--extendReads", "200"]),
        ({"mode": "narrow", "ext": "yes"}, ["--extendReads", "200"]),
        ({"assay": "atac", "mode": "narrow"}, ["--extendReads", "200"]),
        ({"assay": "mnase", "layout": "PE", "mode": "nucleosome"}, []),
    ],
)
def test_non_model_paths_do_not_gain_peak_dependency(
    tmp_path, project, snakemake_executable, options, extension
):
    config = project(**options)
    result = invoke(tmp_path, snakemake_executable, config)
    assert_coverage(tmp_path, extension)
    assert PEAKS not in bam_inputs(result)
    assert [c["tool"] for c in calls(tmp_path)] == ["bamCoverage"]
    assert not (tmp_path / LOG).exists()
    if options.get("mode") == "narrow" and options.get("assay", "cuttag") == "cuttag":
        # Explicit peaks still retain the original --nomodel/shift/extsize policy.
        invoke(tmp_path, snakemake_executable, config, peaks=True, prioritize=PEAKS)
        macs = [c["argv"] for c in calls(tmp_path) if c["tool"] == "macs3"]
        assert len(macs) == 1
        argv = macs[0]
        assert "--nomodel" in argv and argv[argv.index("--shift") + 1] == "-100"
        assert argv[argv.index("--extsize") + 1] == "200"
        assert "predicted fragment length is" not in (tmp_path / LOG).read_text()
        assert_coverage(tmp_path, extension)


@pytest.mark.parametrize("prediction", [None, 40], ids=["no-prediction", "too-small"])
def test_completed_model_without_usable_prediction_keeps_fallback(
    tmp_path, project, snakemake_executable, prediction
):
    config = project(prediction=prediction)
    result = invoke(tmp_path, snakemake_executable, config)
    assert [c["tool"] for c in calls(tmp_path)] == ["macs3", "bamCoverage"]
    assert PEAKS in bam_inputs(result)
    assert_coverage(tmp_path, ["--extendReads", "200"])
    assert (
        "MACS3 fragment size unavailable for S1, using --extendReads 200 (fallback)."
        in result.stderr
    )


def test_changed_prediction_triggers_params_then_converges(
    tmp_path, project, snakemake_executable
):
    config = project()
    # Works even on the old graph: producer priority gives a recorded 150.
    invoke(tmp_path, snakemake_executable, config, peaks=True, prioritize=PEAKS)
    assert_coverage(tmp_path, ["--extendReads", "150"])
    initial = calls(tmp_path)
    (tmp_path / LOG).write_text("predicted fragment length is 250\n")
    planned = invoke(tmp_path, snakemake_executable, config, dry_run=True)
    assert "params have changed" in (planned.stdout + planned.stderr).lower()
    assert calls(tmp_path) == initial
    invoke(tmp_path, snakemake_executable, config)
    assert_coverage(tmp_path, ["--extendReads", "250"])
    assert calls(tmp_path) == initial + [
        {"tool": "bamCoverage", "argv": json.loads((tmp_path / BW).read_text())}
    ]
    final = calls(tmp_path)
    unchanged = invoke(tmp_path, snakemake_executable, config)
    assert "Nothing to be done" in unchanged.stdout + unchanged.stderr
    assert calls(tmp_path) == final
