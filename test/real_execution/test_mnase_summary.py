"""Original MNase/Picard rules, real samtools counts, and a Picard tool stub.

Prebuilt fragment BAMs and placeholder bigWigs isolate dependency scheduling;
these tests do not validate Picard insert-size distributions or signal tracks.
"""

import csv
import json
import os
import re
import shutil
import sys

import pytest

from _complexity import REPO, run
from _tool_resolver import resolve_tool


pytestmark = pytest.mark.real_execution

SUMMARY = "results/M1/01_qc/M1.mnase_qc_summary.tsv"
METRICS = "results/M1/05_qc/picard/M1.insert_size_metrics"
HEADER = (
    "sample assay peak_mode sub_min sub_max mono_min mono_max di_min di_max "
    "dyad_min dyad_max sub_bam mono_bam di_bam dyad_bigwig mono_bigwig "
    "insert_size_metrics caller_danpos3_enabled caller_inps_enabled "
    "caller_sem_enabled sub_reads mono_reads di_reads"
).split()


def _bam(workdir, samtools, name, pairs, length):
    """Complete 99/147 pairs, matching mates, 50M and TLEN; sorted by samtools."""
    sam = workdir / f"{name}.sam"
    records = []
    for index in range(pairs):
        left = 101 + 1000 * index
        right = left + length - 50
        for flag, pos, mate, tlen in (
            (99, left, right, length),
            (147, right, left, -length),
        ):
            records.append(
                "\t".join(
                    map(
                        str,
                        (
                            f"{name}_{index}",
                            flag,
                            "chr1",
                            pos,
                            60,
                            "50M",
                            "=",
                            mate,
                            tlen,
                            "A" * 50,
                            "I" * 50,
                        ),
                    )
                )
            )
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n" + "\n".join(records) + "\n"
    )
    bam = workdir / f"{name}.bam"
    run(workdir, [samtools, "sort", "--no-PG", "-o", bam, sam])
    run(workdir, [samtools, "index", bam])
    run(workdir, [samtools, "quickcheck", bam])
    assert run(workdir, [samtools, "view", "-c", bam]).stdout.strip() == str(2 * pairs)
    return bam


@pytest.fixture
def project(tmp_path, monkeypatch, run_validator):
    samtools = shutil.which(resolve_tool("samtools", "SAMTOOLS"))
    assert samtools, "This real-execution test requires samtools, with no count stub"
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "samtools").symlink_to(samtools)
    (bindir / "python3").symlink_to(sys.executable)
    picard = bindir / "picard"
    picard.write_text(
        f"#!{sys.executable}\n"
        + """import json, sys
from pathlib import Path
args = sys.argv[1:]
assert args[0] == "CollectMultipleMetrics"
with Path("picard-calls.jsonl").open("a") as handle:
    handle.write(json.dumps(args) + "\\n")
values = dict(arg.split("=", 1) for arg in args[1:])
assert Path(values["I"]).is_file() and Path(values["R"]).is_file()
mode = Path("picard-mode.txt").read_text()
if mode == "fail":
    print("PICARD_STUB_FAILURE", file=sys.stderr)
    sys.exit(23)
for suffix in ("alignment_summary_metrics", "insert_size_metrics",
               "quality_distribution_metrics", "insert_size_histogram.pdf"):
    if mode == "missing" and suffix == "insert_size_metrics":
        continue
    Path(values["O"] + "." + suffix).write_text("PICARD_STUB_OUTPUT\\n")
print("PICARD_STUB_COMPLETE")
"""
    )
    picard.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    (tmp_path / "tool-identity.json").write_text(json.dumps({"samtools": samtools}))

    def prepare(
        *,
        picard=True,
        summary=False,
        old_metrics=False,
        reference=True,
        mode="success",
        mixed=False,
        validate=True,
    ):
        config_dir = tmp_path / "config"
        config_dir.mkdir()
        fastqs = [tmp_path / f"R{mate}.fastq" for mate in (1, 2)]
        for mate, path in enumerate(fastqs, 1):
            path.write_text(f"@pair/{mate}\n" + "A" * 50 + "\n+\n" + "I" * 50 + "\n")
        ids = [("M1", "mnase", "treatment", "MC" if mixed else "")]
        if mixed:
            ids += [
                ("MC", "mnase", "control", ""),
                ("T1", "chipseq", "treatment", "CC"),
                ("CC", "chipseq", "control", ""),
            ]
        samples = [
            dict(
                sample=sid,
                fastq_1=str(fastqs[0]),
                fastq_2=str(fastqs[1]),
                layout="PE",
                assay=assay,
                target="H3",
                genome="tiny",
                peak_mode="nucleosome" if assay == "mnase" else "narrow",
                bowtie2_index=str(tmp_path / "index"),
                experiment=sid,
                biological_replicate=1,
                role=role,
                control_sample=control,
            )
            for sid, assay, role, control in ids
        ]
        with (config_dir / "samples.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=samples[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(samples)
        ref = tmp_path / "reference.fa"
        if reference:
            ref.write_text(">chr1\n" + "A" * 10000 + "\n")
            run(tmp_path, [samtools, "faidx", ref])
        qc = dict(
            blacklist_filter=False,
            frip=False,
            library_complexity=False,
            nrf_pbc=False,
            signal_tracks=False,
            summary=summary,
            cuttag_fragment_size=False,
        )
        if picard is not None:
            qc["picard_metrics"] = picard
        config = dict(
            samples=str(config_dir / "samples.tsv"),
            outdir="results",
            threads=1,
            trim=False,
            use_control=mixed,
            multiqc=False,
            replicate_analysis=False,
            qc=qc,
            genome_resources={
                "tiny": {
                    "effective_genome_size": 10000,
                    "reference_fasta": str(ref) if reference else "",
                }
            },
            mnase={
                "mono_range": [130, 205],
                "fragments": {"sub": [80, 120], "mono": [140, 210], "di": [300, 410]},
                "dyad_range": [135, 195],
                "callers": {"danpos3": False, "inps": False, "sem": False},
            },
        )
        config_path = config_dir / "config.yaml"
        config_path.write_text(json.dumps(config))
        result = run_validator(config_path)
        (tmp_path / "validation.json").write_text(json.dumps(vars(result)))
        if validate:
            assert result.returncode == 0, result.stdout + result.stderr
        for name, pairs, length in (("sub", 1, 100), ("mono", 2, 150), ("di", 3, 350)):
            bam = _bam(tmp_path, samtools, name, pairs, length)
            dest = tmp_path / f"results/M1/03_fragments/M1.{name}.bam"
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(bam, dest)
            shutil.copyfile(str(bam) + ".bai", str(dest) + ".bai")
        final = tmp_path / "results/M1/02_align/M1.final.bam"
        final.parent.mkdir(parents=True)
        shutil.copyfile(tmp_path / "mono.bam", final)
        shutil.copyfile(tmp_path / "mono.bam.bai", str(final) + ".bai")
        for kind in ("mono", "dyad"):
            path = tmp_path / f"results/M1/04_signal/M1.{kind}.CPM.bw"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("Prebuilt path placeholder, not a scientific bigWig\n")
        (tmp_path / "picard-mode.txt").write_text(mode)
        if old_metrics:
            path = tmp_path / METRICS
            path.parent.mkdir(parents=True)
            path.write_text("Existing metrics marker\n")
        return config_path

    return prepare


def invoke(
    workdir,
    executable,
    config,
    *,
    producer=False,
    priority=SUMMARY,
    default_targets=False,
    check=True,
):
    args = [
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
        "--latency-wait",
        "1",
    ]
    if default_targets:
        args.append("--dry-run")
    else:
        args += [
            "--allowed-rules",
            "mnase_qc_summary",
            "picard_collect_multiple_metrics",
            "--prioritize",
            priority,
            "--",
            SUMMARY,
            *([METRICS] if producer else []),
        ]
    return run(workdir, args, check=check)


def assert_summary(workdir, metrics):
    with (workdir / SUMMARY).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == HEADER
    expected = dict(
        sample="M1",
        assay="mnase",
        peak_mode="nucleosome",
        sub_min="80",
        sub_max="120",
        mono_min="140",
        mono_max="210",
        di_min="300",
        di_max="410",
        dyad_min="135",
        dyad_max="195",
        insert_size_metrics=metrics,
        caller_danpos3_enabled="false",
        caller_inps_enabled="false",
        caller_sem_enabled="false",
        sub_reads="2",
        mono_reads="4",
        di_reads="6",
    )
    for kind in ("sub", "mono", "di"):
        expected[f"{kind}_bam"] = f"results/M1/03_fragments/M1.{kind}.bam"
    for kind in ("dyad", "mono"):
        expected[f"{kind}_bigwig"] = f"results/M1/04_signal/M1.{kind}.CPM.bw"
    assert rows == [expected]  # samtools read-record counts, not PE fragments.


def assert_dependency(result):
    output = result.stdout + result.stderr
    block = re.search(
        r"(?:local)?rule mnase_qc_summary:\n((?:[ \t]+[^\n]*\n)+)", output
    )
    assert block, output
    assert METRICS in re.search(r"(?m)^    input: (.+)$", block[1])[1].split(", ")
    assert output.index("localrule picard_collect_multiple_metrics:") < output.index(
        "localrule mnase_qc_summary:"
    )


@pytest.mark.parametrize("picard", [True, "true"])
@pytest.mark.parametrize(
    "priority", [SUMMARY, METRICS], ids=["summary-first", "producer-first"]
)
def test_first_run_waits_for_picard_with_either_priority(
    tmp_path, project, snakemake_executable, picard, priority
):
    config = project(picard=picard, summary=False)
    result = invoke(
        tmp_path, snakemake_executable, config, producer=True, priority=priority
    )
    assert_summary(tmp_path, METRICS)
    assert_dependency(result)


@pytest.mark.parametrize("summary", [False, True])
def test_summary_alone_pulls_picard_and_then_is_up_to_date(
    tmp_path, project, snakemake_executable, summary
):
    config = project(summary=summary)
    result = invoke(tmp_path, snakemake_executable, config)
    assert_summary(tmp_path, METRICS)
    assert_dependency(result)
    before = (tmp_path / "picard-calls.jsonl").read_bytes()
    content = (tmp_path / SUMMARY).read_bytes()
    again = invoke(tmp_path, snakemake_executable, config)
    assert "Nothing to be done" in again.stdout + again.stderr
    assert (tmp_path / "picard-calls.jsonl").read_bytes() == before
    assert (tmp_path / SUMMARY).read_bytes() == content


@pytest.mark.parametrize(
    "picard,old", [(None, False), (False, False), ("false", False), (False, True)]
)
def test_disabled_picard_keeps_existence_based_field_without_reference(
    tmp_path, project, snakemake_executable, picard, old
):
    config = project(picard=picard, old_metrics=old, reference=False)
    result = invoke(tmp_path, snakemake_executable, config)
    assert_summary(tmp_path, METRICS if old else "NA")
    assert (
        "localrule picard_collect_multiple_metrics:"
        not in result.stdout + result.stderr
    )
    assert not (tmp_path / "picard-calls.jsonl").exists()


@pytest.mark.parametrize("mode", ["fail", "missing"])
def test_picard_failure_blocks_summary_in_a_fresh_directory(
    tmp_path, project, snakemake_executable, mode
):
    config = project(mode=mode)
    result = invoke(tmp_path, snakemake_executable, config, producer=True, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert (
        "PICARD_STUB_FAILURE" if mode == "fail" else "MissingOutputException"
    ) in output
    assert not (tmp_path / SUMMARY).exists()
    assert "localrule mnase_qc_summary:" not in output


def test_enabled_picard_still_rejects_missing_reference(
    tmp_path, project, snakemake_executable
):
    config = project(reference=False, validate=False)
    validation = json.loads((tmp_path / "validation.json").read_text())
    assert validation["returncode"] != 0
    assert "reference_fasta is missing" in validation["stderr"] + validation["stdout"]
    result = invoke(tmp_path, snakemake_executable, config, check=False)
    assert result.returncode != 0
    assert "reference_fasta is missing" in result.stdout + result.stderr
    assert not (tmp_path / SUMMARY).exists()
    assert not (tmp_path / "picard-calls.jsonl").exists()


def test_default_targets_stay_treatment_and_assay_scoped(
    tmp_path, project, snakemake_executable
):
    config = project(mixed=True)
    result = invoke(tmp_path, snakemake_executable, config, default_targets=True)
    output = result.stdout + result.stderr
    for rule, samples in (
        ("mnase_qc_summary", {"M1"}),
        ("picard_collect_multiple_metrics", {"M1", "T1"}),
    ):
        blocks = re.findall(rf"(?:local)?rule {rule}:\n((?:[ \t]+[^\n]*\n)+)", output)
        observed = {re.search(r"wildcards: sample=(\w+)", block)[1] for block in blocks}
        assert observed == samples
