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


@pytest.mark.parametrize("profile", ["legacy", "updated"])
def test_production_cross_correlation_and_summary(
    tmp_path, cc_input, tmp_config, snakemake_executable, profile
):
    from conftest import seed_run_identity

    bam, _, _ = cc_input
    _, config, _ = tmp_config(
        config={
            "config_schema_version": "2.0.0",
            "analysis_profile": profile,
            "outdir": "results",
            "threads": 1,
            "use_control": False,
            "blacklist_mode": "off",
            "organelle_mode": "keep",
            "multiqc": False,
            "replicate_analysis": False,
            "qc": {"cross_correlation": True},
            "genome_resources": {"tiny": {"effective_genome_size": 1000000}},
        },
        samples=(
            "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\trole\n"
            "S1\tunused.fq\t\tSE\tchipseq\tCTCF\tnarrow\ttiny\tindex\ttreatment\n"
        ),
    )
    seed_run_identity(
        tmp_path, snakemake_executable, REPO / "workflow/Snakefile", config
    )
    analysis = tmp_path / "results/S1/02_align/S1.analysis.bam"
    analysis.parent.mkdir(parents=True)
    shutil.copyfile(bam, analysis)
    summary = tmp_path / "results/multiqc/cross_correlation_summary.tsv"
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
            "--cores",
            "1",
            "--printshellcmds",
            "--allowed-rules",
            "cross_correlation",
            "cross_correlation_summary",
            "--",
            summary.relative_to(tmp_path),
        ],
    )
    output = result.stdout + result.stderr
    assert "-p=1" in output and "-rf" in output and "-speak=" not in output
    detail = tmp_path / "results/S1/05_qc/cross_correlation/S1.cross_correlation_qc.tsv"
    with detail.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["execution_status"] == "reported"
    assert row["interpretation_status"] == "unknown_target"
    assert row["threshold_hint"] == "not_evaluated"
    assert row["input_stage"] == "analysis_se" and row["window_upper"] == "100"
    assert int(row["candidate_count"]) >= 1
    raw_path = detail.with_name("S1.cc.qc")
    raw = raw_path.read_text().strip().split("\t")
    assert float(row["estimated_fragment_length"]) == float(raw[2].split(",")[0])
    assert float(row["nsc"]) == float(raw[8]) and float(row["rsc"]) == float(raw[9])
    with summary.open() as handle:
        project = next(csv.DictReader(handle, delimiter="\t"))
    for field in (
        "nsc",
        "rsc",
        "interpretation_status",
        "threshold_hint",
        "input_stage",
    ):
        assert project[field] == row[field]
    assert raw_path.with_name("S1.cc.plot.pdf").read_bytes().startswith(b"%PDF-")
    assert analysis.read_bytes() == bam.read_bytes()


def test_real_spp_phantom_and_non_degenerate_positive(tmp_path, cc_tools):
    from scripts.parse_cross_correlation import (
        parse_cc_qc_record,
        interpret_cross_correlation,
    )

    rng = random.Random(881)
    results = {}
    for kind, mix in [
        ("phantom", [(50, 20000)]),
        ("positive", [(200, 18000), (50, 2000)]),
    ]:
        tags = []
        for shift, count in mix:
            for _ in range(count):
                start = rng.randint(1000, 900000)
                tags.extend(
                    [(start, start + 50, "+"), (start + shift - 50, start + shift, "-")]
                )
        path = tmp_path / f"{kind}.tagAlign"
        path.write_text(
            "".join(
                f"chr1\t{s}\t{e}\tN\t1000\t{strand}\n" for s, e, strand in sorted(tags)
            )
        )
        for upper in (15, 60):
            qc = tmp_path / f"{kind}_{upper}.cc.qc"
            run(
                tmp_path,
                [cc_tools[0], f"-c={path}", f"-out={qc}", "-rf", f"-x=-500:{upper}"],
            )
            parsed = parse_cc_qc_record(qc)
            interpreted = interpret_cross_correlation(
                parsed,
                metadata={
                    "assay": "chipseq",
                    "target_type": "tf",
                    "read_length": 50,
                    "variable_read_length": False,
                },
            )
            results[kind, upper] = interpreted
    assert results["phantom", 15]["estimated_fragment_length"] == 50
    assert results["phantom", 15]["interpretation_status"] == "phantom_selected"
    assert results["phantom", 15]["threshold_hint"] == "not_evaluated"
    assert results["phantom", 60]["threshold_hint"] == "local_fail"
    assert results["phantom", 60]["estimated_fragment_length"] > 60
    for upper in (15, 60):
        assert results["positive", upper]["estimated_fragment_length"] == 200
        assert results["positive", upper]["rsc"] > 0.8
        assert results["positive", upper]["threshold_hint"] == "local_pass"


@pytest.mark.parametrize("profile", ["legacy", "updated"])
def test_empty_bam_production_rule_has_detail_without_fake_tool_artifacts(
    tmp_path, tmp_config, cc_tools, snakemake_executable, profile
):
    from conftest import seed_run_identity

    _, config, _ = tmp_config(
        config={
            "config_schema_version": "2.0.0",
            "analysis_profile": profile,
            "outdir": "results",
            "threads": 1,
            "use_control": False,
            "replicate_analysis": False,
            "blacklist_mode": "off",
            "organelle_mode": "keep",
            "multiqc": False,
            "qc": {"cross_correlation": True},
            "genome_resources": {"tiny": {"effective_genome_size": 1000000}},
        },
        samples="sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\tunused.fq\t\tSE\tchipseq\tCTCF\tnarrow\ttiny\tindex\n",
    )
    seed_run_identity(
        tmp_path, snakemake_executable, REPO / "workflow/Snakefile", config
    )
    bam = tmp_path / "results/S1/02_align/S1.analysis.bam"
    bam.parent.mkdir(parents=True)
    sam = tmp_path / "empty.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n")
    run(tmp_path, [cc_tools[1], "view", "-b", "-o", bam, sam])
    detail = tmp_path / "results/S1/05_qc/cross_correlation/S1.cross_correlation_qc.tsv"
    argv = [
        snakemake_executable,
        "-s",
        REPO / "workflow/Snakefile",
        "--workflow-profile",
        "none",
        "--configfile",
        config,
        "--cores",
        "1",
        "--allowed-rules",
        "cross_correlation",
        "--",
        detail.relative_to(tmp_path),
    ]
    run(tmp_path, argv)
    with detail.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["execution_status"] == row["interpretation_status"] == "no_data"
    assert row["nsc"] == row["rsc"] == "NA"
    assert not detail.with_name("S1.cc.qc").exists()
    assert not detail.with_name("S1.cc.plot.pdf").exists()
    # Snakemake can resume and can rebuild a deleted required detail normally.
    run(tmp_path, argv)
    detail.unlink()
    run(tmp_path, argv)
    assert detail.is_file()


@pytest.mark.parametrize("profile", ["legacy", "updated"])
def test_production_detail_policy_summary_manifest_adapter_chain(
    tmp_path, cc_tools, snakemake_executable, profile
):
    from test_qc_optional_argv_rules import (
        _write_project,
        _seed_identity,
        _write_summary_inputs,
        _write_stage_metrics,
    )
    from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
    from encode_pipeline.platform.adapters import (
        WorkflowInputs,
        QcSourceArtifact,
        QcSourceDocument,
    )
    from encode_pipeline.qc.contract import MITO_QC_V2_COLUMNS
    from scripts.qc_common import write_exact_tsv

    config = _write_project(
        tmp_path,
        assay="chipseq",
        cuttag_fragment_size=False,
        cross_correlation=True,
        profile=profile,
    )
    cfg = json.loads(config.read_text())
    cfg["organelle_mode"] = "keep"
    cfg["qc"]["summary"] = True
    config.write_text(json.dumps(cfg))
    samples = tmp_path / "config/samples.tsv"
    with samples.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        row = next(reader)
    row["layout"] = "SE"
    row["fastq_2"] = ""
    with samples.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    _seed_identity(tmp_path, snakemake_executable, config)
    _write_summary_inputs(tmp_path, "chipseq", profile=profile)
    base = tmp_path / "results/S1"
    stage = base / "05_qc/S1.blacklist_dropped.stage_metrics.json"
    _write_stage_metrics(stage)
    metrics = json.loads(stage.read_text())
    metrics.update(
        total_alignments=40000,
        kept_alignments=40000,
        dropped_alignments=0,
        dropped_templates=0,
    )
    stage.write_text(json.dumps(metrics))
    write_exact_tsv(
        base / "05_qc/S1.mito_fraction.tsv",
        MITO_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": "S1",
            "organelle_reads": "NA",
            "total_mapped_reads": "NA",
            "fraction": "NA",
            "mode": "keep",
            "status": "not_applicable",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "",
            "organelle_chroms_explicit": "false",
        },
    )
    # Execute the policy producer rather than leaving the historical staged policy.
    (base / "05_qc/S1.analysis_policy.json").unlink()
    rng = random.Random(881)
    records = []
    for shift, count in [(200, 18000), (50, 2000)]:
        for i in range(count):
            start = rng.randint(1000, 900000)
            records.extend(
                [
                    (start + 1, 0, f"{shift}-{i}a"),
                    (start + shift - 49, 16, f"{shift}-{i}b"),
                ]
            )
    sam = tmp_path / "positive.sam"
    with sam.open("w") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n")
        for pos, flag, name in sorted(records):
            handle.write(
                f"{name}\t{flag}\tchr1\t{pos}\t60\t50M\t*\t0\t0\t"
                + "A" * 50
                + "\t"
                + "I" * 50
                + "\n"
            )
    run(
        tmp_path,
        [cc_tools[1], "view", "-b", "-o", base / "02_align/S1.analysis.bam", sam],
    )
    run(
        tmp_path,
        [
            snakemake_executable,
            "-s",
            REPO / "workflow/Snakefile",
            "--workflow-profile",
            "none",
            "--configfile",
            config,
            "--cores",
            "1",
            "--allowed-rules",
            "analysis_template_metrics",
            "analysis_policy",
            "cross_correlation",
            "qc_summary",
            "cross_correlation_summary",
            "project_qc_summary",
            "--",
            "results/S1/01_qc/S1.qc_summary.tsv",
            "results/multiqc/cross_correlation_summary.tsv",
            "results/multiqc/project_qc_summary.tsv",
        ],
    )
    summary = base / "01_qc/S1.qc_summary.tsv"
    with summary.open() as handle:
        result = next(csv.DictReader(handle, delimiter="\t"))
    assert result["cross_correlation_status"] == "ok"
    assert float(result["cross_correlation_rsc"]) > 0.8
    policy = json.loads((base / "05_qc/S1.analysis_policy.json").read_text())
    assert policy["scientific_policy_revision"] == 6 and policy["profile"] == profile
    adapter = EncodeStyleWorkflowAdapter()
    inputs = WorkflowInputs(config={})
    artifacts = adapter.extract_artifacts(inputs, tmp_path)
    assert artifacts.is_success, artifacts.issues
    detail = [a for a in artifacts.value if a.output_type == "cross_correlation_qc"]
    assert (
        len(detail) == 1
        and detail[0].relative_path
        == "results/S1/05_qc/cross_correlation/S1.cross_correlation_qc.tsv"
    )
    source = QcSourceArtifact(
        artifact_id="summary",
        output_type="qc_summary",
        relative_path="results/S1/01_qc/S1.qc_summary.tsv",
        metadata={"scope": "sample", "sample_id": "S1", "assay": "chipseq"},
    )
    extracted = adapter.extract_qc_metrics(
        inputs, (QcSourceDocument(source=source, content=summary.read_bytes()),)
    )
    assert extracted.is_success, extracted.issues
    values = {m.metric_key: float(m.value) for m in extracted.value}
    assert values["cross_correlation.nsc"] == float(result["cross_correlation_nsc"])
    assert values["cross_correlation.rsc"] == float(result["cross_correlation_rsc"])


@pytest.mark.parametrize(
    "spans,expected", [([50, 51] * 250 + [100], 50), ([51, 52] * 250 + [100], 52)]
)
def test_native_r_round_and_late_variable_length(tmp_path, cc_tools, spans, expected):
    from scripts.make_cross_correlation_tagalign import _read_stats

    path = tmp_path / "variable.tagAlign"
    path.write_text(
        "".join(
            f"chr1\t{i * 200}\t{i * 200 + n}\tN\t1000\t+\n" for i, n in enumerate(spans)
        )
    )
    rscript = Path(cc_tools[0]).resolve().parent / "Rscript"
    result = run(
        tmp_path,
        [
            rscript,
            "-e",
            "d<-read.table(commandArgs(TRUE)[1],nrows=500); cat(round(median(d$V3-d$V2)))",
            path,
        ],
    )
    length, distribution, variable = _read_stats(path)
    assert int(result.stdout) == length == expected
    assert variable and sum(distribution.values()) == 501 and distribution["100"] == 1


def test_pe_r1_from_real_bam_retains_both_orientations(tmp_path, cc_tools):
    sam = tmp_path / "pairs.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        "a\t99\tchr1\t100\t60\t50M\t=\t250\t200\t*\t*\n"
        "a\t147\tchr1\t250\t60\t50M\t=\t100\t-200\t*\t*\n"
        "b\t163\tchr1\t250\t60\t50M\t=\t400\t200\t*\t*\n"
        "b\t83\tchr1\t400\t60\t50M\t=\t250\t-200\t*\t*\n"
    )
    bam = tmp_path / "pairs.bam"
    run(tmp_path, [cc_tools[1], "view", "-b", "-o", bam, sam])
    out = tmp_path / "tags.tsv"
    meta = tmp_path / "meta.json"
    run(
        tmp_path,
        [
            sys.executable,
            REPO / "scripts/make_cross_correlation_tagalign.py",
            "--analysis-bam",
            bam,
            "--samtools",
            cc_tools[1],
            "--layout",
            "PE",
            "--assay",
            "atac",
            "--target",
            "ATAC",
            "--target-type",
            "tf",
            "--sample",
            "S1",
            "--tagalign",
            out,
            "--metadata",
            meta,
        ],
    )
    assert out.read_text().splitlines() == [
        "chr1\t99\t149\tN\t1000\t+",
        "chr1\t399\t449\tN\t1000\t-",
    ]
    metadata = json.loads(meta.read_text())
    assert metadata["valid_templates"] == 2 and metadata["input_stage"] == "analysis_r1"
    assert metadata["target_type"] == "unknown" and metadata["window_upper"] == 100
