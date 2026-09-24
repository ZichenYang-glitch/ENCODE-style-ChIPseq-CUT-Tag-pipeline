"""Explicit Hi-TrAC tool qualification; never discovered by the generic fast job.

Requires reviewed runtime/reference bindings and tiny_inputs.py output through
HELIXWEAVE_HITRAC_* coordinates. Missing prerequisites fail, never skip. Fault
cases replace just one child executable after real runtime admission, without
changing the original tracPre2 script, command, assertions or scientific bytes.
"""

from collections import Counter
from dataclasses import replace
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from encode_pipeline.adapters.hitrac_preprocess import qualification as q
from encode_pipeline.adapters.hitrac_preprocess.admission import (
    Sample,
    load_runtime_binding,
)
from encode_pipeline.adapters.hitrac_preprocess.calls import verify_calls
from encode_pipeline.adapters.hitrac_preprocess.outputs import summary_columns

pytestmark = pytest.mark.real_execution


@pytest.fixture(scope="module")
def coordinates():
    keys = ["RUNTIME_BINDING", "REFERENCE_BINDING", "REFERENCE_SHA256", "TINY_INPUTS"]
    values = {key: os.environ["HELIXWEAVE_HITRAC_" + key] for key in keys}
    values["runtime"] = load_runtime_binding(Path(values["RUNTIME_BINDING"]))
    values["design"] = json.loads(
        (Path(values["TINY_INPUTS"]) / "design.json").read_text()
    )
    return values


def arguments(coordinates, tmp_path, scenario="positive", mapq=10):
    directory = Path(coordinates["TINY_INPUTS"]) / scenario / "fastq"
    samples = [
        Sample(
            "sample " + p.name.split("_R1")[0],
            p,
            directory / p.name.replace("_R1", "_R2"),
        )
        for p in sorted(directory.glob("*_R1.fastq.gz"))
    ]
    return dict(
        runtime_binding=Path(coordinates["RUNTIME_BINDING"]),
        reference_binding=Path(coordinates["REFERENCE_BINDING"]),
        reference_sha256=coordinates["REFERENCE_SHA256"],
        samples=samples,
        attempt=tmp_path / "attempt",
        threads=2,
        mapq=mapq,
        timeout=120,
    )


def save_result(tmp_path, result):
    (tmp_path / "test-observation.json").write_text(json.dumps(result, indent=2) + "\n")


def rows(path):
    with gzip.open(path, "rt") as stream:
        return [line.rstrip("\n").split("\t") for line in stream]


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expected_rows(design):
    result = []
    for number, pair in enumerate(
        design["scenarios"]["positive"]["samples"][0]["pairs"], 1
    ):
        a, b = pair["source_endpoints"]
        l1, l2 = pair["expected_linker_positions"]
        result.append(
            tuple(
                map(
                    str,
                    [
                        a["chrom"],
                        a["start"],
                        a["end"],
                        b["chrom"],
                        b["start"],
                        b["end"],
                        f"{number}_{l1}_{l2}",
                        42,
                        a["strand"],
                        b["strand"],
                    ],
                )
            )
        )
    return result


@pytest.mark.parametrize("mapq", [10, 17])
def test_original_cli_double_sample_preserves_pet_multiset_and_all_metrics(
    coordinates, tmp_path, mapq
):
    args = arguments(coordinates, tmp_path, mapq=mapq)
    request = args.copy()
    request.pop("attempt")
    request["samples"] = [
        {"id": s.id, "r1": str(s.r1), "r2": str(s.r2)} for s in args["samples"]
    ]
    request = {k: str(v) if isinstance(v, Path) else v for k, v in request.items()}
    request_path = tmp_path / "request.json"
    request_path.write_text(json.dumps(request))
    source = Path(__file__).resolve().parents[2]
    argv = [
        sys.executable,
        "-I",
        "-S",
        "-B",
        str(source / "scripts/qualify_hitrac_preprocess.py"),
        "--request",
        str(request_path),
        "--attempt",
        str(args["attempt"]),
    ]
    before = {str(p): file_hash(p) for s in args["samples"] for p in (s.r1, s.r2)}
    run = subprocess.run(
        argv, cwd=tmp_path, capture_output=True, text=True, timeout=180
    )
    (tmp_path / "cli.stdout").write_text(run.stdout)
    (tmp_path / "cli.stderr").write_text(run.stderr)
    save_result(tmp_path, {"argv": argv, "exit": run.returncode})
    assert run.returncode == 0, run.stdout + run.stderr
    outcome = json.loads(run.stdout)
    assert outcome["status"] == "complete"
    complete = json.loads((args["attempt"] / "complete.json").read_text())
    assert set(complete["results"]["samples"]) == {"s000001", "s000002"}
    expected = expected_rows(coordinates["design"])
    wanted = coordinates["design"]["scenarios"]["positive"][
        "conditional_expected_summary"
    ]
    for token, result in complete["results"]["samples"].items():
        root = args["attempt"] / "output" / token
        assert Counter(map(tuple, rows(root / (token + "_all.bedpe.gz")))) == Counter(
            expected
        )
        assert Counter(
            map(tuple, rows(root / (token + "_unique.bedpe.gz")))
        ) == Counter(expected[i] for i in (0, 3, 4, 5, 6))
        assert (
            result["all"] == 8 and result["noBg"] == 5 and result["qc_all_unique"] == 7
        )
        assert result["metrics"] == pytest.approx(wanted, rel=0, abs=1e-12)
        assert (root / (token + ".bam")).stat().st_size > 0
        assert all((root / (token + f"_R{i}.fastq.gz")).is_file() for i in (1, 2))
        # All artifact candidates are original BEDPE, never private BAM/FASTQ.
        assert all(a["path"].endswith(".bedpe.gz") for a in result["artifacts"])
    with (args["attempt"] / "output/tracPre_summary.txt").open() as stream:
        table = list(csv.reader(stream, delimiter="\t"))
    assert table[0] == [""] + summary_columns(mapq) and all(
        len(row) == 16 for row in table
    )
    receipts = verify_calls(args["attempt"] / "private/call-plan.json")
    assert len(receipts) == 15 and all(r["returncode"] == 0 for r in receipts)
    assert all(file_hash(Path(p)) == sha for p, sha in before.items())


@pytest.mark.parametrize(
    "scenario,reason,sample,collection",
    [
        ("noBg_empty", "empty_pet_set", "s000001", "noBg"),
        ("all_trans", "no_cis_denominator", "s000001", "all"),
        ("all_trans_no_linker", "no_cis_denominator", "s000001", "all"),
        ("all_unmapped", "empty_pet_set", "s000001", "all"),
        ("only_low_mapq", "empty_pet_set", "s000001", "all"),
        ("multisample_one_empty", "empty_pet_set", "s000002", "noBg"),
    ],
)
def test_real_zero_policy_refuses_whole_attempt_and_keeps_original_outputs(
    coordinates, tmp_path, scenario, reason, sample, collection
):
    args = arguments(coordinates, tmp_path, scenario)
    result = q.qualify(**args)
    save_result(tmp_path, result)
    assert result == dict(
        status="rejected", reason_code=reason, sample=sample, collection=collection
    )
    assert not (args["attempt"] / "complete.json").exists()
    assert (
        json.loads((args["attempt"] / "private/execution.json").read_text())[
            "returncode"
        ]
        == 0
    )
    assert (args["attempt"] / "output/tracPre_summary.txt").is_file()
    for index in range(1, len(args["samples"]) + 1):
        token = f"s{index:06d}"
        root = args["attempt"] / "output" / token
        assert (root / (token + ".bam")).is_file()
        assert (root / (token + "_all.bedpe.gz")).is_file()
        assert (root / (token + "_R1.fastq.gz")).is_file()


@pytest.mark.parametrize("scenario", ["empty", "all_trimmed"])
def test_original_zero_read_parser_failure_remains_failure(
    coordinates, tmp_path, scenario
):
    args = arguments(coordinates, tmp_path, scenario)
    result = q.qualify(**args)
    save_result(tmp_path, result)
    assert result["status"] == "rejected"
    assert result["reason_code"] == "call_missing"
    assert not (args["attempt"] / "complete.json").exists()
    assert (
        json.loads((args["attempt"] / "private/execution.json").read_text())[
            "returncode"
        ]
        == 1
    )
    assert (args["attempt"] / "output/s000001/s000001_R1.fastq.gz").is_file()


def test_real_singlemate_counterexample_is_refused_without_rewriting_science(
    coordinates, tmp_path
):
    args = arguments(coordinates, tmp_path, "singlemate_only")
    result = q.qualify(**args)
    save_result(tmp_path, result)
    assert result == dict(
        status="rejected",
        reason_code="pair_support_missing",
        sample="s000001",
        stage="pairing",
        collection="all",
    )
    assert not (args["attempt"] / "complete.json").exists()
    root = args["attempt"] / "output/s000001"
    assert rows(root / "s000001_all.bedpe.gz") == [
        ["chrA", "75000", "75120", "chrA", "75000", "75120", "2_-1_120", "42", "-", "-"]
    ]
    sam = subprocess.run(
        [
            str(coordinates["runtime"].tools["samtools"]),
            "view",
            str(root / "s000001.bam"),
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    (tmp_path / "singlemate.sam").write_text(sam)
    records = [line.split("\t") for line in sam.splitlines()]
    assert len(records) == 2 and len({row[0] for row in records}) == 2
    assert {int(row[1]) for row in records} == {73, 153}
    assert (
        json.loads((args["attempt"] / "private/execution.json").read_text())[
            "returncode"
        ]
        == 0
    )
    with (args["attempt"] / "output/tracPre_summary.txt").open() as stream:
        table = list(csv.reader(stream, delimiter="\t"))
    assert float(table[1][4]) == 1 and float(table[1][10]) == 1


def test_real_linker_distance_mapping_boundaries_remain_recorded(coordinates, tmp_path):
    args = arguments(coordinates, tmp_path, "boundaries")
    result = q.qualify(**args)
    save_result(tmp_path, result)
    # Third sample contains the real orphan artifact; the whole attempt rejects.
    assert result == dict(
        status="rejected",
        reason_code="pair_support_missing",
        sample="s000003",
        stage="pairing",
        collection="all",
    )
    output = args["attempt"] / "output"
    with gzip.open(output / "s000001/s000001_R1.fastq.gz", "rt") as stream:
        lines = stream.read().splitlines()
    trimmed = {lines[i][1:]: len(lines[i + 1]) for i in range(0, len(lines), 4)}
    assert "8_9_-1" not in trimmed and trimmed["9_10_-1"] == 10
    assert trimmed["4_120_-1"] == 120 and trimmed["5_120_-1"] == 120
    assert trimmed["6_-1_-1"] == 129 and trimmed["7_-1_-1"] == 129
    kept = {
        row[6].split("_")[0] for row in rows(output / "s000002/s000002_unique.bedpe.gz")
    }
    assert kept == {"1", "3", "4", "5", "6", "8"}
    all_rows = rows(output / "s000002/s000002_all.bedpe.gz")
    assert len(all_rows) == 8
    assert not (args["attempt"] / "complete.json").exists()


def inject_tool(monkeypatch, coordinates, tmp_path, stage, partial=False):
    runtime = coordinates["runtime"]
    tool = {
        "align": "bowtie2",
        "view": "samtools",
        "sort": "samtools",
        "bedpe": "bamToBed",
        "compress": "gzip",
        "qc_all": "cLoops2",
        "qc_noBg": "cLoops2",
        "remove_sam": "rm",
        "remove_qc": "rm",
    }[stage]
    target = tmp_path / "fault-tool"
    body = (
        """import json,pathlib,subprocess,sys
args=sys.argv[1:]
real=REAL
stage=STAGE
match=(stage=='align' or stage=='bedpe' or stage=='compress'
       or stage in ('view','sort') and args[0]==stage
       or stage=='remove_sam' and args[0].endswith('.sam')
       or stage=='remove_qc' and args[0].endswith('_bedpeQc.txt')
       or stage=='qc_all' and args[args.index('-o')+1].endswith('/allBedpeQc')
       or stage=='qc_noBg' and args[args.index('-o')+1].endswith('/uniNonBgBedpeQc'))
if stage not in ('qc_all','qc_noBg','remove_qc'):
    match=match and any('/s000001/' in arg for arg in args)
code=subprocess.call([real,*args])
if match and code==0:
    if PARTIAL:
        bam=pathlib.Path(args[args.index('-o')+1])
        raw=subprocess.run([real,'view','-h',str(bam)],capture_output=True,text=True,check=True).stdout
        removed=[line for line in raw.splitlines(True) if not line.startswith('@') and line.split('\\t')[0]=='5_-1_-1']
        assert len(removed)==2
        EVIDENCE.write_bytes(bam.read_bytes())
        data=''.join(line for line in raw.splitlines(True) if line not in removed)
        subprocess.run([real,'view','-b','-o',str(bam),'-'],input=data,text=True,check=True)
    sys.exit(73)
sys.exit(code)
""".replace("REAL", repr(str(runtime.tools[tool])))
        .replace("STAGE", repr(stage))
        .replace("PARTIAL", repr(partial))
        .replace("EVIDENCE", f"pathlib.Path({str(tmp_path / 'before-drop.bam')!r})")
    )
    target.write_text("#!" + str(runtime.python) + "\n" + body)
    target.chmod(0o700)
    altered = replace(runtime, tools=runtime.tools | {tool: target})
    # Fault injection is confined to the admitted executable dependency, not the
    # actual tool, whole route, upstream script or completion assertions.
    monkeypatch.setattr(q, "load_runtime_binding", lambda _: altered)


@pytest.mark.parametrize(
    "stage",
    [
        "align",
        "view",
        "sort",
        "bedpe",
        "compress",
        "qc_all",
        "qc_noBg",
        "remove_sam",
        "remove_qc",
    ],
)
def test_each_original_child_nonzero_is_rejected_even_with_complete_output(
    coordinates, tmp_path, monkeypatch, stage
):
    inject_tool(monkeypatch, coordinates, tmp_path, stage)
    args = arguments(coordinates, tmp_path)
    result = q.qualify(**args)
    save_result(tmp_path, result)
    expected = dict(status="rejected", reason_code="child_failed", stage=stage)
    if stage not in ("qc_all", "qc_noBg", "remove_qc"):
        expected["sample"] = "s000001"
    assert result == expected
    assert not (args["attempt"] / "complete.json").exists()
    assert (args["attempt"] / "output/tracPre_summary.txt").exists()


def test_original_partial_bam_exit73_counterexample_now_refused(
    coordinates, tmp_path, monkeypatch
):
    inject_tool(monkeypatch, coordinates, tmp_path, "view", partial=True)
    args = arguments(coordinates, tmp_path)
    result = q.qualify(**args)
    save_result(tmp_path, result)
    assert result == dict(
        status="rejected", reason_code="child_failed", sample="s000001", stage="view"
    )
    assert not (args["attempt"] / "complete.json").exists()
    assert len(rows(args["attempt"] / "output/s000001/s000001_all.bedpe.gz")) == 7
    assert len(rows(args["attempt"] / "output/s000001/s000001_unique.bedpe.gz")) == 4
    assert (tmp_path / "before-drop.bam").is_file()
    assert (
        json.loads((args["attempt"] / "private/execution.json").read_text())[
            "returncode"
        ]
        == 0
    )


@pytest.mark.parametrize("mode", ["cancel", "timeout"])
def test_qualification_retains_evidence_after_cancel_or_timeout(
    coordinates, tmp_path, mode
):
    args = arguments(coordinates, tmp_path)
    if mode == "timeout":
        args["timeout"] = 0.01
    else:
        # Event-driven: wait for original child launch, not arbitrary sleep order.
        args["cancelled"] = lambda: bool(
            list((args["attempt"] / "private/calls").glob("*.start.json"))
        )
    result = q.qualify(**args)
    save_result(tmp_path, result)
    assert result == dict(
        status="rejected", reason_code="cancelled" if mode == "cancel" else "timed_out"
    )
    assert not (args["attempt"] / "complete.json").exists()
    assert (args["attempt"] / "input").is_dir() or (args["attempt"] / "inputs").is_dir()
    record = json.loads((args["attempt"] / "private/execution.json").read_text())
    with pytest.raises(ProcessLookupError):
        os.killpg(record["pid"], 0)
