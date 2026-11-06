"""Unit tests for scripts/calc_organelle_metrics.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts/calc_organelle_metrics.py"
_SPEC = importlib.util.spec_from_file_location("calc_organelle_metrics", _SCRIPT)
om = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(om)


def test_parse_args():
    args = om.parse_args(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "exclude",
            "--chroms",
            "chrM",
            "chrMT",
            "--is-explicit",
            "--output",
            "out.tsv",
        ]
    )
    assert args.bam == "test.bam"
    assert args.sample == "S1"
    assert args.mode == "exclude"
    assert args.chroms == ["chrM", "chrMT"]
    assert args.is_explicit is True
    assert args.output == "out.tsv"


class MockPopenCommunicate:
    def __init__(self, stdout: str, stderr: str = "", returncode: int = 0):
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode

    def communicate(self):
        return self.stdout, self.stderr


def test_run_idxstats_mock(monkeypatch):
    idxstats_text = "chr1\t1000\t500\t10\nchrM\t200\t50\t0\n*\t0\t0\t5\n"
    monkeypatch.setattr(
        om.subprocess, "Popen", lambda *args, **kw: MockPopenCommunicate(idxstats_text)
    )
    counts = om.run_idxstats("mock.bam")
    assert counts == {"chr1": 500, "chrM": 50}


def test_calc_organelle_metrics_modes(tmp_path: Path, monkeypatch):
    # Mock run_idxstats
    def mock_idxstats(bam):
        return {"chr1": 800, "chr2": 200, "chrM": 100}

    monkeypatch.setattr(om, "run_idxstats", mock_idxstats)

    # 1. Mode report: chrM resolved (100 / 1100)
    out1 = tmp_path / "out1.tsv"
    ret1 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "report",
            "--chroms",
            "chrM",
            "--output",
            str(out1),
        ]
    )
    assert ret1 == 0
    lines1 = out1.read_text(encoding="utf-8").strip().splitlines()
    assert lines1[1].split("\t") == [
        "2",
        "S1",
        "100",
        "1100",
        "0.090909",
        "report",
        "ok",
        "pre_duplicate_mapq_bam",
        "chrM",
        "false",
    ]

    # 2. Mode report: unresolved contigs
    out2 = tmp_path / "out2.tsv"
    ret2 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "report",
            "--chroms",
            "chrMT",
            "--output",
            str(out2),
        ]
    )
    assert ret2 == 0
    lines2 = out2.read_text(encoding="utf-8").strip().splitlines()
    assert lines2[1].split("\t") == [
        "2",
        "S1",
        "NA",
        "1100",
        "NA",
        "report",
        "unresolved_contigs",
        "pre_duplicate_mapq_bam",
        "NA",
        "false",
    ]

    # 3. Mode exclude: success
    out3 = tmp_path / "out3.tsv"
    ret3 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "exclude",
            "--chroms",
            "chrM",
            "--output",
            str(out3),
        ]
    )
    assert ret3 == 0
    lines3 = out3.read_text(encoding="utf-8").strip().splitlines()
    assert lines3[1].split("\t") == [
        "2",
        "S1",
        "100",
        "1100",
        "0.090909",
        "exclude",
        "ok",
        "pre_duplicate_mapq_bam",
        "chrM",
        "false",
    ]

    # 4. Mode exclude: no matches fails closed
    out4 = tmp_path / "out4.tsv"
    ret4 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "exclude",
            "--chroms",
            "chrMT_nonexistent",
            "--output",
            str(out4),
        ]
    )
    assert ret4 == 1

    # 5. Mode exclude: explicit missing contig fails closed
    ret5 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "exclude",
            "--chroms",
            "chrM_missing",
            "--is-explicit",
            "--output",
            str(out4),
        ]
    )
    assert ret5 == 1

    # 6. Mode keep: resolved
    out6 = tmp_path / "out6.tsv"
    ret6 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "keep",
            "--chroms",
            "chrM",
            "--output",
            str(out6),
        ]
    )
    assert ret6 == 0
    lines6 = out6.read_text(encoding="utf-8").strip().splitlines()
    assert lines6[1].split("\t") == [
        "2",
        "S1",
        "100",
        "1100",
        "0.090909",
        "keep",
        "ok",
        "pre_duplicate_mapq_bam",
        "chrM",
        "false",
    ]

    # 7. Zero total mapped reads
    monkeypatch.setattr(om, "run_idxstats", lambda bam: {"chr1": 0, "chrM": 0})
    out7 = tmp_path / "out7.tsv"
    ret7 = om.main(
        [
            "--bam",
            "test.bam",
            "--sample",
            "S1",
            "--mode",
            "report",
            "--chroms",
            "chrM",
            "--output",
            str(out7),
        ]
    )
    assert ret7 == 0
    lines7 = out7.read_text(encoding="utf-8").strip().splitlines()
    assert lines7[1].split("\t") == [
        "2",
        "S1",
        "NA",
        "0",
        "NA",
        "report",
        "no_mapped_reads",
        "pre_duplicate_mapq_bam",
        "chrM",
        "false",
    ]
