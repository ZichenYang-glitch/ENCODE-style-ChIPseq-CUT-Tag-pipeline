"""Unit tests for scripts/filter_blacklist_templates.py."""

from __future__ import annotations

import importlib.util
import io
from pathlib import Path
import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts/filter_blacklist_templates.py"
_SPEC = importlib.util.spec_from_file_location("filter_blacklist_templates", _SCRIPT)
fb = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fb)


def test_cigar_ref_length():
    assert fb.cigar_ref_length("50M") == 50
    assert fb.cigar_ref_length("20M5D30M") == 55
    assert fb.cigar_ref_length("20M5I30M") == 50
    assert fb.cigar_ref_length("10S40M10S") == 40
    assert fb.cigar_ref_length("100=") == 100
    assert fb.cigar_ref_length("10X40M") == 50


def test_validate_reference_contig_lengths(tmp_path: Path):
    sq = {"chr1": 1000, "chr2": 500}

    # Matching chrom_sizes
    cs_good = tmp_path / "good.sizes"
    cs_good.write_text("chr1\t1000\nchr2\t500\n", encoding="utf-8")
    fb.validate_reference_contig_lengths(sq, chrom_sizes_path=str(cs_good), sample="S1")

    # Mismatched chrom_sizes
    cs_bad = tmp_path / "bad.sizes"
    cs_bad.write_text("chr1\t999\nchr2\t500\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="BAM @SQ contig length mismatch for 'chr1'"):
        fb.validate_reference_contig_lengths(
            sq, chrom_sizes_path=str(cs_bad), sample="S1"
        )

    # Matching FAI
    fai_good = tmp_path / "good.fai"
    fai_good.write_text(
        "chr1\t1000\t0\t50\t51\nchr2\t500\t0\t50\t51\n", encoding="utf-8"
    )
    fb.validate_reference_contig_lengths(sq, fai_path=str(fai_good), sample="S1")

    # Mismatched FAI
    fai_bad = tmp_path / "bad.fai"
    fai_bad.write_text(
        "chr1\t1001\t0\t50\t51\nchr2\t500\t0\t50\t51\n", encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="BAM @SQ contig length mismatch for 'chr1'"):
        fb.validate_reference_contig_lengths(sq, fai_path=str(fai_bad), sample="S1")


def test_validate_organelle_contigs():
    sq = {"chr1": 1000, "chrM": 200}

    # Mode keep returns empty set
    assert fb.validate_organelle_contigs(sq, ["chrM"], False, "keep", "S1") == set()

    # Mode report returns matched contigs
    assert fb.validate_organelle_contigs(
        sq, ["chrM", "chrMT"], False, "report", "S1"
    ) == {"chrM"}

    # Mode report with no match returns empty set without failing
    assert fb.validate_organelle_contigs(sq, ["chrMT"], False, "report", "S1") == set()

    # Mode exclude with match returns matched contigs
    assert fb.validate_organelle_contigs(
        sq, ["chrM", "chrMT"], False, "exclude", "S1"
    ) == {"chrM"}

    # Mode exclude with no match fails closed
    with pytest.raises(RuntimeError, match="none of candidate contigs"):
        fb.validate_organelle_contigs(sq, ["chrMT"], False, "exclude", "S1")

    # Explicit candidate missing fails closed
    with pytest.raises(
        RuntimeError, match="explicit organelle contig 'chrM_missing' not found"
    ):
        fb.validate_organelle_contigs(
            sq, ["chrM", "chrM_missing"], True, "exclude", "S1"
        )


def test_validate_and_load_blacklist(tmp_path: Path):
    sq = {"chr1": 1000, "chr2": 500}

    # Valid BED (tab and space)
    bed = tmp_path / "valid.bed"
    bed.write_text("chr1\t10\t50\nchr2 100 200\n", encoding="utf-8")
    bl = fb.validate_and_load_blacklist(str(bed), sq, "S1")
    assert bl.total_intervals == 2

    # Nonexistent
    with pytest.raises(RuntimeError, match="blacklist file not found"):
        fb.validate_and_load_blacklist(str(tmp_path / "missing.bed"), sq, "S1")

    # Fewer than 3 fields
    bad_fields = tmp_path / "bad_fields.bed"
    bad_fields.write_text("chr1\t10\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="fewer than 3 fields"):
        fb.validate_and_load_blacklist(str(bad_fields), sq, "S1")

    # Non-integer coords
    bad_int = tmp_path / "bad_int.bed"
    bad_int.write_text("chr1\tfoo\t50\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="non-integer coordinates"):
        fb.validate_and_load_blacklist(str(bad_int), sq, "S1")

    # start >= end
    bad_span = tmp_path / "bad_span.bed"
    bad_span.write_text("chr1\t50\t50\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="invalid interval start=50 >= end=50"):
        fb.validate_and_load_blacklist(str(bad_span), sq, "S1")

    # Contig not in @SQ
    bad_chrom = tmp_path / "bad_chrom.bed"
    bad_chrom.write_text("chrX\t10\t50\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="not found in BAM @SQ"):
        fb.validate_and_load_blacklist(str(bad_chrom), sq, "S1")

    # Interval exceeds contig length
    bad_len = tmp_path / "bad_len.bed"
    bad_len.write_text("chr1\t900\t1100\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="exceeds contig length"):
        fb.validate_and_load_blacklist(str(bad_len), sq, "S1")

    # Empty BED
    empty = tmp_path / "empty.bed"
    empty.write_text("# comment\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="empty or contains no valid intervals"):
        fb.validate_and_load_blacklist(str(empty), sq, "S1")


def test_parse_args():
    args = fb.parse_args(
        [
            "--bam",
            "input.bam",
            "--output",
            "output.bam",
            "--dropped-count",
            "drop.tsv",
            "--sample",
            "S1",
            "--layout",
            "PE",
            "--exclude-dup",
            "--organelle-mode",
            "exclude",
            "--organelle-chroms",
            "chrM,chrMT",
            "--blacklist",
            "bl.bed",
            "--chrom-sizes",
            "chrom.sizes",
            "--fai",
            "ref.fa.fai",
            "--threads",
            "4",
        ]
    )
    assert args.bam == "input.bam"
    assert args.output == "output.bam"
    assert args.layout == "PE"
    assert args.exclude_dup is True
    assert args.organelle_mode == "exclude"
    assert args.organelle_chroms == ["chrM,chrMT"]
    assert args.blacklist == "bl.bed"
    assert args.chrom_sizes == "chrom.sizes"
    assert args.fai == "ref.fa.fai"
    assert args.threads == 4


class MockProcess:
    def __init__(self, lines: list[str], returncode: int = 0):
        self.stdout = lines
        self.stderr = io.StringIO("")
        self.returncode = returncode

    def wait(self):
        return self.returncode


class MockPopenCommunicate:
    def __init__(self, stdout: str, stderr: str = "", returncode: int = 0):
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode

    def communicate(self):
        return self.stdout, self.stderr


def test_read_bam_header_mock(monkeypatch):
    header_text = (
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n@PG\tID:bwa\n"
    )
    monkeypatch.setattr(
        fb.subprocess, "Popen", lambda *args, **kw: MockPopenCommunicate(header_text)
    )
    sq = fb.read_bam_header("mock.bam")
    assert sq == {"chr1": 1000, "chr2": 500}


class MockWriterProcess:
    def __init__(self):
        self.buf: list[str] = []
        self.stderr = io.StringIO("")

        class InStream:
            def __init__(inner, outer):
                inner.outer = outer

            def write(inner, text):
                inner.outer.buf.append(text)

            def close(inner):
                pass

        self.stdin = InStream(self)

    def wait(self):
        return 0


def test_write_filtered_bam_mock(monkeypatch, tmp_path: Path):
    reader_lines = [
        "@HD\tVN:1.6\n",
        "keep1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tA\tI\n",
        "drop1\t0\tchr1\t200\t60\t50M\t*\t0\t0\tA\tI\n",
    ]
    writer = MockWriterProcess()

    def mock_popen(cmd, *args, **kw):
        if (
            len(cmd) >= 3
            and cmd[0] == "samtools"
            and cmd[1] == "view"
            and cmd[2] == "-h"
        ):
            return MockProcess(reader_lines)
        return writer

    monkeypatch.setattr(fb.subprocess, "Popen", mock_popen)
    dropped, kept, by_reason = fb.write_filtered_bam(
        bam="in.bam",
        output=str(tmp_path / "out.bam"),
        drop_qnames={"drop1"},
        drop_reasons={"drop1": "blacklist"},
    )
    assert dropped == 1
    assert kept == 1
    assert by_reason == {"blacklist": 1}
    written = "".join(writer.buf)
    assert "keep1" in written
    assert "drop1" not in written


def test_scan_drop_qnames_se(monkeypatch, tmp_path: Path):
    bed = tmp_path / "bl.bed"
    bed.write_text("chr1\t390\t450\n", encoding="utf-8")
    bl = fb.Blacklist(str(bed))

    # SE mock alignments
    lines = [
        "@HD\tVN:1.6\n",
        "normal\t0\tchr1\t100\t60\t50M\t*\t0\t0\tAAAAA\tIIIII\n",
        "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAA\tIIIII\n",
        "dup\t1024\tchr1\t200\t60\t50M\t*\t0\t0\tAAAAA\tIIIII\n",
        "mito\t0\tchrM\t300\t60\t50M\t*\t0\t0\tAAAAA\tIIIII\n",
        "bl_hit\t0\tchr1\t400\t60\t50M\t*\t0\t0\tAAAAA\tIIIII\n",
    ]
    monkeypatch.setattr(
        fb.subprocess, "Popen", lambda *args, **kwargs: MockProcess(lines)
    )

    drop_qnames, drop_reasons, total = fb.scan_drop_qnames(
        bam="mock.bam",
        layout="SE",
        exclude_dup=True,
        matched_organelle={"chrM"},
        blacklist=bl,
    )
    assert total == 5
    assert "normal" not in drop_qnames
    assert drop_reasons["unmapped"] == "invalid"
    assert drop_reasons["dup"] == "duplicate"
    assert drop_reasons["mito"] == "organelle"
    assert drop_reasons["bl_hit"] == "blacklist"


def test_scan_drop_qnames_pe(monkeypatch, tmp_path: Path):
    bed = tmp_path / "bl.bed"
    bed.write_text("chr1\t420\t430\n", encoding="utf-8")
    bl = fb.Blacklist(str(bed))

    lines = [
        "@HD\tVN:1.6\n",
        # normal pair
        "normal\t99\tchr1\t100\t60\t50M\t=\t200\t150\tA\tI\n",
        "normal\t147\tchr1\t200\t60\t50M\t=\t100\t-150\tA\tI\n",
        # duplicate pair
        "dup\t1123\tchr1\t100\t60\t50M\t=\t200\t150\tA\tI\n",
        "dup\t1171\tchr1\t200\t60\t50M\t=\t100\t-150\tA\tI\n",
        # organelle pair
        "mito\t99\tchrM\t100\t60\t50M\t=\t200\t150\tA\tI\n",
        "mito\t147\tchrM\t200\t60\t50M\t=\t100\t-150\tA\tI\n",
        # blacklist hit
        "bl_hit\t99\tchr1\t400\t60\t50M\t=\t500\t150\tA\tI\n",
        "bl_hit\t147\tchr1\t500\t60\t50M\t=\t400\t-150\tA\tI\n",
        # orphan / incomplete mate
        "orphan\t73\tchr1\t600\t60\t50M\t*\t0\t0\tA\tI\n",
        # coord mismatch
        "bad_coord\t99\tchr1\t700\t60\t50M\t=\t850\t150\tA\tI\n",
        "bad_coord\t147\tchr1\t800\t60\t50M\t=\t700\t-150\tA\tI\n",
    ]
    monkeypatch.setattr(
        fb.subprocess, "Popen", lambda *args, **kwargs: MockProcess(lines)
    )

    drop_qnames, drop_reasons, total = fb.scan_drop_qnames(
        bam="mock.bam",
        layout="PE",
        exclude_dup=True,
        matched_organelle={"chrM"},
        blacklist=bl,
    )
    assert total == 11
    assert "normal" not in drop_qnames
    assert drop_reasons["dup"] == "duplicate"
    assert drop_reasons["mito"] == "organelle"
    assert drop_reasons["bl_hit"] == "blacklist"
    assert drop_reasons["orphan"] == "invalid"
    assert drop_reasons["bad_coord"] == "invalid"


def test_main_workflow(tmp_path: Path, monkeypatch):
    bam = tmp_path / "input.bam"
    bam.touch()
    out = tmp_path / "out.bam"
    drop_tsv = tmp_path / "drop.tsv"
    bl_file = tmp_path / "bl.bed"
    bl_file.write_text("chr1\t100\t200\n", encoding="utf-8")

    monkeypatch.setattr(fb, "read_bam_header", lambda b: {"chr1": 1000, "chrM": 200})
    monkeypatch.setattr(fb, "validate_reference_contig_lengths", lambda **kw: None)
    monkeypatch.setattr(
        fb,
        "scan_drop_qnames",
        lambda **kw: ({"drop1"}, {"drop1": "blacklist"}, 4),
    )
    monkeypatch.setattr(
        fb,
        "write_filtered_bam",
        lambda **kw: (2, 2, {"blacklist": 2}),
    )

    ret = fb.main(
        [
            "--bam",
            str(bam),
            "--output",
            str(out),
            "--dropped-count",
            str(drop_tsv),
            "--sample",
            "S1",
            "--layout",
            "PE",
            "--blacklist",
            str(bl_file),
        ]
    )
    assert ret == 0
    assert drop_tsv.exists()
    assert (tmp_path / "drop.stage_metrics.json").exists()
