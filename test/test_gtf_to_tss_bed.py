#!/usr/bin/env python3
"""Unit tests for scripts/gtf_to_tss_bed.py."""

import os
import shutil
import sys
import tempfile

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO, "scripts"))

from gtf_to_tss_bed import collect_tss_records, main  # noqa: E402


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def test_transcript_tss_plus_and_minus():
    td = tempfile.mkdtemp(prefix="tss_unit_")
    try:
        gtf = os.path.join(td, "anno.gtf")
        _write(
            gtf,
            "\n".join([
                'chr1\tsrc\ttranscript\t10\t50\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";',
                'chr1\tsrc\ttranscript\t100\t180\t.\t-\t.\tgene_id "g2"; transcript_id "tx2";',
                "",
            ]),
        )
        assert collect_tss_records(gtf) == [
            ("chr1", 9, 10, "tx1", "0", "+"),
            ("chr1", 179, 180, "tx2", "0", "-"),
        ]
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_gene_fallback_when_no_transcripts():
    td = tempfile.mkdtemp(prefix="tss_unit_")
    try:
        gtf = os.path.join(td, "anno.gtf")
        _write(
            gtf,
            'chr2\tsrc\tgene\t5\t25\t.\t-\t.\tgene_id "geneA";\n',
        )
        assert collect_tss_records(gtf) == [
            ("chr2", 24, 25, "geneA", "0", "-"),
        ]
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_transcripts_preferred_over_gene_records():
    td = tempfile.mkdtemp(prefix="tss_unit_")
    try:
        gtf = os.path.join(td, "anno.gtf")
        _write(
            gtf,
            "\n".join([
                'chr1\tsrc\tgene\t1\t90\t.\t+\t.\tgene_id "geneA";',
                'chr1\tsrc\ttranscript\t20\t90\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";',
                "",
            ]),
        )
        assert collect_tss_records(gtf) == [
            ("chr1", 19, 20, "txA", "0", "+"),
        ]
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_cli_writes_bed6():
    td = tempfile.mkdtemp(prefix="tss_unit_")
    try:
        gtf = os.path.join(td, "anno.gtf")
        out = os.path.join(td, "tss.bed")
        _write(
            gtf,
            'chr1\tsrc\ttranscript\t10\t50\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n',
        )
        assert main(["--gtf", gtf, "--output", out]) == 0
        with open(out) as fh:
            assert fh.read() == "chr1\t9\t10\ttx1\t0\t+\n"
    finally:
        shutil.rmtree(td, ignore_errors=True)


def run_all():
    tests = [
        test_transcript_tss_plus_and_minus,
        test_gene_fallback_when_no_transcripts,
        test_transcripts_preferred_over_gene_records,
        test_cli_writes_bed6,
    ]
    for test in tests:
        test()
    print(f"OK: {len(tests)}/{len(tests)} gtf_to_tss_bed tests passed")


if __name__ == "__main__":
    run_all()
