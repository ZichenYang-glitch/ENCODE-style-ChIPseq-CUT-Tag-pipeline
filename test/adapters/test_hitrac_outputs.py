"""Fast synthetic contracts for the private Hi-TrAC output consumer.

BAM bytes are placeholders here: independent pair-source checks and real-tool
qualification exercise BAM validity elsewhere. Nothing in this fixture runs QC.
"""

import csv
import gzip
import hashlib

import pytest

from encode_pipeline.adapters.hitrac_preprocess.outputs import (
    OutputRejected,
    iter_bedpe,
    summary_columns,
    verify_outputs,
)


ROWS = [
    ("chrA", "10", "20", "chrA", "2010", "2020", "1_10_-1", "42", "+", "-"),
    ("chrA", "10", "20", "chrA", "2010", "2020", "2_10_-1", "42", "+", "-"),
    ("chrB", "10", "20", "chrA", "4010", "4020", "3_10_-1", "42", "+", "-"),
    ("chrA", "20", "30", "chrA", "520", "530", "4_10_-1", "42", "+", "+"),
]
CONTIGS = {"chrA": 100000, "chrB": 100000}
METRICS = [4, 4, 0.75, 4, 0.25, 2 / 3, 0.5, 0.5, 0, 3, 0.75, 2 / 3, 0.5, 0.5, 0]


def _bedpe(path, rows):
    with gzip.open(path, "wt", encoding="ascii") as stream:
        stream.writelines("\t".join(row) + "\n" for row in rows)


def _summary(root, data, mapq=10):
    with (root / "tracPre_summary.txt").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow([""] + summary_columns(mapq))
        for sample, metrics in data.items():
            writer.writerow([sample] + metrics)


def _fixture(tmp_path, samples=("s000001",), mapq=10):
    root = tmp_path / "output"
    root.mkdir()
    work = tmp_path / "private"
    work.mkdir()
    log = []
    for sample in samples:
        directory = root / sample
        directory.mkdir()
        (directory / f"{sample}.bam").write_bytes(b"unit-test-only")
        for mate in (1, 2):
            with gzip.open(directory / f"{sample}_R{mate}.fastq.gz", "wt") as stream:
                for row in ROWS:
                    stream.write(f"@{row[6]}\nACGTACGTAC\n+\nIIIIIIIIII\n")
        rows = [(*row[:7], str(max(42, mapq)), *row[8:]) for row in ROWS]
        _bedpe(directory / f"{sample}_all.bedpe.gz", rows)
        _bedpe(directory / f"{sample}_unique.bedpe.gz", [rows[i] for i in (0, 2, 3)])
        log.append(
            f"2000-01-01 cLoops2 INFO FLAG_A:{sample}\n"
            "4 reads; of these:\n"
            "  4 (100.00%) were paired; of these:\n"
            "75.00% overall alignment rate\nFLAG_A\n"
        )
    (root / "2000-01-01_tracPre2.py.log").write_text("\n".join(log))
    _summary(root, {sample: METRICS for sample in samples}, mapq)
    return root, work, {sample: 4 for sample in samples}


def _verify(fixture, mapq=10):
    root, work, counts = fixture
    return verify_outputs(root, counts, CONTIGS, mapq, work)


def _file(root, suffix, sample="s000001"):
    return root / sample / f"{sample}_{suffix}.bedpe.gz"


def _fingerprints(root):
    # The synthetic fixture tree contains only this test's outputs.
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*")
        if path.is_file()
    }


@pytest.mark.parametrize("mapq", [0, 10, 30, 255])
def test_final_metrics_dynamic_header_and_private_preservation(tmp_path, mapq):
    fixture = _fixture(tmp_path, samples=("s000001", "s000002"), mapq=mapq)
    before = _fingerprints(fixture[0])
    result = _verify(fixture, mapq)
    assert set(result["samples"]) == {"s000001", "s000002"}
    for sample, record in result["samples"].items():
        assert record["metrics"] == METRICS
        assert (record["all"], record["noBg"], record["qc_all_unique"]) == (4, 3, 3)
        assert [item["path"] for item in record["artifacts"]] == [
            f"{sample}/{sample}_all.bedpe.gz",
            f"{sample}/{sample}_unique.bedpe.gz",
        ]
        assert all(len(item["sha256"]) == 64 for item in record["artifacts"])
    assert result["summary"]["path"] == "tracPre_summary.txt"
    assert before == _fingerprints(fixture[0])
    assert not list(fixture[0].glob("*complete*"))


@pytest.mark.parametrize("mapq", [True, "10", -1, 256, 10.0])
def test_invalid_mapq_is_not_coerced(mapq):
    with pytest.raises(OutputRejected, match="^invalid_mapq$"):
        summary_columns(mapq)


@pytest.mark.parametrize(
    "case", ["intermediate_header", "wrong_mapq", "extra_sample", "duplicate"]
)
def test_final_summary_shape_rejects_intermediate_or_wrong_sample_set(tmp_path, case):
    fixture = _fixture(tmp_path)
    summary = fixture[0] / "tracPre_summary.txt"
    text = summary.read_text()
    if case == "intermediate_header":
        text = "\tTotalRawReads\tMappingRatio(%s)\ns000001\t4\t75\n"
    elif case == "wrong_mapq":
        text = text.replace("mapq>=10", "mapq>=11")
    elif case == "extra_sample":
        text += text.splitlines(True)[1].replace("s000001", "s000003")
    else:
        text += text.splitlines(True)[1]
    summary.write_text(text)
    with pytest.raises(OutputRejected) as caught:
        _verify(fixture)
    assert caught.value.reason_code in {"final_summary_header", "summary_samples"}


@pytest.mark.parametrize("column,value", [(0, 5), (2, 0.9), (3, 3), (5, 1), (10, 1)])
def test_count_ratio_and_mapping_log_consistency(tmp_path, column, value):
    fixture = _fixture(tmp_path)
    metrics = METRICS.copy()
    metrics[column] = value
    _summary(fixture[0], {"s000001": metrics})
    with pytest.raises(OutputRejected, match="^summary_count_ratio_mismatch$"):
        _verify(fixture)


@pytest.mark.parametrize("bad", ["nan", "inf", "-inf"])
def test_summary_rejects_nonfinite(tmp_path, bad):
    fixture = _fixture(tmp_path)
    metrics = METRICS.copy()
    metrics[2] = bad
    _summary(fixture[0], {"s000001": metrics})
    with pytest.raises(OutputRejected, match="^summary_nonfinite$"):
        _verify(fixture)


def test_display_rounding_cannot_hide_fractional_count(tmp_path):
    fixture = _fixture(tmp_path)
    metrics = METRICS.copy()
    metrics[3] = "4.0000000000000001"
    _summary(fixture[0], {"s000001": metrics})
    with pytest.raises(OutputRejected, match="^summary_count$"):
        _verify(fixture)


@pytest.mark.parametrize(
    "collection,rows,reason",
    [
        ("all", [], "empty_pet_set"),
        ("unique", [], "empty_pet_set"),
        ("all", [ROWS[2]], "no_cis_denominator"),
    ],
)
def test_policy_a_rejects_whole_batch_with_sample_collection(
    tmp_path, collection, rows, reason
):
    fixture = _fixture(tmp_path, samples=("s000001", "s000002"))
    _bedpe(_file(fixture[0], collection, "s000002"), rows)
    before = _fingerprints(fixture[0])
    with pytest.raises(OutputRejected) as caught:
        _verify(fixture)
    assert (caught.value.reason_code, caught.value.sample, caught.value.collection) == (
        reason,
        "s000002",
        "noBg" if collection == "unique" else "all",
    )
    assert _fingerprints(fixture[0]) == before


@pytest.mark.parametrize(
    "case", ["duplicate", "wrong_first", "missing", "altered_strand"]
)
def test_no_bg_exact_coordinate_dedup_and_first_row(tmp_path, case):
    fixture = _fixture(tmp_path)
    rows = [ROWS[i] for i in (0, 2, 3)]
    if case == "duplicate":
        rows.append(ROWS[0])
    elif case == "wrong_first":
        rows[0] = ROWS[1]
    elif case == "missing":
        rows.pop()
    else:
        rows[0] = (*rows[0][:8], "-", "+")
    _bedpe(_file(fixture[0], "unique"), rows)
    with pytest.raises(OutputRejected) as caught:
        _verify(fixture)
    assert caught.value.reason_code == (
        "noBg_missing_pet" if case == "missing" else "noBg_content_mismatch"
    )


def test_background_float_midpoint_and_qc_truncated_midpoint_are_distinct(tmp_path):
    fixture = _fixture(tmp_path)
    root = fixture[0]
    # First PET: exact 999.5 bp is removed; second: 1000.5 retained, QC floors
    # each midpoint and reports 1000 (close). Other rows retain cis/trans support.
    rows = [
        ("chrA", "10", "21", "chrA", "1010", "1020", "1_-1_-1", "42", "+", "-"),
        ("chrA", "10", "21", "chrA", "1011", "1021", "2_-1_-1", "42", "+", "-"),
        ROWS[2],
        ROWS[3],
    ]
    # Make second distance 1000.5, with QC integer centers 15 and 1016 => 1001.
    # Use centers 15.0 and 1015.5 instead: float 1000.5 but integer 1000.
    rows[1] = ("chrA", "10", "20", "chrA", "1010", "1021", "2_-1_-1", "42", "+", "-")
    for mate in (1, 2):
        with gzip.open(root / "s000001" / f"s000001_R{mate}.fastq.gz", "wt") as stream:
            for row in rows:
                stream.write(f"@{row[6]}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    _bedpe(_file(root, "all"), rows)
    _bedpe(_file(root, "unique"), rows[1:])
    metrics = [4, 4, 0.75, 4, 0, 0.75, 1, 0, 0, 3, 0.75, 2 / 3, 1, 0, 0]
    _summary(root, {"s000001": metrics})
    assert _verify(fixture)["samples"]["s000001"]["metrics"] == metrics


@pytest.mark.parametrize(
    "position,value,reason",
    [
        (0, "unknown", "bedpe_unknown_contig"),
        (1, "-1", "bedpe_coordinate"),
        (2, "100001", "bedpe_reference_bounds"),
        (8, "?", "bedpe_strand"),
        (7, "256", "bedpe_mapq"),
        (6, "private/path", "invalid_internal_read_id"),
    ],
)
def test_bedpe_structure(tmp_path, position, value, reason):
    path = tmp_path / "a.gz"
    row = list(ROWS[0])
    row[position] = value
    _bedpe(path, [row])
    with pytest.raises(OutputRejected, match=f"^{reason}$"):
        list(iter_bedpe(path, CONTIGS))


def test_bad_gzip_and_symlink_are_rejected(tmp_path):
    broken = tmp_path / "broken.gz"
    broken.write_bytes(b"not gzip private-marker")
    with pytest.raises(OutputRejected, match="^bedpe_unreadable$"):
        list(iter_bedpe(broken, CONTIGS))
    linked = tmp_path / "linked.gz"
    linked.symlink_to(broken)
    with pytest.raises(OutputRejected, match="^output_not_regular$"):
        list(iter_bedpe(linked, CONTIGS))


def test_trimmed_fastq_pair_mismatch_and_unknown_bedpe_read(tmp_path):
    fixture = _fixture(tmp_path)
    root = fixture[0]
    path = root / "s000001/s000001_R2.fastq.gz"
    with gzip.open(path, "rt") as stream:
        data = stream.read()
    with gzip.open(path, "wt") as stream:
        stream.write(data.replace("@1_10_-1", "@9_10_-1"))
    with pytest.raises(OutputRejected, match="^trimmed_mate_id$"):
        _verify(fixture)
    with gzip.open(path, "wt") as stream:
        stream.write(data)
    row = (*ROWS[0][:6], "5_10_-1", *ROWS[0][7:])
    _bedpe(_file(root, "all"), [row])
    with pytest.raises(OutputRejected, match="^bedpe_unknown_trimmed_read$"):
        _verify(fixture)


def test_missing_and_duplicate_mapping_blocks_reject(tmp_path):
    fixture = _fixture(tmp_path)
    path = fixture[0] / "2000-01-01_tracPre2.py.log"
    original = path.read_text()
    path.write_text(original.replace("FLAG_A\n", ""))
    with pytest.raises(OutputRejected, match="^mapping_log_incomplete$"):
        _verify(fixture)
    path.write_text(original + original)
    with pytest.raises(OutputRejected, match="^mapping_log_sample$"):
        _verify(fixture)


def test_output_below_actual_mapq_is_rejected(tmp_path):
    fixture = _fixture(tmp_path, mapq=30)
    low = (*ROWS[0][:7], "29", *ROWS[0][8:])
    _bedpe(_file(fixture[0], "all"), [low, *ROWS[1:]])
    with pytest.raises(OutputRejected, match="^bedpe_below_mapq$"):
        _verify(fixture, 30)


def test_no_bg_can_have_no_cis_while_all_contains_cis(tmp_path):
    fixture = _fixture(tmp_path)
    root = fixture[0]
    rows = [
        ("chrA", "10", "20", "chrA", "510", "520", "1_-1_-1", "42", "+", "-"),
        ROWS[2],
    ]
    for mate in (1, 2):
        with gzip.open(root / "s000001" / f"s000001_R{mate}.fastq.gz", "wt") as stream:
            for row in [rows[0], ROWS[1], ROWS[2], ROWS[3]]:
                stream.write(f"@{row[6]}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    _bedpe(_file(root, "all"), rows)
    _bedpe(_file(root, "unique"), [ROWS[2]])
    with pytest.raises(OutputRejected) as caught:
        _verify(fixture)
    assert (caught.value.reason_code, caught.value.sample, caught.value.collection) == (
        "no_cis_denominator",
        "s000001",
        "noBg",
    )
