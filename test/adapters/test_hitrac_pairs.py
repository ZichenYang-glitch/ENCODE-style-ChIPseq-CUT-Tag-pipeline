"""Independent PET provenance, without a converter or a scientific tool fixture."""

from contextlib import closing
from pathlib import Path
import sqlite3

import pytest

from encode_pipeline.adapters.hitrac_preprocess.pairs import PairError, verify_pairs


def sam(qname, flag, start, *, contig="chr1", cigar="20M"):
    return f"{qname}\t{flag}\t{contig}\t{start + 1}\t0\t{cigar}\t*\t0\t0\t*\t*\n"


def pet(qname, first=("chr1", 10, 30, "+"), second=("chr1", 100, 120, "-")):
    return tuple(
        str(value) for value in (*first[:3], *second[:3], qname, 0, first[3], second[3])
    )


def verify(tmp_path, lines, rows):
    verify_pairs(iter(lines), iter(rows), tmp_path / "pairs.sqlite")


@pytest.mark.parametrize("flags", [(99, 147), (97, 145), (65, 145), (321, 2193)])
def test_roles_not_proper_pair_mapq_or_secondary_status_control_support(
    tmp_path, flags
):
    verify(tmp_path, [sam("q", flags[0], 10), sam("q", flags[1], 100)], [pet("q")])


def test_whole_endpoint_swap_and_trans(tmp_path):
    left = ("chr2", 100, 120, "-")
    right = ("chr1", 10, 30, "+")
    verify(
        tmp_path,
        [sam("q", 65, 10), sam("q", 145, 100, contig="chr2")],
        [pet("q", left, right)],
    )


def test_same_position_and_strand_still_needs_two_mate_records(tmp_path):
    endpoint = ("chr1", 10, 30, "+")
    verify(
        tmp_path, [sam("q", 65, 10), sam("q", 129, 10)], [pet("q", endpoint, endpoint)]
    )


def test_reference_consuming_cigar_and_soft_clipping(tmp_path):
    verify(
        tmp_path,
        [sam("q", 65, 10, cigar="5S10M2I5M3D4N2=1X5H"), sam("q", 145, 100)],
        [pet("q", ("chr1", 10, 35, "+"))],
    )


def test_natural_and_noncontiguous_qname_order_and_duplicate_endpoints(tmp_path):
    lines = [sam(q, 65, 10) for q in ("q1", "q2", "q10")]
    lines += [sam(q, 145, 100) for q in ("q10", "q1", "q2")]
    verify(tmp_path, lines, [pet(q) for q in ("q10", "q2", "q1")])
    with closing(sqlite3.connect(tmp_path / "pairs.sqlite")) as database:
        assert database.execute("SELECT COUNT(*) FROM alignments").fetchone()[0] == 6


def test_multi_alignment_search_is_not_first_record_and_ignores_unrelated_orphans(
    tmp_path,
):
    lines = [
        sam("q", 65, 800),
        sam("q", 145, 900),
        sam("q", 321, 10),
        sam("q", 2193, 100),
        sam("orphan", 73, 400),
        sam("ambiguous", 193, 450),
        sam("neither", 1, 500),
        "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",
    ]
    verify(tmp_path, lines, [pet("q")])


def test_duplicate_rows_are_checked_without_claiming_unique_consumption(tmp_path):
    observed = []

    def rows():
        for index in range(3):
            observed.append(index)
            yield pet("q")

    verify_pairs(
        [sam("q", 65, 10), sam("q", 145, 100)], rows(), tmp_path / "pairs.sqlite"
    )
    assert observed == [0, 1, 2]


@pytest.mark.parametrize(
    "lines, row",
    [
        ([sam("q", 65, 10)], pet("q", ("chr1", 10, 30, "+"), ("chr1", 10, 30, "+"))),
        ([sam("q", 193, 10)], pet("q", ("chr1", 10, 30, "+"), ("chr1", 10, 30, "+"))),
        ([sam("q", 1, 10), sam("q", 145, 100)], pet("q")),
        ([sam("q1", 65, 10), sam("q2", 145, 100)], pet("q1")),
        ([sam("q", 65, 10), sam("q", 145, 100)], pet("q", ("chr1", 10, 31, "+"))),
        ([sam("q", 65, 10), sam("q", 145, 100)], pet("q", ("chr1", 10, 30, "-"))),
        (
            [sam("q", 65, 10), sam("q", 145, 100)],
            pet("q", ("chr1", 100, 120, "+"), ("chr1", 10, 30, "-")),
        ),
    ],
    ids=[
        "single_record",
        "both_role_bits",
        "missing_role",
        "cross_qname",
        "wrong_end",
        "wrong_strand",
        "independently_sorted_ends",
    ],
)
def test_unsupported_rows_are_rejected_without_raw_details(tmp_path, lines, row):
    with pytest.raises(PairError, match="^pair_support_missing$") as error:
        verify(tmp_path, lines, [row])
    assert error.value.reason_code == "pair_support_missing"


def test_samples_cannot_reuse_one_anothers_records(tmp_path):
    verify_pairs(
        [sam("q", 65, 10), sam("q", 145, 100)], [pet("q")], tmp_path / "sample1.sqlite"
    )
    with pytest.raises(PairError, match="^pair_support_missing$"):
        verify_pairs([sam("q", 145, 100)], [pet("q")], tmp_path / "sample2.sqlite")


@pytest.mark.parametrize("cigar", ["*", "0M", "20Mbad", "5S", "20Z"])
def test_invalid_cigar_is_rejected(tmp_path, cigar):
    with pytest.raises(PairError, match="^pair_sam_invalid$"):
        verify(tmp_path, [sam("q", 65, 10, cigar=cigar)], [])


def test_old_database_cannot_supply_evidence(tmp_path):
    database = tmp_path / "pairs.sqlite"
    verify_pairs([sam("q", 65, 10), sam("q", 145, 100)], [pet("q")], database)
    original = database.read_bytes()
    with pytest.raises(PairError, match="^pair_index_unavailable$"):
        verify_pairs([], [pet("q")], database)
    assert database.read_bytes() == original


def test_failure_keeps_committed_evidence_and_releases_database(tmp_path):
    database = tmp_path / "pairs.sqlite"
    with pytest.raises(PairError, match="^pair_support_missing$"):
        verify_pairs([sam("q", 65, 10)], [pet("q")], database)
    assert database.stat().st_mode & 0o777 == 0o600
    with closing(sqlite3.connect(database, timeout=0)) as connection:
        assert connection.execute("SELECT COUNT(*) FROM alignments").fetchone()[0] == 1
        connection.execute("BEGIN EXCLUSIVE")
    assert not Path(str(database) + "-journal").exists()


def test_large_stream_does_not_require_sorted_qname_groups(tmp_path):
    def records():
        yield "@HD\tVN:1.6\tSO:queryname\n"
        for role, start in ((65, 10), (145, 100)):
            for index in range(2049):
                yield sam(f"q{index}", role, start)

    verify_pairs(
        records(),
        (pet(f"q{index}") for index in range(2049)),
        tmp_path / "pairs.sqlite",
    )
