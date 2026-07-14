"""Behavior tests for the replicate-consensus peak engine."""

from __future__ import annotations

import csv
import json

import pytest

import scripts.compute_consensus as compute_consensus


SUMMARY_COLUMNS = [
    "experiment",
    "assay",
    "peak_mode",
    "caller",
    "n_bioreps",
    "min_replicates",
    "reciprocal_overlap",
    "consensus_peak_count",
    "support_distribution",
    "biorep_labels",
    "source_peak_files",
    "final_method",
    "final_output",
]


def _narrowpeak_row(
    chrom,
    start,
    end,
    name=".",
    score=100,
    strand=".",
    signal=5.0,
    pvalue="3.0",
    qvalue="2.0",
    summit=100,
):
    return "\t".join(
        str(value)
        for value in (
            chrom,
            start,
            end,
            name,
            score,
            strand,
            signal,
            pvalue,
            qvalue,
            summit,
        )
    )


def _broadpeak_row(
    chrom,
    start,
    end,
    name=".",
    score=100,
    strand=".",
    signal=5.0,
    pvalue="3.0",
    qvalue="2.0",
):
    return "\t".join(
        str(value)
        for value in (
            chrom,
            start,
            end,
            name,
            score,
            strand,
            signal,
            pvalue,
            qvalue,
        )
    )


def _bed_row(chrom, start, end, name=".", score=100, strand="."):
    return "\t".join(str(value) for value in (chrom, start, end, name, score, strand))


def _bed3_row(chrom, start, end):
    return "\t".join(str(value) for value in (chrom, start, end))


def _nonempty_lines(path):
    return [line for line in path.read_text(encoding="utf-8").splitlines() if line]


def _summary(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == SUMMARY_COLUMNS
    assert len(rows) == 1
    return rows[0]


@pytest.fixture
def consensus_case(tmp_path):
    """Build CLI arguments and disposable peak files for one engine invocation."""

    invocation = 0

    def _make(
        peak_rows,
        *,
        bioreps=None,
        peak_format="narrowPeak",
        min_replicates=2,
        reciprocal_overlap=0.5,
        experiment="test_exp",
        assay="chipseq",
        caller="macs3",
        peak_mode="narrow",
        final_method="",
        final_output="",
    ):
        nonlocal invocation
        invocation += 1
        input_dir = tmp_path / f"inputs_{invocation}"
        input_dir.mkdir()

        suffix = {
            "narrowPeak": ".narrowPeak",
            "broadPeak": ".broadPeak",
            "bed": ".bed",
        }.get(peak_format, ".peaks")
        peak_paths = []
        for index, rows in enumerate(peak_rows, start=1):
            path = input_dir / f"biorep_{index}{suffix}"
            path.write_text("\n".join(rows) + "\n", encoding="utf-8")
            peak_paths.append(path)

        if bioreps is None:
            bioreps = [str(index) for index in range(1, len(peak_paths) + 1)]

        output = tmp_path / f"outputs_{invocation}" / f"consensus{suffix}"
        summary = tmp_path / f"outputs_{invocation}" / "consensus.summary.tsv"
        argv = [
            "--peaks",
            *(str(path) for path in peak_paths),
            "--bioreps",
            *bioreps,
            "--format",
            peak_format,
            "--min-replicates",
            str(min_replicates),
            "--reciprocal-overlap",
            str(reciprocal_overlap),
            "--output",
            str(output),
            "--summary",
            str(summary),
            "--experiment",
            experiment,
            "--assay",
            assay,
            "--caller",
            caller,
            "--peak-mode",
            peak_mode,
            "--final-method",
            final_method,
            "--final-output",
            final_output,
        ]
        return argv, output, summary, peak_paths

    return _make


@pytest.mark.parametrize(
    ("support", "total", "expected_score"),
    [(2, 2, "1000"), (2, 3, "667"), (3, 3, "1000"), (3, 4, "750")],
)
def test_support_fraction_sets_consensus_score(
    consensus_case, support, total, expected_score
):
    peak_rows = []
    for index in range(total):
        if index < support:
            peak_rows.append([_narrowpeak_row("chr1", 100, 500, summit=200)])
        else:
            peak_rows.append([_narrowpeak_row(f"chr{index + 2}", 600, 1000)])

    argv, output, _, _ = consensus_case(peak_rows)
    assert compute_consensus.main(argv) == 0

    lines = _nonempty_lines(output)
    assert len(lines) == 1
    assert lines[0].split("\t")[4] == expected_score


@pytest.mark.parametrize(
    ("left", "right", "expected_count"),
    [
        ((100, 600), (550, 1050), 0),
        ((100, 600), (200, 700), 1),
        ((100, 300), (200, 400), 1),
    ],
)
def test_reciprocal_overlap_threshold_controls_clustering(
    consensus_case, left, right, expected_count
):
    argv, output, _, _ = consensus_case(
        [
            [_narrowpeak_row("chr1", *left)],
            [_narrowpeak_row("chr1", *right)],
        ]
    )
    assert compute_consensus.main(argv) == 0
    assert len(_nonempty_lines(output)) == expected_count


def test_multiple_peaks_from_one_biorep_count_as_one_support(consensus_case):
    argv, output, _, _ = consensus_case(
        [
            [
                _narrowpeak_row("chr1", 100, 500, signal=5.0, summit=200),
                _narrowpeak_row("chr1", 120, 480, signal=4.0, summit=150),
            ],
            [_narrowpeak_row("chr1", 100, 500, signal=3.0, summit=200)],
        ]
    )
    assert compute_consensus.main(argv) == 0

    lines = _nonempty_lines(output)
    assert len(lines) == 1
    assert lines[0].split("\t")[4] == "1000"


@pytest.mark.parametrize(
    ("peak_rows", "expected_interval"),
    [
        (
            [
                [_narrowpeak_row("chr1", 100, 300)],
                [_narrowpeak_row("chr1", 200, 400)],
                [_narrowpeak_row("chr1", 300, 500)],
            ],
            ("chr1", "100", "500"),
        ),
        (
            [
                [_narrowpeak_row("chr1", 100, 200)],
                [_narrowpeak_row("chr1", 200, 300)],
            ],
            None,
        ),
        (
            [
                [_narrowpeak_row("chr1", 100, 200)],
                [_narrowpeak_row("chr2", 100, 200)],
            ],
            None,
        ),
    ],
    ids=["chained", "adjacent", "different-chromosomes"],
)
def test_overlap_components_respect_interval_topology(
    consensus_case, peak_rows, expected_interval
):
    argv, output, _, _ = consensus_case(peak_rows)
    assert compute_consensus.main(argv) == 0

    lines = _nonempty_lines(output)
    if expected_interval is None:
        assert lines == []
    else:
        assert len(lines) == 1
        assert tuple(lines[0].split("\t")[:3]) == expected_interval
        assert lines[0].split("\t")[4] == "1000"


def test_consensus_output_is_coordinate_sorted(consensus_case):
    rows = [
        _narrowpeak_row("chr2", 100, 300),
        _narrowpeak_row("chr1", 500, 700),
    ]
    argv, output, _, _ = consensus_case([rows, rows])
    assert compute_consensus.main(argv) == 0
    assert [line.split("\t")[0] for line in _nonempty_lines(output)] == [
        "chr1",
        "chr2",
    ]


@pytest.mark.parametrize(
    ("peak_count", "bioreps", "min_replicates", "error_fragment"),
    [
        (2, ["1"], 2, "peak file(s) but"),
        (2, ["1", "1"], 2, "duplicate biorep label"),
        (2, ["1", "2"], 1, "min-replicates must be >= 2"),
        (2, ["1", "2"], 3, "exceeds number of biological replicates"),
        (1, ["1"], 2, "at least two peak files"),
    ],
    ids=[
        "count-mismatch",
        "duplicate-label",
        "minimum-too-low",
        "minimum-too-high",
        "one-file",
    ],
)
def test_invalid_replicate_configuration_exits(
    consensus_case,
    capsys,
    peak_count,
    bioreps,
    min_replicates,
    error_fragment,
):
    peak_rows = [[_narrowpeak_row("chr1", 100, 200)] for _ in range(peak_count)]
    argv, _, _, _ = consensus_case(
        peak_rows,
        bioreps=bioreps,
        min_replicates=min_replicates,
    )

    with pytest.raises(SystemExit) as exc_info:
        compute_consensus.main(argv)

    assert exc_info.value.code == 1
    assert error_fragment in capsys.readouterr().err


@pytest.mark.parametrize(
    ("peak_format", "row_factory", "expected_columns", "expected_name"),
    [
        ("narrowPeak", _narrowpeak_row, 10, None),
        ("broadPeak", _broadpeak_row, 9, None),
        ("bed", _bed_row, 6, None),
        ("bed", _bed3_row, 6, "consensus_peak_1"),
    ],
    ids=["narrowpeak", "broadpeak", "bed6", "bed3-to-bed6"],
)
def test_output_format_contract(
    consensus_case, peak_format, row_factory, expected_columns, expected_name
):
    row = row_factory("chr1", 100, 500)
    peak_mode = "broad" if peak_format == "broadPeak" else "narrow"
    argv, output, _, _ = consensus_case(
        [[row], [row]], peak_format=peak_format, peak_mode=peak_mode
    )
    assert compute_consensus.main(argv) == 0

    fields = _nonempty_lines(output)[0].split("\t")
    assert len(fields) == expected_columns
    if peak_format == "narrowPeak":
        assert 0 <= int(fields[9]) < int(fields[2]) - int(fields[1])
    if expected_name is not None:
        assert fields[3] == expected_name


def test_empty_consensus_writes_empty_peak_and_zero_summary(consensus_case):
    argv, output, summary, _ = consensus_case(
        [
            [_narrowpeak_row("chr1", 100, 200)],
            [_narrowpeak_row("chr2", 100, 200)],
        ]
    )
    assert compute_consensus.main(argv) == 0
    assert output.read_text(encoding="utf-8") == ""

    row = _summary(summary)
    assert row["consensus_peak_count"] == "0"
    assert json.loads(row["support_distribution"]) == {}


def test_malformed_peak_reports_physical_line_number(consensus_case, capsys):
    argv, _, _, _ = consensus_case(
        [
            ["# header", "", _narrowpeak_row("chr1", "bad_start", 200)],
            [_narrowpeak_row("chr1", 100, 200)],
        ]
    )

    with pytest.raises(SystemExit) as exc_info:
        compute_consensus.main(argv)

    assert exc_info.value.code == 1
    assert "line 3" in capsys.readouterr().err


def test_peak_end_must_exceed_start(consensus_case, capsys):
    argv, _, _, _ = consensus_case(
        [
            [_narrowpeak_row("chr1", 200, 100)],
            [_narrowpeak_row("chr1", 100, 200)],
        ]
    )

    with pytest.raises(SystemExit) as exc_info:
        compute_consensus.main(argv)

    assert exc_info.value.code == 1
    assert "must be greater than start" in capsys.readouterr().err


@pytest.mark.parametrize("invalid_summit", [9999, "not_an_int"])
def test_invalid_source_summit_uses_midpoint(consensus_case, invalid_summit):
    argv, output, _, _ = consensus_case(
        [
            [_narrowpeak_row("chr1", 100, 200, signal=5.0, summit=invalid_summit)],
            [_narrowpeak_row("chr1", 100, 200, signal=3.0, summit=50)],
        ]
    )
    assert compute_consensus.main(argv) == 0
    assert _nonempty_lines(output)[0].split("\t")[9] == "50"


def test_summary_schema_and_values(consensus_case):
    final_output = "../06_idr/final/conservative.narrowPeak"
    argv, _, summary, peak_paths = consensus_case(
        [
            [_narrowpeak_row("chr1", 100, 200)],
            [_narrowpeak_row("chr1", 100, 200)],
        ],
        bioreps=["B", "A"],
        experiment="myexp",
        assay="cuttag",
        caller="macs3",
        peak_mode="narrow",
        final_method="idr",
        final_output=final_output,
    )
    assert compute_consensus.main(argv) == 0

    row = _summary(summary)
    assert row == {
        "experiment": "myexp",
        "assay": "cuttag",
        "peak_mode": "narrow",
        "caller": "macs3",
        "n_bioreps": "2",
        "min_replicates": "2",
        "reciprocal_overlap": "0.5",
        "consensus_peak_count": "1",
        "support_distribution": '{"2": 1}',
        "biorep_labels": "A,B",
        "source_peak_files": json.dumps([str(path.resolve()) for path in peak_paths]),
        "final_method": "idr",
        "final_output": final_output,
    }
    assert json.loads(row["support_distribution"]) == {"2": 1}
    assert json.loads(row["source_peak_files"]) == [
        str(path.resolve()) for path in peak_paths
    ]


@pytest.mark.parametrize("reciprocal_overlap", [0, 1.5])
def test_reciprocal_overlap_must_be_in_valid_range(
    consensus_case, capsys, reciprocal_overlap
):
    row = _narrowpeak_row("chr1", 100, 200)
    argv, _, _, _ = consensus_case(
        [[row], [row]], reciprocal_overlap=reciprocal_overlap
    )

    with pytest.raises(SystemExit) as exc_info:
        compute_consensus.main(argv)

    assert exc_info.value.code == 1
    assert "reciprocal-overlap must be in (0, 1]" in capsys.readouterr().err


def test_argparse_rejects_unknown_peak_format(consensus_case, capsys):
    row = _narrowpeak_row("chr1", 100, 200)
    argv, _, _, _ = consensus_case([[row], [row]], peak_format="gff")

    with pytest.raises(SystemExit) as exc_info:
        compute_consensus.main(argv)

    assert exc_info.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_output_parent_directories_are_created(consensus_case, tmp_path):
    row = _narrowpeak_row("chr1", 100, 200)
    argv, _, _, _ = consensus_case([[row], [row]])
    output = tmp_path / "nested" / "peaks" / "consensus.narrowPeak"
    summary = tmp_path / "nested" / "summary" / "consensus.tsv"
    argv[argv.index("--output") + 1] = str(output)
    argv[argv.index("--summary") + 1] = str(summary)

    assert not output.parent.exists()
    assert not summary.parent.exists()
    assert compute_consensus.main(argv) == 0
    assert output.is_file()
    assert summary.is_file()
