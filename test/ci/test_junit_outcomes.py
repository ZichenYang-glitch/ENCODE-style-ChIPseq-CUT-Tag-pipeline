"""Contracts for the CI JUnit outcome gate."""

from pathlib import Path

from scripts.check_junit_outcomes import main, summarize


def _write_report(path: Path, testcase_xml: str) -> None:
    path.write_text(
        f'<testsuites><testsuite name="pytest">{testcase_xml}</testsuite></testsuites>',
        encoding="utf-8",
    )


def test_accepts_a_nonempty_all_passing_report(tmp_path):
    report = tmp_path / "report.xml"
    summary = tmp_path / "summary.md"
    _write_report(
        report,
        '<testcase classname="test_example" name="test_ok" time="1.25" />',
    )

    assert (
        main([str(report), "--label", "PR fast", "--summary-file", str(summary)]) == 0
    )
    assert summarize(report) == {
        "collected": 1,
        "passed": 1,
        "failures": 0,
        "errors": 0,
        "skipped": 0,
        "xfailed": 0,
        "seconds": 1.25,
    }
    assert "PR fast pytest outcomes" in summary.read_text(encoding="utf-8")


def test_rejects_skip_and_xfail_outcomes(tmp_path):
    report = tmp_path / "report.xml"
    _write_report(
        report,
        """
        <testcase name="test_skip"><skipped type="pytest.skip" /></testcase>
        <testcase name="test_xfail"><skipped type="pytest.xfail" /></testcase>
        """,
    )

    assert main([str(report)]) == 1
    counts = summarize(report)
    assert counts["skipped"] == 2
    assert counts["xfailed"] == 1


def test_rejects_failures_errors_empty_and_malformed_reports(tmp_path):
    failing = tmp_path / "failing.xml"
    empty = tmp_path / "empty.xml"
    malformed = tmp_path / "malformed.xml"
    _write_report(
        failing,
        """
        <testcase name="test_failure"><failure /></testcase>
        <testcase name="test_error"><error /></testcase>
        """,
    )
    _write_report(empty, "")
    malformed.write_text("<testsuite>", encoding="utf-8")

    assert main([str(failing)]) == 1
    assert main([str(empty)]) == 1
    assert main([str(malformed)]) == 1
