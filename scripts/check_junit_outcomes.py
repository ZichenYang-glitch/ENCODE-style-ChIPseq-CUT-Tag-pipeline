#!/usr/bin/env python3
"""Fail a CI tier when its pytest JUnit report contains non-passing outcomes."""

from __future__ import annotations

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def summarize(report: Path) -> dict[str, int | float]:
    """Return outcome counts derived from testcase nodes in a JUnit report."""
    root = ET.parse(report).getroot()
    cases = root.findall(".//testcase")
    failures = sum(case.find("failure") is not None for case in cases)
    errors = sum(case.find("error") is not None for case in cases)
    skipped_nodes = [case.find("skipped") for case in cases]
    skipped_nodes = [node for node in skipped_nodes if node is not None]
    xfailed = sum(
        "xfail" in f"{node.get('type', '')} {node.get('message', '')}".lower()
        for node in skipped_nodes
    )
    skipped = len(skipped_nodes)
    passed = len(cases) - failures - errors - skipped
    duration = sum(float(case.get("time", "0") or 0) for case in cases)
    return {
        "collected": len(cases),
        "passed": passed,
        "failures": failures,
        "errors": errors,
        "skipped": skipped,
        "xfailed": xfailed,
        "seconds": duration,
    }


def _summary_markdown(label: str, counts: dict[str, int | float]) -> str:
    return (
        f"### {label} pytest outcomes\n\n"
        "| Collected | Passed | Failed | Errors | Skipped | Xfailed | Test time |\n"
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n"
        f"| {counts['collected']} | {counts['passed']} | {counts['failures']} | "
        f"{counts['errors']} | {counts['skipped']} | {counts['xfailed']} | "
        f"{counts['seconds']:.2f}s |\n"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--label", default="Python tier")
    parser.add_argument(
        "--summary-file",
        type=Path,
        default=os.environ.get("GITHUB_STEP_SUMMARY"),
    )
    args = parser.parse_args(argv)

    try:
        counts = summarize(args.report)
    except (OSError, ET.ParseError, ValueError) as error:
        print(f"Unable to read JUnit report {args.report}: {error}", file=sys.stderr)
        return 1

    markdown = _summary_markdown(args.label, counts)
    print(markdown, end="")
    if args.summary_file:
        with args.summary_file.open("a", encoding="utf-8") as handle:
            handle.write(markdown)

    invalid = (
        counts["collected"] == 0
        or counts["failures"] > 0
        or counts["errors"] > 0
        or counts["skipped"] > 0
    )
    if invalid:
        print(
            "JUnit gate requires at least one test and zero failures, errors, skips, "
            "or xfails.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
