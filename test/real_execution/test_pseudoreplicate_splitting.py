"""Real-tool contracts for deterministic pseudoreplicate splitting."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

TEST_ROOT = Path(__file__).resolve().parents[1]


REPO_ROOT = TEST_ROOT.parent
SPLIT_SCRIPT = REPO_ROOT / "scripts" / "split_pseudoreps.py"

pytestmark = pytest.mark.real_execution


@pytest.fixture(scope="module")
def samtools_path():
    candidate = resolve_tool("samtools", "SAMTOOLS")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            resolved = candidate
    if not resolved:
        require_external_tools(
            ["samtools"],
            "pseudoreplicate real execution",
        )
    return resolved


def _run(command, samtools_path):
    env = os.environ.copy()
    env["SAMTOOLS"] = samtools_path
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-500:]}\n"
        f"stderr:\n{result.stderr[-500:]}"
    )
    return result


def _read_names(bam_path, samtools_path):
    result = _run([samtools_path, "view", bam_path], samtools_path)
    return {line.split("\t", 1)[0] for line in result.stdout.splitlines() if line}


def test_split_is_complementary_mate_preserving_indexed_and_deterministic(
    tmp_path, samtools_path
):
    sam_path = tmp_path / "input.sam"
    input_bam = tmp_path / "input.bam"
    first_pr1 = tmp_path / "first.pr1.bam"
    first_pr2 = tmp_path / "first.pr2.bam"
    second_pr1 = tmp_path / "second.pr1.bam"
    second_pr2 = tmp_path / "second.pr2.bam"

    sam_path.write_text(
        "\n".join(
            [
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:1000",
                "read1/1\t99\tchr1\t1\t60\t100M\t=\t100\t99\t*\t*",
                "read1/2\t147\tchr1\t100\t60\t100M\t=\t1\t-99\t*\t*",
                "read2\t0\tchr1\t500\t60\t100M\t*\t0\t0\t*\t*",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _run(
        [samtools_path, "view", "-bS", "-o", input_bam, sam_path],
        samtools_path,
    )

    for out1, out2 in (
        (first_pr1, first_pr2),
        (second_pr1, second_pr2),
    ):
        _run(
            [
                sys.executable,
                SPLIT_SCRIPT,
                "--input",
                input_bam,
                "--out1",
                out1,
                "--out2",
                out2,
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            samtools_path,
        )

    assert first_pr1.is_file()
    assert first_pr2.is_file()
    assert first_pr1.with_suffix(first_pr1.suffix + ".bai").stat().st_size > 0
    assert first_pr2.with_suffix(first_pr2.suffix + ".bai").stat().st_size > 0

    input_count = int(
        _run([samtools_path, "view", "-c", input_bam], samtools_path).stdout
    )
    first_count = int(
        _run([samtools_path, "view", "-c", first_pr1], samtools_path).stdout
    )
    second_count = int(
        _run([samtools_path, "view", "-c", first_pr2], samtools_path).stdout
    )
    assert first_count + second_count == input_count

    pr1_names = _read_names(first_pr1, samtools_path)
    pr2_names = _read_names(first_pr2, samtools_path)
    assert pr1_names.isdisjoint(pr2_names)
    assert ({"read1/1", "read1/2"} <= pr1_names) or (
        {"read1/1", "read1/2"} <= pr2_names
    )

    first_view = _run([samtools_path, "view", first_pr1], samtools_path).stdout
    repeated_view = _run([samtools_path, "view", second_pr1], samtools_path).stdout
    assert repeated_view == first_view
