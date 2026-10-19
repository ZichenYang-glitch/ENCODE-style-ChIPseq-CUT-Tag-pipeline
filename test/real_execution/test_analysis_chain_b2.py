"""Real-tool contracts for B2c organelle filtering and mito fraction."""

import os
import shutil
import subprocess

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:1000",
    "@SQ\tSN:chrM\tLN:100",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _samtools():
    candidate = resolve_tool("samtools", "SAMTOOLS")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools(["samtools"], "B2c analysis-chain real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    return _samtools()


def _run(command, *, env=None):
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-1000:]}\n"
        f"stderr:\n{result.stderr[-1000:]}"
    )
    return result


def _record(name, rname, pos):
    return "\t".join(
        [
            name,
            "0",
            rname,
            str(pos),
            "60",
            "50M",
            "*",
            "0",
            "0",
            SEQ,
            QUAL,
        ]
    )


def _make_bam(tmp_path, samtools_path):
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "input.bam"
    sam.write_text(
        "\n".join(
            HEADER
            + [
                _record("r1", "chr1", 100),
                _record("r2", "chr1", 200),
                _record("m1", "chrM", 10),
                _record("m2", "chrM", 50),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    return bam


def test_organelle_filter_drops_chrm(tmp_path, samtools_path):
    bam = _make_bam(tmp_path, samtools_path)
    out = tmp_path / "analysis.bam"
    script = (
        'samtools view -h "$1" '
        '| awk -v omit="chrM MT chrMT" '
        '\'BEGIN { n=split(omit,a," "); for (i=1;i<=n;i++) drop[a[i]]=1 } '
        "/^@/ { print; next } !($3 in drop) { print }' "
        '| samtools view -b -o "$2" -'
    )
    _run(["bash", "-c", script, "bash", str(bam), str(out)])
    records = _run([samtools_path, "view", str(out)]).stdout.splitlines()
    assert len(records) == 2
    assert {line.split("\t")[2] for line in records} == {"chr1"}


def test_mito_fraction_matches_idxstats(tmp_path, samtools_path):
    bam = _make_bam(tmp_path, samtools_path)
    out = tmp_path / "mito.tsv"
    script = """
samtools idxstats "$1" \
| awk -v omit="chrM MT chrMT" -v sample="S1" -v mode="exclude" '
BEGIN { n=split(omit,a," "); for (i=1;i<=n;i++) drop[a[i]]=1 }
{ total += $3; if ($1 in drop) organelle += $3 }
END {
    printf "sample\\torganelle_reads\\ttotal_mapped_reads\\tfraction\\tmode\\n";
    frac = (total > 0) ? organelle / total : 0;
    printf "%s\\t%d\\t%d\\t%.6f\\t%s\\n", sample, organelle, total, frac, mode
}' > "$2"
"""
    _run(["bash", "-c", script, "bash", str(bam), str(out)])
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "sample\torganelle_reads\ttotal_mapped_reads\tfraction\tmode"
    fields = lines[1].split("\t")
    assert fields[0] == "S1"
    assert int(fields[1]) == 2
    assert int(fields[2]) == 4
    assert abs(float(fields[3]) - 0.5) < 1e-6
    assert fields[4] == "exclude"
