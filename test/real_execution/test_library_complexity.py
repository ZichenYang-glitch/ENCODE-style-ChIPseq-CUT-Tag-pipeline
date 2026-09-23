"""Real samtools/preseq and original-rule complexity contracts on tiny BAMs."""

import csv
import os
import re
import shutil
import sys

import pytest

from _complexity import REPO, project, run, snakemake
from _tool_resolver import resolve_tool


pytestmark = pytest.mark.real_execution


@pytest.fixture
def tools(tmp_path, monkeypatch):
    resolved = {}
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name in ("samtools", "preseq"):
        path = shutil.which(resolve_tool(name, name.upper()))
        assert path, f"Real complexity tests require {name}; no simulated substitute"
        resolved[name] = path
        (bindir / name).symlink_to(path)
    (bindir / "python3").symlink_to(sys.executable)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    return resolved


def alignment(name, pos, *, flag=0, mate=0, tlen=0, cigar="50M"):
    return [
        name,
        str(flag),
        "chr1",
        str(pos),
        "60",
        cigar,
        "=" if mate else "*",
        str(mate),
        str(tlen),
        "A" * 50,
        "I" * 50,
    ]


def fragment(name, pos, layout, *, length=150, duplicate=False):
    extra = 0x400 if duplicate else 0
    if layout == "SE":
        return [alignment(name, pos, flag=extra)]
    end = pos + length - 50
    return [
        alignment(name, pos, flag=99 | extra, mate=end, tlen=length),
        alignment(name, end, flag=147 | extra, mate=pos, tlen=-length),
    ]


def write_bam(tmp_path, tools, rows, name="input"):
    # Validate fixture records independently of the metric's fragment-key code.
    pairs = {}
    for row in rows:
        query_length = sum(
            int(n)
            for n, op in re.findall(r"(\d+)([MIDNSHP=X])", row[5])
            if op in "MIS=X"
        )
        assert query_length == len(row[9]) == len(row[10])
        if int(row[1]) & 1:
            pairs.setdefault(row[0], []).append(row)
    for first, second in pairs.values():
        assert int(first[1]) & 0x43 == 0x43
        assert int(second[1]) & 0x83 == 0x83
        assert first[2] == second[2] == "chr1"
        assert first[6] == second[6] == "="
        assert first[3] == second[7] and first[7] == second[3]
        assert int(first[8]) == -int(second[8]) == int(second[3]) + 50 - int(first[3])
    sam = tmp_path / f"{name}.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n"
        + "".join("\t".join(row) + "\n" for row in rows)
    )
    bam = tmp_path / f"{name}.bam"
    run(tmp_path, [tools["samtools"], "sort", "-o", bam, sam])
    run(tmp_path, [tools["samtools"], "quickcheck", bam])
    return bam


def dedup(tmp_path, tools, bam):
    name, fix, pos, final = [
        tmp_path / f"dedup-{suffix}.bam" for suffix in ("name", "fix", "pos", "final")
    ]
    for args in (
        ["sort", "-n", "-o", name, bam],
        ["fixmate", "-m", name, fix],
        ["sort", "-o", pos, fix],
        ["markdup", "-r", pos, final],
    ):
        run(tmp_path, [tools["samtools"], *args])
    return final


def metrics(tmp_path, bam, name):
    output = tmp_path / f"{name}.tsv"
    run(
        tmp_path,
        [
            sys.executable,
            REPO / "scripts/calc_nrf_pbc.py",
            "--sample",
            "S1",
            "--bam",
            bam,
            "--output",
            output,
        ],
    )
    with output.open() as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def assert_metrics(row, total, distinct, one, two):
    assert row == dict(
        sample="S1",
        total_fragments=str(total),
        distinct_fragments=str(distinct),
        one_read_fragments=str(one),
        two_read_fragments=str(two),
        nrf=f"{distinct / total:.6f}" if total else "NA",
        pbc1=f"{one / distinct:.6f}" if distinct else "NA",
        pbc2=f"{one / two:.6f}" if two else "NA",
    )


@pytest.mark.parametrize("layout", ["SE", "PE"])
def test_nrf_before_and_after_real_duplicate_removal(tmp_path, tools, layout):
    rows = [
        row
        for pos, count in ((100, 1), (500, 2), (900, 3))
        for copy in range(count)
        for row in fragment(f"r{pos}_{copy}", pos, layout, duplicate=copy > 0)
    ]
    bam = write_bam(tmp_path, tools, rows)
    flagged = run(tmp_path, [tools["samtools"], "view", "-c", "-f", "1024", bam])
    assert int(flagged.stdout) == (6 if layout == "PE" else 3)
    assert_metrics(metrics(tmp_path, bam, "before"), 6, 3, 1, 1)
    final = dedup(tmp_path, tools, bam)
    assert_metrics(metrics(tmp_path, final, "after"), 3, 3, 3, 0)


@pytest.mark.parametrize("case", ["empty", "unique", "clipping"])
def test_nrf_legal_boundaries_and_clipping_counterexample(tmp_path, tools, case):
    rows = [] if case == "empty" else [alignment("a", 100), alignment("b", 500)]
    if case == "clipping":
        rows = [alignment("a", 100), alignment("b", 100, cigar="5S45M")]
    bam = write_bam(tmp_path, tools, rows)
    if case == "clipping":
        bam = dedup(tmp_path, tools, bam)
    expected = {"empty": (0, 0, 0, 0), "unique": (2, 2, 2, 0), "clipping": (2, 1, 0, 1)}
    assert_metrics(metrics(tmp_path, bam, case), *expected[case])


def spectrum(layout):
    """658 loci, 1070 molecules; paired loci share mate1 but differ in mate2."""
    rows = []
    group = 0
    for count, loci in ((1, 400), (2, 160), (3, 60), (4, 24), (5, 10), (6, 4)):
        for _ in range(loci // 2):
            for alternative in (0, 1):
                pos = 100 + group * 500 + (alternative * 200 if layout == "SE" else 0)
                for copy in range(count):
                    rows.extend(
                        fragment(
                            f"g{group}a{alternative}c{copy}",
                            pos,
                            layout,
                            length=150 + alternative * 30,
                            duplicate=copy > 0,
                        )
                    )
            group += 1
    return rows


def test_preseq_pair_mode_is_the_only_difference(tmp_path, tools):
    bam = write_bam(tmp_path, tools, spectrum("PE"))
    results = []
    for paired in (False, True):
        output = tmp_path / f"paired-{paired}.txt"
        result = run(
            tmp_path,
            [
                tools["preseq"],
                "lc_extrap",
                "-B",
                *(["-P"] if paired else []),
                "-v",
                "-o",
                output,
                bam,
            ],
            check=False,
        )
        results.append(result)
    assert results[0].returncode == 1
    assert "Unable to extrapolate" in results[0].stderr
    assert results[1].returncode == 0, results[1].stderr
    for result, distinct, one in ((results[0], 329, 0), (results[1], 658, 400)):
        assert re.search(r"TOTAL READS\s*= 1070\b", result.stderr)
        assert re.search(rf"DISTINCT READS\s*= {distinct}\b", result.stderr)
        assert re.search(rf"COUNTS OF 1\s*= {one}\b", result.stderr)
    assert "MERGED PAIRED END READS = 1070" in results[1].stderr


@pytest.mark.parametrize("layout", ["SE", "PE"])
@pytest.mark.parametrize(
    "policy,peak_mode",
    [
        ("yes", "narrow"),
        ("no", "narrow"),
        ("auto", "narrow"),
        ("auto", "broad"),
    ],
)
def test_original_metric_rules_execute_without_final_bam_or_qc_directories(
    tmp_path, tools, snakemake_executable, layout, policy, peak_mode
):
    config = project(tmp_path, layout, remove_dup=policy, peak_mode=peak_mode)
    rows = spectrum(layout)
    # The original MAPQ/primary filter, but not duplicate exclusion, must apply.
    for name, mapq, excluded_flag in (
        ("low", 0, 0),
        ("secondary", 60, 0x100),
        ("supplementary", 60, 0x800),
    ):
        extra = fragment(name, 300000, layout)
        for row in extra:
            row[4] = str(mapq)
            row[1] = str(int(row[1]) | excluded_flag)
        rows.extend(extra)
    bam = write_bam(tmp_path, tools, rows)
    sorted_bam = tmp_path / "results/S1/02_align/S1.sorted.bam"
    sorted_bam.parent.mkdir(parents=True)
    shutil.copyfile(bam, sorted_bam)
    nrf = "results/S1/01_qc/S1.nrf_pbc.tsv"
    preseq = "results/S1/05_qc/preseq/S1.preseq.txt"
    assert not (tmp_path / "results/S1/logs").exists()
    assert not (tmp_path / "results/S1/01_qc").exists()
    result = snakemake(tmp_path, snakemake_executable, config, [nrf], dry_run=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / nrf).is_file()
    assert (tmp_path / "results/S1/logs/S1.nrf_pbc.log").is_file()
    assert not (tmp_path / "results/S1/05_qc/preseq").exists()
    result = snakemake(tmp_path, snakemake_executable, config, [preseq], dry_run=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert not (tmp_path / "results/S1/02_align/S1.final.bam").exists()
    assert not (tmp_path / "results/S1/01_qc/S1.dup_metrics.txt").exists()
    with (tmp_path / nrf).open() as handle:
        assert_metrics(
            next(csv.DictReader(handle, delimiter="\t")), 1070, 658, 400, 160
        )
    lines = (tmp_path / preseq).read_text().splitlines()
    assert lines[0].startswith("TOTAL_READS\tEXPECTED_DISTINCT")
    assert len(lines) > 2
    assert (tmp_path / "results/S1/logs/S1.nrf_pbc.log").is_file()
    assert (tmp_path / "results/S1/logs/S1.preseq.log").is_file()
    assert (
        int(
            run(
                tmp_path,
                [
                    tools["samtools"],
                    "view",
                    "-c",
                    "-f",
                    "1024",
                    tmp_path / "results/S1/02_align/S1.mapq17.bam",
                ],
            ).stdout
        )
        > 0
    )


def test_uninformative_preseq_input_fails_without_fabricating_a_curve(
    tmp_path, tools, snakemake_executable
):
    config = project(tmp_path, "SE", remove_dup="no")
    bam = write_bam(
        tmp_path, tools, [alignment(f"r{i}", 100 + i * 300) for i in range(3)]
    )
    direct = run(
        tmp_path,
        [tools["preseq"], "lc_extrap", "-B", "-o", tmp_path / "direct.txt", bam],
        check=False,
    )
    assert direct.returncode == 1
    assert "max count before zero is less than min required count" in direct.stderr
    sorted_bam = tmp_path / "results/S1/02_align/S1.sorted.bam"
    sorted_bam.parent.mkdir(parents=True)
    shutil.copyfile(bam, sorted_bam)
    target = "results/S1/05_qc/preseq/S1.preseq.txt"
    result = snakemake(tmp_path, snakemake_executable, config, [target], dry_run=False)
    assert result.returncode != 0
    assert (
        "max count before zero is less than min required count"
        in (tmp_path / "results/S1/logs/S1.preseq.log").read_text()
    )
    assert not (tmp_path / target).exists()
