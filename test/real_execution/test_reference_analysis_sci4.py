"""Comprehensive real-tool and DAG contracts for SCI-4 reference analysis remediation.

Covers:
1. PE template closure (proper pair, mate reciprocity, TLEN symmetry, negative TLEN, orphan removal).
2. Resource admission validation (BED syntax, contig bounds, FAI vs chrom_sizes consistency).
3. Organelle admission and metric calculation (calc_organelle_metrics.py).
4. Dual-profile canonical analysis BAM production and downstream routing.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

from encode_pipeline.config.validator import ValidationError
from encode_pipeline.samples.admission import (
    _get_reference_contigs,
    _validate_blacklist_bed,
    validate_reference_control_admission,
)

pytestmark = pytest.mark.real_execution

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
FILTER_SCRIPT = REPO_ROOT / "scripts" / "filter_blacklist_templates.py"
ORGANELLE_SCRIPT = REPO_ROOT / "scripts" / "calc_organelle_metrics.py"

HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:10000",
    "@SQ\tSN:chr2\tLN:5000",
    "@SQ\tSN:chrM\tLN:200",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _tool(name: str, env_var: str) -> str:
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "SCI-4 real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


def _run(command: list[str], *, cwd: Path | None = None, env: dict | None = None):
    result = subprocess.run(command, capture_output=True, text=True, cwd=cwd, env=env)
    assert result.returncode == 0, (
        f"Command failed ({result.returncode}): {command!r}\n"
        f"stdout:\n{result.stdout[-2000:]}\n"
        f"stderr:\n{result.stderr[-2000:]}"
    )
    return result


def _pair(
    name: str,
    chrom: str,
    start: int,
    mate_start: int,
    tlen: int,
    flag1: int = 99,
    flag2: int = 147,
    chrom2: str | None = None,
) -> list[str]:
    c2 = chrom2 or chrom
    rnext1 = "=" if c2 == chrom else c2
    rnext2 = "=" if c2 == chrom else chrom
    return [
        "\t".join(
            [
                name,
                str(flag1),
                chrom,
                str(start),
                "60",
                "50M",
                rnext1,
                str(mate_start),
                str(tlen),
                SEQ,
                QUAL,
            ]
        ),
        "\t".join(
            [
                name,
                str(flag2),
                c2,
                str(mate_start),
                "60",
                "50M",
                rnext2,
                str(start),
                str(-tlen),
                SEQ,
                QUAL,
            ]
        ),
    ]


def _make_bam(tmp_path: Path, samtools_path: str, records: list[str]) -> Path:
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "sorted.bam"
    sam.write_text("\n".join(HEADER + records) + "\n", encoding="utf-8")
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    return bam


# ==============================================================================
# 1. Template closure & PE reciprocity tests
# ==============================================================================


def test_template_closure_drops_inconsistent_templates_as_a_pair(
    tmp_path: Path, samtools_path: str
):
    """Verify complete template closure: inconsistent or orphan records drop as a pair."""
    records = []
    # 1. Valid proper pair (retained)
    records += _pair("valid_pair", "chr1", 100, 300, 250)

    # 2. Valid proper pair where R1 is to the right (negative TLEN on R1, legal abs(TLEN))
    records += _pair("valid_neg_tlen", "chr1", 600, 400, -250, flag1=163, flag2=83)

    # 3. Orphan R1 (mate unmapped) -> dropped
    records.append(
        "\t".join(
            [
                "orphan_r1",
                "73",
                "chr1",
                "800",
                "60",
                "50M",
                "*",
                "0",
                "0",
                SEQ,
                QUAL,
            ]
        )
    )

    # 4. Incomplete template: only R2 present without R1 -> dropped
    records.append(
        "\t".join(
            [
                "incomplete_r2",
                "147",
                "chr1",
                "900",
                "60",
                "50M",
                "=",
                "750",
                "-200",
                SEQ,
                QUAL,
            ]
        )
    )

    # 5. R1/R1 pair (both reads flagged as read1, flag 99 and 99) -> dropped as a pair
    records.append(
        "\t".join(
            [
                "r1_r1_pair",
                "99",
                "chr1",
                "1200",
                "60",
                "50M",
                "=",
                "1400",
                "250",
                SEQ,
                QUAL,
            ]
        )
    )
    records.append(
        "\t".join(
            [
                "r1_r1_pair",
                "99",
                "chr1",
                "1400",
                "60",
                "50M",
                "=",
                "1200",
                "-250",
                SEQ,
                QUAL,
            ]
        )
    )

    # 6. Mates on different contigs (chr1 vs chr2) -> dropped as a pair
    records += _pair("cross_contig", "chr1", 1600, 200, 250, chrom2="chr2")

    # 7. Coordinate mismatch: R1 says mate is at 2100, but R2 is at 2200 -> dropped as a pair
    records.append(
        "\t".join(
            [
                "coord_mismatch",
                "99",
                "chr1",
                "1900",
                "60",
                "50M",
                "=",
                "2100",
                "250",
                SEQ,
                QUAL,
            ]
        )
    )
    records.append(
        "\t".join(
            [
                "coord_mismatch",
                "147",
                "chr1",
                "2200",
                "60",
                "50M",
                "=",
                "1900",
                "-250",
                SEQ,
                QUAL,
            ]
        )
    )

    # 8. TLEN mismatch: R1 says 250, R2 says -300 -> dropped as a pair
    records.append(
        "\t".join(
            [
                "tlen_mismatch",
                "99",
                "chr1",
                "2500",
                "60",
                "50M",
                "=",
                "2700",
                "250",
                SEQ,
                QUAL,
            ]
        )
    )
    records.append(
        "\t".join(
            [
                "tlen_mismatch",
                "147",
                "chr1",
                "2700",
                "60",
                "50M",
                "=",
                "2500",
                "-300",
                SEQ,
                QUAL,
            ]
        )
    )

    # 9. Zero TLEN -> dropped as a pair
    records.append(
        "\t".join(
            [
                "zero_tlen",
                "99",
                "chr1",
                "3000",
                "60",
                "50M",
                "=",
                "3200",
                "0",
                SEQ,
                QUAL,
            ]
        )
    )
    records.append(
        "\t".join(
            [
                "zero_tlen",
                "147",
                "chr1",
                "3200",
                "60",
                "50M",
                "=",
                "3000",
                "0",
                SEQ,
                QUAL,
            ]
        )
    )

    bam = _make_bam(tmp_path, samtools_path, records)
    out_bam = tmp_path / "analysis.bam"
    dropped_tsv = tmp_path / "dropped.tsv"

    _run(
        [
            "python3",
            str(FILTER_SCRIPT),
            "--bam",
            str(bam),
            "--layout",
            "PE",
            "--output",
            str(out_bam),
            "--dropped-count",
            str(dropped_tsv),
            "--sample",
            "S1",
            "--organelle-mode",
            "keep",
        ]
    )

    kept = _run([samtools_path, "view", str(out_bam)]).stdout.strip().splitlines()
    kept_names = {line.split("\t")[0] for line in kept}
    assert kept_names == {"valid_pair", "valid_neg_tlen"}, (
        f"Unexpected kept records: {kept_names}"
    )
    assert len(kept) == 4


# ==============================================================================
# 2. Resource validation (BED, contig bounds, FAI vs chrom_sizes)
# ==============================================================================


def test_blacklist_bed_validation_rules(tmp_path: Path):
    """Test blacklist BED coordinate, bounds, and syntax validation."""
    ref_contigs = {"chr1": 10000, "chr2": 5000}

    # Empty BED
    empty_bed = tmp_path / "empty.bed"
    empty_bed.write_text("# header only\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="contains no valid intervals"):
        _validate_blacklist_bed(str(empty_bed), ref_contigs=ref_contigs)

    # Negative start
    neg_bed = tmp_path / "neg.bed"
    neg_bed.write_text("chr1\t-5\t100\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="start coordinate -5 is negative"):
        _validate_blacklist_bed(str(neg_bed), ref_contigs=ref_contigs)

    # start >= end
    inv_bed = tmp_path / "inv.bed"
    inv_bed.write_text("chr1\t200\t200\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="invalid interval start=200 >= end=200"):
        _validate_blacklist_bed(str(inv_bed), ref_contigs=ref_contigs)

    # Unknown contig not in reference
    unknown_bed = tmp_path / "unknown.bed"
    unknown_bed.write_text("chrX\t100\t200\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="is absent from reference"):
        _validate_blacklist_bed(str(unknown_bed), ref_contigs=ref_contigs)

    # Out-of-bounds end coordinate
    oob_bed = tmp_path / "oob.bed"
    oob_bed.write_text("chr1\t9900\t10050\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="exceeds contig length 10000"):
        _validate_blacklist_bed(str(oob_bed), ref_contigs=ref_contigs)

    # Valid BED with comments and valid coordinates passes
    valid_bed = tmp_path / "valid.bed"
    valid_bed.write_text(
        "# track name=blacklist\nchr1\t100\t200\nchr2\t500\t1000\n",
        encoding="utf-8",
    )
    _validate_blacklist_bed(str(valid_bed), ref_contigs=ref_contigs)


def test_fai_and_chrom_sizes_agreement(tmp_path: Path):
    """Test FAI and chrom_sizes length agreement checking."""
    fai_file = tmp_path / "ref.fa.fai"
    fai_file.write_text(
        "chr1\t10000\t0\t50\t51\nchr2\t5000\t0\t50\t51\n", encoding="utf-8"
    )
    fa_file = tmp_path / "ref.fa"
    fa_file.write_text(">chr1\nACGT\n", encoding="utf-8")

    # Matching chrom.sizes
    cs_good = tmp_path / "chrom_good.sizes"
    cs_good.write_text("chr1\t10000\nchr2\t5000\n", encoding="utf-8")
    resource_good = {
        "reference_fasta": str(fa_file),
        "chrom_sizes": str(cs_good),
    }
    sizes = _get_reference_contigs(resource_good, genome="hs")
    assert sizes == {"chr1": 10000, "chr2": 5000}

    # Mismatched length on chr1
    cs_bad = tmp_path / "chrom_bad.sizes"
    cs_bad.write_text("chr1\t9999\nchr2\t5000\n", encoding="utf-8")
    resource_bad = {
        "reference_fasta": str(fa_file),
        "chrom_sizes": str(cs_bad),
    }
    with pytest.raises(ValidationError, match="length mismatch"):
        _get_reference_contigs(resource_bad, genome="hs")


def test_explicit_organelle_contig_missing_from_reference_fails(tmp_path: Path):
    """Explicitly configured organelle_chroms missing from reference must fail admission."""
    cs = tmp_path / "chrom.sizes"
    cs.write_text("chr1\t10000\nchr2\t5000\n", encoding="utf-8")
    config = {
        "use_control": False,
        "blacklist_mode": "off",
        "organelle_mode": "auto",
        "genome_resources": {
            "hs": {
                "chrom_sizes": str(cs),
                "organelle_chroms": ["chrM_typo"],
            }
        },
    }
    samples = [{"id": "S1", "genome": "hs", "role": "treatment", "assay": "atac"}]
    with pytest.raises(
        ValidationError, match="explicit organelle contig 'chrM_typo' not found"
    ):
        validate_reference_control_admission(
            samples,
            config,
            reference_mode="strict",
            external_control_policy="warn",
        )


def test_organelle_exclude_without_candidate_match_fails_closed(tmp_path: Path):
    """organelle_mode=exclude with no matching candidate contig must fail closed."""
    cs = tmp_path / "chrom.sizes"
    cs.write_text("chr1\t10000\nchr2\t5000\n", encoding="utf-8")
    config = {
        "use_control": False,
        "blacklist_mode": "off",
        "organelle_mode": "exclude",
        "genome_resources": {
            "hs": {
                "chrom_sizes": str(cs),
            }
        },
    }
    samples = [{"id": "S1", "genome": "hs", "role": "treatment", "assay": "chipseq"}]
    with pytest.raises(
        ValidationError, match="none of candidate contigs .* exist in reference"
    ):
        validate_reference_control_admission(
            samples,
            config,
            reference_mode="strict",
            external_control_policy="warn",
        )


# ==============================================================================
# 3. Organelle metrics script contracts
# ==============================================================================


def test_calc_organelle_metrics_statuses(tmp_path: Path, samtools_path: str):
    """Test calc_organelle_metrics.py reporting, NA handling, and unresolved status."""
    records = []
    records += _pair("n1", "chr1", 100, 300, 250)
    records += _pair("m1", "chrM", 10, 80, 80)
    bam = _make_bam(tmp_path, samtools_path, records)

    env = os.environ.copy()
    env["PATH"] = f"{Path(samtools_path).parent}:{env.get('PATH', '')}"

    # 1. Normal matching: chrM resolved, 2/4 reads = 0.5 fraction
    tsv1 = tmp_path / "mito1.tsv"
    _run(
        [
            "python3",
            str(ORGANELLE_SCRIPT),
            "--bam",
            str(bam),
            "--sample",
            "S1",
            "--mode",
            "report",
            "--chroms",
            "chrM",
            "--output",
            str(tsv1),
        ],
        env=env,
    )
    rows1 = tsv1.read_text(encoding="utf-8").strip().splitlines()
    assert rows1[1].split("\t") == ["S1", "2", "4", "0.500000", "report", "ok"]

    # 2. Unresolved candidate: chrMT not in @SQ -> unresolved_contigs and NA
    tsv2 = tmp_path / "mito2.tsv"
    _run(
        [
            "python3",
            str(ORGANELLE_SCRIPT),
            "--bam",
            str(bam),
            "--sample",
            "S1",
            "--mode",
            "report",
            "--chroms",
            "chrMT_nonexistent",
            "--output",
            str(tsv2),
        ],
        env=env,
    )
    rows2 = tsv2.read_text(encoding="utf-8").strip().splitlines()
    assert rows2[1].split("\t") == [
        "S1",
        "NA",
        "4",
        "NA",
        "report",
        "unresolved_contigs",
    ]

    # 3. Mode exclude fails closed when candidate not in @SQ
    tsv3 = tmp_path / "mito3.tsv"
    res = subprocess.run(
        [
            "python3",
            str(ORGANELLE_SCRIPT),
            "--bam",
            str(bam),
            "--sample",
            "S1",
            "--mode",
            "exclude",
            "--chroms",
            "chrMT_nonexistent",
            "--output",
            str(tsv3),
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert res.returncode != 0
    assert "exclude" in res.stderr and "none of candidate contigs" in res.stderr


# ==============================================================================
# 4. Snakemake execution & dual-profile analysis.bam downstream wiring
# ==============================================================================


@pytest.mark.parametrize("profile", ["legacy", "updated"])
def test_snakemake_dual_profile_analysis_bam_and_downstream_routing(
    tmp_path: Path, snakemake_path: str, samtools_path: str, profile: str
):
    """Verify both legacy and updated profiles produce canonical analysis.bam and route downstream."""
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t10\t50\n", encoding="utf-8")
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr1\t10000\nchr2\t5000\nchrM\t200\n", encoding="utf-8")

    (tmp_path / "samples.tsv").write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "S1\tR1.fq\tR2.fq\tPE\tcuttag\tH3K4me3\tnarrow\ths\tbt2_idx\t"
        "exp1\tH3K4me3\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )

    outdir_name = "results"
    cfg_lines = [
        "samples: samples.tsv",
        f"outdir: {outdir_name}",
        "threads: 1",
        "use_control: false",
        "multiqc: false",
        "replicate_analysis: false",
        "genome_resources:",
        "  hs:",
        '    effective_genome_size: "hs"',
        f'    chrom_sizes: "{chrom_sizes}"',
        f'    blacklist: "{blacklist}"',
        '    organelle_chroms: ["chrM"]',
        "config_schema_version: 2.0.0",
        f"analysis_profile: {profile}",
        "blacklist_mode: analysis",
        "organelle_mode: auto",
        "duplicate_analysis: exclude",
        "reference_mode: legacy",
        "external_control_policy: warn",
        'fragment_cpm: "off"',
        "tss_fragment_profile: false",
        "qc:",
        "  summary: false",
        "  frip: false",
        "  signal_tracks: false",
    ]
    if profile == "updated":
        cfg_lines.append("duplicate_backend: picard")
        cfg_lines.append("control_mode: matched")

    (tmp_path / "config.yaml").write_text("\n".join(cfg_lines) + "\n", encoding="utf-8")

    # 1. Initialize analysis identity
    _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            "config.yaml",
            "--workflow-profile",
            "none",
            "-j",
            "1",
            "results/.helixweave/analysis_identity.json",
        ],
        cwd=tmp_path,
    )

    # 2. Write mock final.bam with normal read, duplicate, blacklist, and mito reads
    mock_final_bam = tmp_path / "results" / "S1" / "02_align" / "S1.final.bam"
    mock_final_bam.parent.mkdir(parents=True, exist_ok=True)
    records = []
    records += _pair("normal", "chr1", 500, 600, 150)
    records += _pair("dup", "chr1", 700, 800, 150, flag1=99 + 1024, flag2=147 + 1024)
    records += _pair("bl_hit", "chr1", 20, 100, 120)
    records += _pair("mito_read", "chrM", 10, 80, 80)
    sam = tmp_path / "mock.sam"
    raw = tmp_path / "mock.raw.bam"
    sam.write_text("\n".join(HEADER + records) + "\n", encoding="utf-8")
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(mock_final_bam), str(raw)])
    _run([samtools_path, "index", str(mock_final_bam)])

    # Target analysis.bam
    target_analysis_bam = "results/analysis_filtered/S1/02_align/S1.analysis.bam"
    _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            "config.yaml",
            "--workflow-profile",
            "none",
            "-j",
            "1",
            target_analysis_bam,
        ],
        cwd=tmp_path,
    )

    analysis_bam_path = tmp_path / target_analysis_bam
    assert analysis_bam_path.exists(), f"analysis.bam missing: {target_analysis_bam}"
    assert (tmp_path / f"{target_analysis_bam}.bai").exists()

    kept_records = (
        _run([samtools_path, "view", str(analysis_bam_path)])
        .stdout.strip()
        .splitlines()
    )
    assert len(kept_records) == 2
    assert {r.split("\t")[0] for r in kept_records} == {"normal"}
