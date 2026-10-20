"""Real-tool contract for fragment-level CPM coverage semantics."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"

HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:10000",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "fragment CPM real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


@pytest.fixture(scope="module")
def bamcoverage_path():
    return _tool("bamCoverage", "BAMCOVERAGE")


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


def _run(command, *, env=None, cwd=None):
    result = subprocess.run(command, capture_output=True, text=True, env=env, cwd=cwd)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-1000:]}\n"
        f"stderr:\n{result.stderr[-1000:]}"
    )
    return result


def _age_inputs(bam):
    """Avoid Snakemake clock-skew false failures in fast fixtures."""
    for path in (bam, Path(str(bam) + ".bai")):
        os.utime(path, (1, 1))


def _pair(name, start, mate_start, tlen):
    return [
        "\t".join(
            [
                name,
                "99",
                "chr1",
                str(start),
                "60",
                "50M",
                "=",
                str(mate_start),
                str(tlen),
                SEQ,
                QUAL,
            ]
        ),
        "\t".join(
            [
                name,
                "147",
                "chr1",
                str(mate_start),
                "60",
                "50M",
                "=",
                str(start),
                str(-tlen),
                SEQ,
                QUAL,
            ]
        ),
    ]


def test_fragment_cpm_counts_one_fragment_per_pair(
    tmp_path, samtools_path, bamcoverage_path
):
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "input.bam"
    bedgraph = tmp_path / "fragment.bedgraph"
    sam.write_text(
        "\n".join(
            HEADER + _pair("frag1", 100, 300, 250) + _pair("frag2", 500, 700, 250)
        )
        + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(
        [str(Path(bamcoverage_path).parent), env.get("PATH", "")]
    )
    _run(
        [
            bamcoverage_path,
            "--bam",
            str(bam),
            "--outFileName",
            str(bedgraph),
            "--outFileFormat",
            "bedgraph",
            "--extendReads",
            "--samFlagInclude",
            "64",
            "--binSize",
            "1",
        ],
        env=env,
    )
    covered = 0
    signal_sum = 0.0
    max_height = 0.0
    for line in bedgraph.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "track")):
            continue
        fields = line.split("\t")
        value = float(fields[3])
        if value <= 0:
            continue
        width = int(fields[2]) - int(fields[1])
        covered += width
        signal_sum += width * value
        max_height = max(max_height, value)
    assert covered == 500
    assert signal_sum == 500
    assert max_height == 1


def test_se_fragment_production_rule_uses_explicit_extension(
    tmp_path, samtools_path, snakemake_path
):
    (tmp_path / "samples.tsv").write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "T1\tR1.fq\t\tSE\tcuttag\tH3K4me3\tnarrow\ths\tbt2_idx\t"
        "exp1\tH3K4me3\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    (tmp_path / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "extend_reads: 150",
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                "config_schema_version: 2.0.0",
                "analysis_profile: updated",
                "fragment_cpm: auto",
                "organelle_mode: auto",
                "duplicate_backend: picard",
                "duplicate_analysis: exclude",
                "control_mode: matched",
                "reference_mode: legacy",
                "external_control_policy: warn",
                "tss_fragment_profile: false",
                "qc:",
                "  summary: false",
                "  frip: false",
                "  signal_tracks: false",
                "",
            ]
        ),
        encoding="utf-8",
    )
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
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "results/T1/02_align/T1.analysis.bam"
    bam.parent.mkdir(parents=True)
    sam.write_text(
        "\n".join(
            HEADER
            + [
                "\t".join(
                    ["read1", "0", "chr1", "100", "60", "50M", "*", "0", "0", SEQ, QUAL]
                )
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    _age_inputs(bam)
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
            "results/T1/03_bigwig/T1.fragment.CPM.bw",
        ],
        cwd=tmp_path,
    )
    output = tmp_path / "results/T1/03_bigwig/T1.fragment.CPM.bw"
    assert output.exists() and output.stat().st_size > 0
    log = (tmp_path / "results/T1/logs/T1.fragment_cpm.log").read_text(encoding="utf-8")
    assert "not paired-end" not in log


def _bigwig_positive_width(path):
    python = shutil.which("python3")
    if not python:
        pytest.skip("python3 unavailable for pyBigWig inspection")
    code = (
        "import pyBigWig,sys; bw=pyBigWig.open(sys.argv[1]); "
        "print(sum(e-s for s,e,v in (bw.intervals('chr1') or []) if v > 0))"
    )
    result = subprocess.run(
        [python, "-c", code, str(path)], capture_output=True, text=True
    )
    if result.returncode != 0:
        pytest.skip(f"pyBigWig unavailable: {result.stderr[-200:]}")
    return int(result.stdout.strip())


def _prepare_se_broad_workflow(tmp_path, samtools_path, snakemake_path, model_length):
    (tmp_path / "samples.tsv").write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "T1\tR1.fq\t\tSE\tcuttag\tH3K27me3\tbroad\ths\tbt2_idx\t"
        "exp1\tH3K27me3\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    (tmp_path / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "extend_reads: auto",
                "binsize: 1",
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                "config_schema_version: 2.0.0",
                "analysis_profile: updated",
                "fragment_cpm: auto",
                "organelle_mode: auto",
                "duplicate_backend: picard",
                "duplicate_analysis: exclude",
                "control_mode: matched",
                "reference_mode: legacy",
                "external_control_policy: warn",
                "tss_fragment_profile: false",
                "qc:",
                "  summary: false",
                "  frip: false",
                "  signal_tracks: false",
                "",
            ]
        ),
        encoding="utf-8",
    )
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
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "results/T1/02_align/T1.analysis.bam"
    bam.parent.mkdir(parents=True)
    sam.write_text(
        "\n".join(
            HEADER
            + [
                "\t".join(
                    ["read1", "0", "chr1", "100", "60", "50M", "*", "0", "0", SEQ, QUAL]
                )
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    _age_inputs(bam)
    peaks = tmp_path / "results/T1/04_peaks/T1"
    peaks.mkdir(parents=True)
    (peaks / "T1_peaks.broadPeak").write_text(
        "chr1\t100\t200\tpeak\t100\t.\t1\t2\t3\n", encoding="utf-8"
    )
    log = tmp_path / "results/T1/logs/T1.macs3.log"
    log.parent.mkdir(parents=True)
    if model_length is not None:
        log.write_text(
            f"predicted fragment length is {model_length}\n", encoding="utf-8"
        )
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
            "results/T1/03_bigwig/T1.fragment.CPM.bw",
        ],
        cwd=tmp_path,
    )
    return tmp_path / "results/T1/03_bigwig/T1.fragment.CPM.bw"


@pytest.mark.parametrize("model_length", [None, 300])
def test_se_broad_fragment_uses_macs3_model_when_available(
    tmp_path, samtools_path, snakemake_path, model_length
):
    output = _prepare_se_broad_workflow(
        tmp_path, samtools_path, snakemake_path, model_length
    )
    width = _bigwig_positive_width(output)
    expected = 200 if model_length is None else 300
    assert width == expected


def test_se_broad_fragment_cold_run_schedules_model_producer(
    tmp_path, samtools_path, snakemake_path
):
    (tmp_path / "samples.tsv").write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "T1\tR1.fq\t\tSE\tcuttag\tH3K27me3\tbroad\ths\tbt2_idx\t"
        "exp1\tH3K27me3\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    (tmp_path / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "extend_reads: auto",
                "binsize: 1",
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                "config_schema_version: 2.0.0",
                "analysis_profile: updated",
                "fragment_cpm: auto",
                "organelle_mode: auto",
                "duplicate_backend: picard",
                "duplicate_analysis: exclude",
                "control_mode: matched",
                "reference_mode: legacy",
                "external_control_policy: warn",
                "tss_fragment_profile: false",
                "qc:",
                "  summary: false",
                "  frip: false",
                "  signal_tracks: false",
                "",
            ]
        ),
        encoding="utf-8",
    )
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
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "results/T1/02_align/T1.analysis.bam"
    bam.parent.mkdir(parents=True)
    sam.write_text(
        "\n".join(
            HEADER
            + [
                "\t".join(
                    ["read1", "0", "chr1", "100", "60", "50M", "*", "0", "0", SEQ, QUAL]
                )
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    _age_inputs(bam)

    bindir = tmp_path / "bin"
    bindir.mkdir()
    macs3 = bindir / "macs3"
    macs3.write_text(
        """#!/bin/bash
set -e
outdir=""
name=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) outdir="$2"; shift 2 ;;
    -n) name="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$outdir"
printf 'chr1\\t100\\t200\\tpeak\\t100\\t.\\t1\\t2\\t3\\n' > "$outdir/${name}_peaks.broadPeak"
echo "predicted fragment length is 300"
""",
        encoding="utf-8",
    )
    macs3.chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(bindir), env.get("PATH", "")])
    target = "results/T1/03_bigwig/T1.fragment.CPM.bw"
    first = _run(
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
            target,
        ],
        cwd=tmp_path,
        env=env,
    )
    assert "macs3_callpeak" in first.stdout + first.stderr
    output = tmp_path / target
    assert _bigwig_positive_width(output) == 300

    second = _run(
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
            target,
        ],
        cwd=tmp_path,
        env=env,
    )
    assert "Nothing to be done" in second.stdout + second.stderr
