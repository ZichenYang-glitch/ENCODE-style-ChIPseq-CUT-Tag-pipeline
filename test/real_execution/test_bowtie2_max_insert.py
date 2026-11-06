"""Real-tool verification for SCI-3 Bowtie2 paired-end max_insert policy.

Verifies against the pinned bowtie2, bowtie2-build, samtools, macs3, and the
production Snakemake workflow rules that:

- Bowtie2 -X controls the concordant proper-pair (0x2) and TLEN boundary at
  500/501/600, 700/701, and 2000/2001 bp in both F1R2 and F2R1 orientations;
- The production bowtie2_align -> samtools_filter -> macs3_callpeak chain
  admits >500 bp templates (e.g. 600 bp in CUT&Tag auto=700 and 1200 bp in
  ChIP-seq/ATAC auto=2000) as proper pairs into BAMPE pileups, whereas
  max_insert=500 leaves them without 0x2 (dropped by updated samtools_filter
  and ignored by MACS3 BAMPE in both profiles);
- Mixed-assay auto, SE non-rendering, control samples, YAML/CLI overrides,
  extra_args -X/--maxins migration (preserving quotes, variables, and
  backslash-escaped spaces), and startup conflict/bounds rejections work on
  the production Snakemake entrypoint.
"""

import gzip
import json
import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"

REF_LEN = 30_000
READ_LEN = 50
SEED = 1337
QUAL = "I" * READ_LEN
COMPLEMENT = str.maketrans("ATGC", "TACG")

SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
    "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
)


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "Bowtie2 max_insert real execution")
    return resolved


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


@pytest.fixture(scope="module")
def bowtie2_path():
    return _tool("bowtie2", "BOWTIE2")


@pytest.fixture(scope="module")
def bowtie2_build_path():
    return _tool("bowtie2-build", "BOWTIE2_BUILD")


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


@pytest.fixture(scope="module")
def macs3_path():
    return _tool("macs3", "MACS3")


def _revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _run(command, *, env=None, cwd=None, check=True):
    result = subprocess.run(command, capture_output=True, text=True, env=env, cwd=cwd)
    if check:
        assert result.returncode == 0, (
            f"Command failed: {command!r}\n"
            f"stdout:\n{result.stdout[-1000:]}\n"
            f"stderr:\n{result.stderr[-1000:]}"
        )
    return result


def _snakemake_env(*tool_paths):
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    dirs = []
    for tool in tool_paths:
        parent = str(Path(tool).parent)
        if parent not in dirs:
            dirs.append(parent)
    conda = shutil.which("conda")
    if conda:
        conda_dir = str(Path(conda).parent)
        if conda_dir not in dirs:
            dirs.append(conda_dir)
    env["PATH"] = ":".join([*dirs, "/usr/bin", "/bin"])
    return env


@pytest.fixture(scope="module")
def reference_bundle(tmp_path_factory, bowtie2_build_path):
    """Build a reusable 30 kb pseudo-random reference and Bowtie2 index."""
    ref_dir = tmp_path_factory.mktemp("bt2_ref")
    rng = random.Random(SEED)
    seq = "".join(rng.choice("ATGC") for _ in range(REF_LEN))
    ref_fa = ref_dir / "ref.fa"
    with open(ref_fa, "w", encoding="utf-8") as handle:
        handle.write(">chr1\n")
        for i in range(0, len(seq), 70):
            handle.write(seq[i : i + 70] + "\n")
    idx_prefix = ref_dir / "bt2_idx"
    _run([bowtie2_build_path, str(ref_fa), str(idx_prefix)], cwd=ref_dir)
    return {"seq": seq, "ref_fa": ref_fa, "idx_prefix": idx_prefix}


def _emit_pair(ref_seq: str, start: int, length: int, orientation: str):
    """Return (r1_seq, r2_seq) for an insert of *length* at 0-based *start*.

    orientation='F1R2': R1 is forward at [start, start+READ_LEN),
                        R2 is reverse-complement at [start+length-READ_LEN, start+length).
    orientation='F2R1': R1 is reverse-complement at [start+length-READ_LEN, start+length),
                        R2 is forward at [start, start+READ_LEN).
    """
    left_fwd = ref_seq[start : start + READ_LEN]
    right_rc = _revcomp(ref_seq[start + length - READ_LEN : start + length])
    if orientation == "F1R2":
        return left_fwd, right_rc
    if orientation == "F2R1":
        return right_rc, left_fwd
    raise ValueError(f"Unsupported orientation: {orientation}")


def _parse_sam_records(sam_text: str):
    records = {}
    for line in sam_text.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        qname = fields[0]
        flag = int(fields[1])
        pos = int(fields[3])
        tlen = int(fields[8])
        yt = next((f[5:] for f in fields[11:] if f.startswith("YT:Z:")), None)
        mate = 1 if (flag & 0x40) else 2
        records[(qname, mate)] = {"flag": flag, "pos": pos, "tlen": tlen, "yt": yt}
    return records


@pytest.mark.parametrize(
    ("max_insert", "concordant_lengths", "discordant_lengths"),
    [
        (500, [500], [501, 600, 700, 701, 2000, 2001]),
        (700, [500, 501, 600, 700], [701, 2000, 2001]),
        (2000, [500, 501, 600, 700, 701, 2000], [2001]),
    ],
)
def test_real_bowtie2_max_insert_boundaries_both_orientations(
    tmp_path,
    bowtie2_path,
    reference_bundle,
    max_insert,
    concordant_lengths,
    discordant_lengths,
):
    """Verify SAM 0x2 proper-pair flags, YT:Z tags, and signed TLEN across
    500/501/600/700/701/2000/2001 in both F1R2 and F2R1 orientations against
    the real bowtie2 binary."""
    ref_seq = reference_bundle["seq"]
    idx_prefix = reference_bundle["idx_prefix"]
    lengths = [500, 501, 600, 700, 701, 2000, 2001]

    r1_path = tmp_path / "R1.fq"
    r2_path = tmp_path / "R2.fq"
    with (
        open(r1_path, "w", encoding="utf-8") as f1,
        open(r2_path, "w", encoding="utf-8") as f2,
    ):
        for idx, length in enumerate(lengths):
            for orient_idx, orientation in enumerate(("F1R2", "F2R1")):
                start = 1_000 + (idx * 2 + orient_idx) * 1_800
                r1_seq, r2_seq = _emit_pair(ref_seq, start, length, orientation)
                qname = f"len{length}_{orientation}"
                f1.write(f"@{qname}/1\n{r1_seq}\n+\n{QUAL}\n")
                f2.write(f"@{qname}/2\n{r2_seq}\n+\n{QUAL}\n")

    sam_out = _run(
        [
            bowtie2_path,
            "-x",
            str(idx_prefix),
            "-1",
            str(r1_path),
            "-2",
            str(r2_path),
            "-X",
            str(max_insert),
            "--no-unal",
        ]
    ).stdout
    records = _parse_sam_records(sam_out)

    for length in concordant_lengths:
        for orientation in ("F1R2", "F2R1"):
            qname = f"len{length}_{orientation}"
            m1 = records[(qname, 1)]
            m2 = records[(qname, 2)]
            assert (m1["flag"] & 0x2) == 0x2, (
                f"{qname} mate1 not proper pair at -X {max_insert}"
            )
            assert (m2["flag"] & 0x2) == 0x2, (
                f"{qname} mate2 not proper pair at -X {max_insert}"
            )
            assert m1["yt"] == "CP"
            assert m2["yt"] == "CP"
            if orientation == "F1R2":
                assert m1["tlen"] == length
                assert m2["tlen"] == -length
            else:
                assert m1["tlen"] == -length
                assert m2["tlen"] == length

    for length in discordant_lengths:
        for orientation in ("F1R2", "F2R1"):
            qname = f"len{length}_{orientation}"
            m1 = records[(qname, 1)]
            m2 = records[(qname, 2)]
            assert (m1["flag"] & 0x2) == 0, (
                f"{qname} mate1 unexpectedly proper pair at -X {max_insert}"
            )
            assert (m2["flag"] & 0x2) == 0, (
                f"{qname} mate2 unexpectedly proper pair at -X {max_insert}"
            )
            assert m1["yt"] == "DP"
            assert m2["yt"] == "DP"
            if orientation == "F1R2":
                assert m1["tlen"] == length
                assert m2["tlen"] == -length
            else:
                assert m1["tlen"] == -length
                assert m2["tlen"] == length


def _pileup_at(bdg_path: Path, position: int) -> float:
    assert bdg_path.exists(), f"missing {bdg_path}"
    for line in bdg_path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 4 and int(fields[1]) <= position < int(fields[2]):
            return float(fields[3])
    return 0.0


@pytest.mark.parametrize("profile", ["legacy", "updated"])
def test_production_workflow_admits_long_pe_fragments_into_macs3_bampe(
    tmp_path,
    snakemake_path,
    bowtie2_path,
    samtools_path,
    macs3_path,
    reference_bundle,
    profile,
):
    """Prove through production Snakemake rules (trim_galore -> bowtie2_align
    -> samtools_filter -> markdup/prepare_analysis_bam -> macs3_callpeak) that:
    - a cluster of 600 bp fragments at pos 5000 and 1200 bp fragments at pos 15000
      both align as proper pairs (0x2) and produce MACS3 BAMPE pileup under
      chipseq auto (max_insert=2000);
    - under cuttag auto (max_insert=700), the 600 bp cluster is proper-paired
      while the 1200 bp cluster is not (and in updated profile is removed from
      mapq30.bam/analysis.bam by samtools_filter -f 2);
    - under explicit max_insert=500 (the old implicit Bowtie2 default), neither
      the 600 bp nor the 1200 bp cluster is proper-paired, and MACS3 BAMPE
      pileup at both loci drops to 0.0.
    """
    ref_seq = reference_bundle["seq"]
    idx_prefix = reference_bundle["idx_prefix"]
    env = _snakemake_env(bowtie2_path, samtools_path, macs3_path)

    pairs = []
    # 30 short 200 bp background pairs around pos 1000 so MACS3 always has input
    for i in range(30):
        start = 500 + i * 5
        r1, r2 = _emit_pair(ref_seq, start, 200, "F1R2")
        pairs.append((f"short_{i}", r1, r2))
    # 30 pairs of 600 bp around pos 5000 (center ~ 5350)
    for i in range(30):
        start = 5_000 + i * 5
        r1, r2 = _emit_pair(ref_seq, start, 600, "F1R2" if i % 2 == 0 else "F2R1")
        pairs.append((f"mid600_{i}", r1, r2))
    # 30 pairs of 1200 bp around pos 15000 (center ~ 15650)
    for i in range(30):
        start = 15_000 + i * 5
        r1, r2 = _emit_pair(ref_seq, start, 1200, "F1R2" if i % 2 == 0 else "F2R1")
        pairs.append((f"long1200_{i}", r1, r2))

    def run_scenario(name: str, assay: str, max_insert_cfg, run_macs3: bool):
        workdir = tmp_path / f"{profile}_{name}"
        workdir.mkdir(parents=True)
        raw1 = workdir / "T1_R1.fq.gz"
        raw2 = workdir / "T1_R2.fq.gz"
        with (
            gzip.open(raw1, "wt", encoding="utf-8") as f1,
            gzip.open(raw2, "wt", encoding="utf-8") as f2,
        ):
            for qname, r1_seq, r2_seq in pairs:
                f1.write(f"@{qname}/1\n{r1_seq}\n+\n{QUAL}\n")
                f2.write(f"@{qname}/2\n{r2_seq}\n+\n{QUAL}\n")

        (workdir / "samples.tsv").write_text(
            SAMPLES_HEADER
            + f"T1\tT1_R1.fq.gz\tT1_R2.fq.gz\tPE\t{assay}\tH3K27ac\tnarrow\ths\t"
            f"{idx_prefix}\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n",
            encoding="utf-8",
        )
        config = {
            "samples": "samples.tsv",
            "outdir": "results",
            "threads": 1,
            "mapq": 30,
            "trim": False,
            "use_control": False,
            "multiqc": False,
            "replicate_analysis": False,
            "chipseq_idr": False,
            "duplicate_analysis": "retain",
            "tool_parameters": {
                "bowtie2": {"max_insert": max_insert_cfg},
                "macs3": {"extra_args": "-B --nomodel"},
            },
            "genome_resources": {
                "hs": {
                    "effective_genome_size": str(REF_LEN),
                    "chrom_sizes": "",
                    "blacklist": "",
                    "gtf": "",
                    "reference_fasta": "",
                }
            },
            "qc": {
                "blacklist_filter": False,
                "frip": False,
                "library_complexity": False,
                "nrf_pbc": False,
                "signal_tracks": False,
                "summary": False,
                "cuttag_fragment_size": False,
            },
        }
        if profile == "updated":
            config.update(
                {
                    "config_schema_version": "2.0.0",
                    "analysis_profile": "updated",
                    "blacklist_mode": "off",
                    "organelle_mode": "keep",
                    "duplicate_backend": "picard",
                    "control_mode": "matched",
                    "reference_mode": "legacy",
                    "external_control_policy": "warn",
                    "fragment_cpm": "off",
                    "tss_fragment_profile": False,
                }
            )
        (workdir / "config.yaml").write_text(
            yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
        )
        if profile == "legacy":
            bam_target = "results/T1/02_align/T1.final.bam"
            targets = [bam_target]
            if run_macs3:
                targets.append("results/T1/04_peaks/T1")
            _run(
                [
                    snakemake_path,
                    "-s",
                    str(SNAKEFILE),
                    "--configfile",
                    str(workdir / "config.yaml"),
                    "-j",
                    "1",
                    *targets,
                ],
                env=env,
                cwd=workdir,
            )
        else:
            # In updated profile, run trim_galore -> bowtie2_align -> samtools_filter (-f 2)
            mapq_bam = "results/T1/02_align/T1.mapq30.bam"
            _run(
                [
                    snakemake_path,
                    "-s",
                    str(SNAKEFILE),
                    "--configfile",
                    str(workdir / "config.yaml"),
                    "-j",
                    "1",
                    mapq_bam,
                ],
                env=env,
                cwd=workdir,
            )
            bam_target = mapq_bam
            if run_macs3:
                # Stage T1.final.bam from T1.mapq30.bam so prepare_analysis_bam -> macs3_callpeak
                # runs without requiring a Picard Java runtime in the minimal test environment.
                final_bam = workdir / "results/T1/02_align/T1.final.bam"
                shutil.copy(workdir / mapq_bam, final_bam)
                _run([samtools_path, "index", str(final_bam)], env=env)
                bam_target = "results/T1/02_align/T1.analysis.bam"
                _run(
                    [
                        snakemake_path,
                        "-s",
                        str(SNAKEFILE),
                        "--configfile",
                        str(workdir / "config.yaml"),
                        "-j",
                        "1",
                        bam_target,
                        "results/T1/04_peaks/T1",
                    ],
                    env=env,
                    cwd=workdir,
                )

        proper_view = _run(
            [samtools_path, "view", "-f", "2", str(workdir / bam_target)], env=env
        ).stdout
        proper_qnames = {
            line.split("\t")[0] for line in proper_view.splitlines() if line
        }
        all_view = _run(
            [samtools_path, "view", str(workdir / bam_target)], env=env
        ).stdout
        all_qnames = {line.split("\t")[0] for line in all_view.splitlines() if line}
        bdg = workdir / "results" / "T1" / "04_peaks" / "T1" / "T1_treat_pileup.bdg"
        return proper_qnames, all_qnames, bdg

    # 1) chipseq auto (2000): admits short (200), mid (600), and long (1200) into proper pairs and MACS3 BAMPE
    proper_2000, _, bdg_2000 = run_scenario("chipseq_auto", "chipseq", "auto", True)
    assert "mid600_0" in proper_2000
    assert "long1200_0" in proper_2000
    assert _pileup_at(bdg_2000, 5_350) >= 25.0
    assert _pileup_at(bdg_2000, 15_650) >= 25.0

    # 2) cuttag auto (700): admits short (200) and mid (600) as proper pairs, excludes long (1200)
    proper_700, all_700, _ = run_scenario("cuttag_auto", "cuttag", "auto", False)
    assert "mid600_0" in proper_700
    assert "long1200_0" not in proper_700
    if profile == "updated":
        # updated samtools_filter enforces -f 2 so non-proper pairs are absent from mapq30.bam
        assert "long1200_0" not in all_700

    # 3) explicit max_insert=500: excludes both mid (600) and long (1200) from proper pairs and MACS3 BAMPE pileup
    proper_500, all_500, bdg_500 = run_scenario("explicit_500", "chipseq", 500, True)
    assert "short_0" in proper_500
    assert "mid600_0" not in proper_500
    assert "long1200_0" not in proper_500
    if profile == "updated":
        assert "mid600_0" not in all_500
        assert "long1200_0" not in all_500
    assert _pileup_at(bdg_500, 5_350) == 0.0
    assert _pileup_at(bdg_500, 15_650) == 0.0


def test_production_bowtie2_rule_mixed_assays_se_control_and_cli_override(
    tmp_path,
    snakemake_path,
    bowtie2_path,
    samtools_path,
    reference_bundle,
):
    """Run the production bowtie2_align rule on a mixed-assay project containing
    all 4 assays (chipseq, atac, cuttag, mnase) plus an SE control sample, then
    verify per-sample proper-pair outcomes and CLI override via --config."""
    ref_seq = reference_bundle["seq"]
    idx_prefix = reference_bundle["idx_prefix"]
    env = _snakemake_env(bowtie2_path, samtools_path)

    # Probe pairs at 500, 600, and 1200 bp
    probe_pairs = [
        ("p500", *_emit_pair(ref_seq, 2_000, 500, "F1R2")),
        ("p600", *_emit_pair(ref_seq, 6_000, 600, "F2R1")),
        ("p1200", *_emit_pair(ref_seq, 12_000, 1200, "F1R2")),
    ]
    workdir = tmp_path / "mixed_real"
    workdir.mkdir()
    r1 = workdir / "R1.fq.gz"
    r2 = workdir / "R2.fq.gz"
    with (
        gzip.open(r1, "wt", encoding="utf-8") as f1,
        gzip.open(r2, "wt", encoding="utf-8") as f2,
    ):
        for qname, s1, s2 in probe_pairs:
            f1.write(f"@{qname}/1\n{s1}\n+\n{QUAL}\n")
            f2.write(f"@{qname}/2\n{s2}\n+\n{QUAL}\n")

    rows = [
        f"S_chip\tR1.fq.gz\tR2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t{idx_prefix}\texp1\tH3K27ac\t1\t1\t1\ttreatment\tC_se\t",
        f"S_atac\tR1.fq.gz\tR2.fq.gz\tPE\tatac\tATAC\tnarrow\ths\t{idx_prefix}\texp2\tATAC\t1\t1\t1\ttreatment\t\t",
        f"S_cut\tR1.fq.gz\tR2.fq.gz\tPE\tcuttag\tH3K4me3\tnarrow\ths\t{idx_prefix}\texp3\tH3K4me3\t1\t1\t1\ttreatment\t\t",
        f"S_mnase\tR1.fq.gz\tR2.fq.gz\tPE\tmnase\tH3\tnucleosome\ths\t{idx_prefix}\texp4\tH3\t1\t1\t1\ttreatment\t\t",
        f"C_se\tR1.fq.gz\t\tSE\tchipseq\tInput\tnarrow\ths\t{idx_prefix}\texp1\tInput\t1\t1\t1\tcontrol\t\t",
    ]
    (workdir / "samples.tsv").write_text(
        SAMPLES_HEADER + "\n".join(rows) + "\n", encoding="utf-8"
    )
    config = {
        "samples": "samples.tsv",
        "outdir": "results",
        "threads": 1,
        "mapq": 30,
        "trim": False,
        "use_control": True,
        "multiqc": False,
        "replicate_analysis": False,
        "chipseq_idr": False,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": str(REF_LEN)}},
    }
    (workdir / "config.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    targets = [
        f"results/{s}/02_align/{s}.sorted.bam"
        for s in ("S_chip", "S_atac", "S_cut", "S_mnase", "C_se")
    ]
    _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(workdir / "config.yaml"),
            "-j",
            "1",
            *targets,
        ],
        env=env,
        cwd=workdir,
    )

    def proper_set(sample: str):
        out = _run(
            [
                samtools_path,
                "view",
                "-f",
                "2",
                str(workdir / f"results/{sample}/02_align/{sample}.sorted.bam"),
            ],
            env=env,
        ).stdout
        return {line.split("\t")[0] for line in out.splitlines() if line}

    assert proper_set("S_chip") == {"p500", "p600", "p1200"}
    assert proper_set("S_atac") == {"p500", "p600", "p1200"}
    assert proper_set("S_cut") == {"p500", "p600"}
    assert proper_set("S_mnase") == {"p500"}
    assert proper_set("C_se") == set()

    # Now verify CLI --config override (max_insert=650) in a fresh outdir
    _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(workdir / "config.yaml"),
            "--config",
            "outdir=results_cli",
            "tool_parameters={'bowtie2': {'max_insert': 650}}",
            "-j",
            "1",
            "results_cli/S_mnase/02_align/S_mnase.sorted.bam",
            "results_cli/S_chip/02_align/S_chip.sorted.bam",
        ],
        env=env,
        cwd=workdir,
    )
    for sample in ("S_mnase", "S_chip"):
        out = _run(
            [
                samtools_path,
                "view",
                "-f",
                "2",
                str(workdir / f"results_cli/{sample}/02_align/{sample}.sorted.bam"),
            ],
            env=env,
        ).stdout
        assert {line.split("\t")[0] for line in out.splitlines() if line} == {
            "p500",
            "p600",
        }


def test_production_bowtie2_rule_migrates_extra_args_and_preserves_quotes_and_vars(
    tmp_path,
    snakemake_path,
    bowtie2_path,
    samtools_path,
    reference_bundle,
):
    """Verify on the real bowtie2_align rule that:
    - a line-continued --max\\nins 800 in extra_args migrates to max_insert=800
      so a 600 bp pair aligns concordantly in an mnase sample (whose auto default is 500);
    - quoted spaces, backslash-escaped trailing spaces, and runtime shell variable
      expansion in extra_args (--rg-id "$RG_ID", --rg 'DS:Desc With Spaces', and
      --rg CN:Center\\ ) survive normalization into the BAM header;
    - conflicting --maxins or line-continued -X\\n 900 vs max_insert=700 fails closed at startup.
    """
    ref_seq = reference_bundle["seq"]
    idx_prefix = reference_bundle["idx_prefix"]
    env = _snakemake_env(bowtie2_path, samtools_path)
    env["RG_ID"] = "RG_FROM_ENV"

    pairs = [
        ("mid600_f1r2", *_emit_pair(ref_seq, 2_000, 600, "F1R2")),
        ("mid600_f2r1", *_emit_pair(ref_seq, 5_000, 600, "F2R1")),
    ]
    workdir = tmp_path / "migration_run"
    workdir.mkdir()
    raw1 = workdir / "M1_R1.fq.gz"
    raw2 = workdir / "M1_R2.fq.gz"
    with (
        gzip.open(raw1, "wt", encoding="utf-8") as f1,
        gzip.open(raw2, "wt", encoding="utf-8") as f2,
    ):
        for qname, r1_seq, r2_seq in pairs:
            f1.write(f"@{qname}/1\n{r1_seq}\n+\n{QUAL}\n")
            f2.write(f"@{qname}/2\n{r2_seq}\n+\n{QUAL}\n")

    (workdir / "samples.tsv").write_text(
        SAMPLES_HEADER
        + f"M1\tM1_R1.fq.gz\tM1_R2.fq.gz\tPE\tmnase\tH3\tnucleosome\ths\t"
        f"{idx_prefix}\texp1\tH3\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    extra_args = (
        '--max\\\nins 800 --rg "DS:Desc With Spaces" '
        '--rg "FO:$RG_ID" --rg CN:Trailing\\ '
    )
    (workdir / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "mapq: 30",
                "trim: false",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "chipseq_idr: false",
                "blacklist_mode: off",
                "tool_parameters:",
                "  bowtie2:",
                f"    extra_args: {json.dumps(extra_args)}",
                "genome_resources:",
                "  hs:",
                f'    effective_genome_size: "{REF_LEN}"',
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
            str(workdir / "config.yaml"),
            "-j",
            "1",
            "results/M1/02_align/M1.mapq30.bam",
        ],
        env=env,
        cwd=workdir,
    )
    header_and_reads = _run(
        [
            samtools_path,
            "view",
            "-h",
            str(workdir / "results/M1/02_align/M1.mapq30.bam"),
        ],
        env=env,
    ).stdout
    assert "DS:Desc With Spaces" in header_and_reads
    assert "FO:RG_FROM_ENV" in header_and_reads
    assert "CN:Trailing " in header_and_reads
    records = _parse_sam_records(header_and_reads)
    assert (records[("mid600_f1r2", 1)]["flag"] & 0x2) == 0x2
    assert (records[("mid600_f2r1", 1)]["flag"] & 0x2) == 0x2

    with open(
        workdir / "results/.helixweave/effective_config.json", encoding="utf-8"
    ) as handle:
        record = json.load(handle)
    assert record["config"]["tool_parameters"]["bowtie2"]["max_insert"] == 800
    assert (
        record["config"]["tool_parameters"]["bowtie2"]["extra_args"]
        == '--rg "DS:Desc With Spaces" --rg "FO:$RG_ID" --rg CN:Trailing\\ '
    )
    assert record["resolved_max_insert"] == {"M1": 800}

    # Now verify a line-continued conflict fails at Snakemake startup
    conflict_dir = tmp_path / "conflict_run"
    conflict_dir.mkdir()
    shutil.copy(workdir / "samples.tsv", conflict_dir / "samples.tsv")
    conflict_extra = json.dumps("-X\\\n 900")
    (conflict_dir / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "use_control: false",
                "multiqc: false",
                "blacklist_mode: off",
                "tool_parameters:",
                "  bowtie2:",
                "    max_insert: 700",
                f"    extra_args: {conflict_extra}",
                "genome_resources:",
                "  hs:",
                f'    effective_genome_size: "{REF_LEN}"',
                "",
            ]
        ),
        encoding="utf-8",
    )
    failed = _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(conflict_dir / "config.yaml"),
            "-j",
            "1",
            "results/.helixweave/analysis_identity.json",
        ],
        env=env,
        cwd=conflict_dir,
        check=False,
    )
    assert failed.returncode != 0
    assert "tool_parameters.bowtie2.max_insert" in (failed.stdout + failed.stderr)


@pytest.mark.parametrize(
    "rg_extra_args,expected_rg_id",
    [
        ("--rg-id '-X500' --rg 'SM:sample 1'", "ID:-X500"),
        ("--rg-id=-X500 --rg 'SM:sample 1'", "ID:-X500"),
        ("--rg-id '-I5000' --rg 'SM:sample 1'", "ID:-I5000"),
        ("--rg-id=-I5000 --rg 'SM:sample 1'", "ID:-I5000"),
    ],
)
def test_production_bowtie2_rule_preserves_rg_id_values_resembling_insert_flags(
    tmp_path,
    snakemake_path,
    bowtie2_path,
    samtools_path,
    reference_bundle,
    rg_extra_args,
    expected_rg_id,
):
    """Verify on the production bowtie2_align rule that option values such as
    --rg-id '-X500' and --rg-id '-I5000' are preserved into @RG headers while
    max_insert=2000 still admits 600, 1200, and 2000 bp pairs as proper pairs.
    """
    ref_seq = reference_bundle["seq"]
    idx_prefix = reference_bundle["idx_prefix"]
    env = _snakemake_env(bowtie2_path, samtools_path)

    pairs = [
        ("p600", *_emit_pair(ref_seq, 2_000, 600, "F1R2")),
        ("p1200", *_emit_pair(ref_seq, 6_000, 1200, "F1R2")),
        ("p2000", *_emit_pair(ref_seq, 12_000, 2000, "F2R1")),
        ("p2001", *_emit_pair(ref_seq, 18_000, 2001, "F1R2")),
    ]
    workdir = tmp_path / "rg_id_run"
    workdir.mkdir()
    r1 = workdir / "S1_R1.fq.gz"
    r2 = workdir / "S1_R2.fq.gz"
    with (
        gzip.open(r1, "wt", encoding="utf-8") as f1,
        gzip.open(r2, "wt", encoding="utf-8") as f2,
    ):
        for qname, s1, s2 in pairs:
            f1.write(f"@{qname}/1\n{s1}\n+\n{QUAL}\n")
            f2.write(f"@{qname}/2\n{s2}\n+\n{QUAL}\n")

    (workdir / "samples.tsv").write_text(
        SAMPLES_HEADER
        + f"S1\tS1_R1.fq.gz\tS1_R2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t"
        f"{idx_prefix}\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    (workdir / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "mapq: 30",
                "trim: false",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "chipseq_idr: false",
                "blacklist_mode: off",
                "tool_parameters:",
                "  bowtie2:",
                "    max_insert: 2000",
                f"    extra_args: {json.dumps(rg_extra_args)}",
                "genome_resources:",
                "  hs:",
                f'    effective_genome_size: "{REF_LEN}"',
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
            str(workdir / "config.yaml"),
            "-j",
            "1",
            "results/S1/02_align/S1.sorted.bam",
        ],
        env=env,
        cwd=workdir,
    )
    sam_text = _run(
        [
            samtools_path,
            "view",
            "-h",
            str(workdir / "results/S1/02_align/S1.sorted.bam"),
        ],
        env=env,
    ).stdout
    assert f"@RG\t{expected_rg_id}\t" in sam_text
    assert "\tSM:sample 1\n" in sam_text
    expected_rg_tag = f"RG:Z:{expected_rg_id.split(':', 1)[1]}"
    assert expected_rg_tag in sam_text
    records = _parse_sam_records(sam_text)
    for qname in ("p600", "p1200", "p2000"):
        assert (records[(qname, 1)]["flag"] & 0x2) == 0x2
    assert (records[("p2001", 1)]["flag"] & 0x2) == 0


def test_bowtie2_read_length_exceeding_max_insert_boundary(
    tmp_path, bowtie2_path, reference_bundle
):
    """When read length (150 bp) exceeds -X (100), Bowtie2 still marks a 150 bp
    template as a proper pair (0x2) because each 150 bp mate spans the full
    template, whereas a 151 bp template is not marked as a proper pair (0x2).
    """
    ref_seq = reference_bundle["seq"]
    idx_prefix = reference_bundle["idx_prefix"]
    read_len_150 = 150
    qual_150 = "I" * read_len_150

    def emit_150(pos0: int, frag_len: int):
        frag = ref_seq[pos0 : pos0 + frag_len]
        return frag[:read_len_150], _revcomp(frag[-read_len_150:])

    r1_fq = tmp_path / "r1_150.fq"
    r2_fq = tmp_path / "r2_150.fq"
    with (
        open(r1_fq, "w", encoding="utf-8") as f1,
        open(r2_fq, "w", encoding="utf-8") as f2,
    ):
        for qname, pos0, flen in (("len150", 2_000, 150), ("len151", 5_000, 151)):
            s1, s2 = emit_150(pos0, flen)
            f1.write(f"@{qname}/1\n{s1}\n+\n{qual_150}\n")
            f2.write(f"@{qname}/2\n{s2}\n+\n{qual_150}\n")

    sam_out = _run(
        [
            bowtie2_path,
            "--very-sensitive",
            "-X",
            "100",
            "-x",
            str(idx_prefix),
            "-1",
            str(r1_fq),
            "-2",
            str(r2_fq),
        ]
    ).stdout
    records = _parse_sam_records(sam_out)
    assert (records[("len150", 1)]["flag"] & 0x2) == 0x2
    assert abs(records[("len150", 1)]["tlen"]) == 150
    assert (records[("len151", 1)]["flag"] & 0x2) == 0
