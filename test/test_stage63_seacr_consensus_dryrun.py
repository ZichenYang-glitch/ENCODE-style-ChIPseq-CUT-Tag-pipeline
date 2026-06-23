"""Stage 63 SEACR consensus DAG dry-run tests.

Verifies config gating, per-mode targets, final semantics, and that
SEACR IDR never appears.
"""

import os
import subprocess
import sys
import tempfile

from _tool_resolver import resolve_tool


SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")


def main():
    passed = 0
    total = 0

    def check(name, fn):
        nonlocal passed, total
        total += 1
        try:
            if fn():
                print("PASS: %s" % name)
                passed += 1
        except Exception as e:
            print("FAIL: %s\n   Unexpected error: %s" % (name, e))

    try:
        version = subprocess.run(
            [SNAKEMAKE, "--version"], capture_output=True, text=True)
    except OSError:
        version = None

    if version is None or version.returncode != 0:
        print("SKIP: snakemake not available (all tests)")
        return 0

    def _write_files(config_yaml, samples_tsv):
        tmp = tempfile.mkdtemp(prefix="stage63_")
        r1 = os.path.join(tmp, "R1.fq")
        r2 = os.path.join(tmp, "R2.fq")
        open(r1, "w").close()
        open(r2, "w").close()
        samples_fmt = samples_tsv.format(r1=r1, r2=r2)
        samples_path = os.path.join(tmp, "samples.tsv")
        with open(samples_path, "w") as f:
            f.write(samples_fmt)
        config_fmt = config_yaml.format(samples=samples_path)
        config_path = os.path.join(tmp, "config.yaml")
        with open(config_path, "w") as f:
            f.write(config_fmt)
        return tmp, config_path

    def _dry_run(name, config_yaml, samples_tsv, expect_target="",
                 expect_no_target="", expect_final="", expect_no_final=""):
        tmp, config_path = _write_files(config_yaml, samples_tsv)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = os.path.join(tmp, ".cache")
        result = subprocess.run(
            [SNAKEMAKE, "-s", SNAKEFILE, "--configfile", config_path, "-n", "-p"],
            capture_output=True, text=True, timeout=120, cwd=os.getcwd(),
            env=env,
        )
        output = result.stdout + result.stderr
        ok = True
        if result.returncode != 0:
            print("FAIL: %s\n   Snakemake dry-run failed:\n%s" % (name, output[-4000:]))
            ok = False
        expect_targets = (
            [expect_target] if isinstance(expect_target, str) else list(expect_target)
        )
        expect_no_targets = (
            [expect_no_target]
            if isinstance(expect_no_target, str)
            else list(expect_no_target)
        )
        for target in expect_targets:
            if target and target not in output:
                print("FAIL: %s\n   Expected target '%s' not found" % (name, target))
                ok = False
        for target in expect_no_targets:
            if target and target in output:
                print("FAIL: %s\n   Unexpected target '%s' found" % (name, target))
                ok = False
        if expect_final and expect_final not in output:
            print("FAIL: %s\n   Expected final target '%s' not found" % (name, expect_final))
            ok = False
        if expect_no_final and expect_no_final in output:
            print("FAIL: %s\n   Unexpected final target '%s' found" % (name, expect_no_final))
            ok = False
        return ok

    HEADER = (
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
        "\tbowtie2_index\texperiment\tbiological_replicate\n"
    )

    def ms(*rows):
        return HEADER + "\n".join(
            "{id}\t{r1}\t{r2}\t{layout}\t{assay}\tT\t{pm}\ths\tidx\t{exp}\t{br}".format(
                id=id, r1="{r1}", r2="{r2}", layout=layout,
                assay=assay, pm=pm, exp=exp, br=br,
            )
            for id, assay, pm, exp, br, layout in rows
        )

    CT_PE_N2 = ms(
        ("T1", "cuttag", "narrow", "exp_ct", "1", "PE"),
        ("T2", "cuttag", "narrow", "exp_ct", "2", "PE"),
    )
    CT_SE_N2 = ms(
        ("TS1", "cuttag", "narrow", "exp_ctse", "1", "SE"),
        ("TS2", "cuttag", "narrow", "exp_ctse", "2", "SE"),
    )
    CT_PE_1 = ms(
        ("T1", "cuttag", "narrow", "exp_one", "1", "PE"),
    )
    CS_PE_N2 = ms(
        ("C1", "chipseq", "narrow", "exp_cs", "1", "PE"),
        ("C2", "chipseq", "narrow", "exp_cs", "2", "PE"),
    )

    BASE = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: true\n"
    CFG_OFF = BASE + "reproducibility:\n  enabled: false\n  consensus:\n    enabled: true\n"
    CFG_CONS_OFF = BASE + "reproducibility:\n  enabled: true\n  consensus:\n    enabled: false\n"
    CFG_ON = BASE + (
        "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        "cuttag:\n  seacr:\n    enabled: true\n    mode: stringent\n"
    )
    CFG_SEACR_OFF = BASE + (
        "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        "cuttag:\n  seacr:\n    enabled: false\n"
    )

    SEACR_CONSENSUS = "06_reproducibility/consensus/exp_ct.cuttag.seacr.stringent.consensus.bed"
    SEACR_FINAL = "06_reproducibility/final/exp_ct.cuttag.seacr.stringent.replicate_validated.consensus.bed"

    # 1: seacr.enabled=false → no SEACR consensus
    check("1: seacr.enabled=false → no SEACR consensus",
          lambda: _dry_run("1", CFG_SEACR_OFF, CT_PE_N2,
                           expect_no_target="seacr.stringent.consensus.bed"))

    # 2: reproducibility.enabled=false → no SEACR consensus
    check("2: reproducibility.enabled=false → no SEACR consensus",
          lambda: _dry_run("2", CFG_OFF, CT_PE_N2,
                           expect_no_target="seacr.stringent.consensus.bed"))

    # 3: consensus.enabled=false → no SEACR consensus
    check("3: consensus.enabled=false → no SEACR consensus",
          lambda: _dry_run("3", CFG_CONS_OFF, CT_PE_N2,
                           expect_no_target="seacr.stringent.consensus.bed"))

    # 4: CUT&Tag PE, 2 bioreps, SEACR on → consensus + summary + final
    check("4: SEACR consensus BED + summary + final",
          lambda: _dry_run("4", CFG_ON, CT_PE_N2,
                           expect_target=(
                               SEACR_CONSENSUS,
                               "06_reproducibility/consensus/exp_ct.cuttag.seacr.stringent.consensus.summary.tsv",
                           ),
                           expect_final=SEACR_FINAL))

    # 5: 1 biorep → excluded
    check("5: 1 biorep → no SEACR consensus",
          lambda: _dry_run("5", CFG_ON, CT_PE_1,
                           expect_no_target="seacr.stringent.consensus.bed"))

    # 6: SE layout → excluded
    check("6: SE layout → no SEACR consensus",
          lambda: _dry_run("6", CFG_ON, CT_SE_N2,
                           expect_no_target="seacr.stringent.consensus.bed"))

    # 7: ChIP-seq with SEACR enabled → no SEACR consensus
    check("7: ChIP-seq → no SEACR consensus",
          lambda: _dry_run("7", CFG_ON, CS_PE_N2,
                           expect_no_target="seacr.stringent.consensus.bed"))

    # 8: dry-run includes --format bed, --caller seacr, --final-method consensus
    check("8: --format bed --caller seacr --final-method consensus",
          lambda: _dry_run("8", CFG_ON, CT_PE_N2,
                           expect_target=(
                               "--format bed",
                               "--caller seacr",
                               "--final-method consensus",
                           )))

    # 9: No SEACR IDR anywhere
    check("9: no SEACR IDR path/rule/key",
          lambda: _dry_run("9", CFG_ON, CT_PE_N2,
                           expect_no_target=(
                               "seacr.idr",
                               "seacr_idr",
                               "seacr_experimental",
                           )))

    # 10: sample-level SEACR sidecar targets still exist
    check("10: sample-level SEACR sidecar targets remain",
          lambda: _dry_run("10", CFG_ON, CT_PE_N2,
                           expect_target=(
                               "seacr_bedgraph",
                               "seacr_call",
                               "04_peaks_seacr",
                           )))

    # 11: MACS3 consensus not replaced by SEACR (CUT&Tag MACS3 still has targets)
    check("11: MACS3 consensus still present for CUT&Tag",
          lambda: _dry_run("11", CFG_ON, CT_PE_N2,
                           expect_target="cuttag.macs3.narrow.consensus.narrowPeak"))

    # 12: Shell safety — all new shell blocks have set -e -o pipefail
    # (verified by Stage 57 test, not repeated here)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
