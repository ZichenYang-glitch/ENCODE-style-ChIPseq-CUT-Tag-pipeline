"""Stage 64 CUT&Tag narrow IDR dry-run tests.

Verifies DAG correctness: rule presence, output paths, consensus interaction,
and that legacy/ATAC IDR remain unchanged.
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
    CT_N1 = ms(
        ("T1", "cuttag", "narrow", "exp_one", "1", "PE"),
    )
    CS_N2 = ms(
        ("C1", "chipseq", "narrow", "exp_cs", "1", "PE"),
        ("C2", "chipseq", "narrow", "exp_cs", "2", "PE"),
    )
    AT_N2 = ms(
        ("A1", "atac", "narrow", "exp_at", "1", "PE"),
        ("A2", "atac", "narrow", "exp_at", "2", "PE"),
    )

    def join(*samples):
        merged = []
        for idx, sample_text in enumerate(samples):
            lines = sample_text.rstrip("\n").split("\n")
            if idx == 0:
                merged.extend(lines)
            else:
                merged.extend(lines[1:])
        return "\n".join(merged) + "\n"

    MIXED = join(CT_PE_N2, CS_N2, AT_N2)

    BASE = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: true\n"
    CFG_ON = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
    )
    CFG_OFF = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: false\n"
    )
    CFG_DISABLED = BASE + (
        "reproducibility:\n  enabled: false\n  idr:\n    cuttag_narrow: true\n"
    )
    CFG_CONS_OFF = CFG_ON + "  consensus:\n    enabled: false\n"
    CFG_S5 = CFG_ON + "stage5: true\n"

    def _write_files(config_yaml, samples_tsv):
        tmp = tempfile.mkdtemp(prefix="stage64_")
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

    def _items(value):
        if not value:
            return []
        if isinstance(value, str):
            return [value]
        return list(value)

    def _dry_run(name, config_yaml, samples_tsv, expect_target="",
                 expect_no_target="", expect_final="", expect_no_final="",
                 expect_fail=False, expected_error=""):
        tmp, config_path = _write_files(config_yaml, samples_tsv)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = os.path.join(tmp, ".cache")
        result = subprocess.run(
            [SNAKEMAKE, "-s", SNAKEFILE, "--configfile", config_path,
             "-n", "-p"],
            capture_output=True, text=True, timeout=120, cwd=os.getcwd(),
            env=env,
        )
        output = result.stdout + result.stderr
        ok = True
        if expect_fail:
            if result.returncode == 0:
                print("FAIL: %s\n   Expected dry-run failure, got success" % name)
                ok = False
            if expected_error and expected_error not in output:
                print("FAIL: %s\n   Expected error '%s' not found"
                      % (name, expected_error))
                ok = False
            return ok

        if result.returncode != 0:
            print("FAIL: %s\n   Snakemake dry-run failed (rc=%d):\n%s"
                  % (name, result.returncode, output[-3000:]))
            ok = False

        for t in _items(expect_target):
            if t and t not in output:
                print("FAIL: %s\n   Expected target '%s' not found" % (name, t))
                ok = False
        for t in _items(expect_no_target):
            if t and t in output:
                print("FAIL: %s\n   Unexpected target '%s' found" % (name, t))
                ok = False
        for t in _items(expect_final):
            if t and t not in output:
                print("FAIL: %s\n   Expected final '%s' not found" % (name, t))
                ok = False
        for t in _items(expect_no_final):
            if t and t in output:
                print("FAIL: %s\n   Unexpected final '%s' found" % (name, t))
                ok = False
        return ok

    # D1: CUT&Tag PE narrow × 2 → all 7 rules + idr path
    check("D1: all 7 rules + idr paths",
          lambda: _dry_run("D1", CFG_ON, CT_PE_N2,
                           expect_target=(
                               "cuttag_macs3_idr_biorep",
                               "cuttag_idr_true_replicates",
                               "cuttag_split_pseudoreps",
                               "cuttag_macs3_idr_pseudorep",
                               "cuttag_idr_self_pseudoreps",
                               "cuttag_idr_pooled_pseudoreps",
                               "cuttag_idr_summary",
                               "06_reproducibility/idr/",
                           )))

    # D2: SE layout → same rules appear
    check("D2: SE layout → IDR rules appear",
          lambda: _dry_run("D2", CFG_ON, CT_SE_N2,
                           expect_target="cuttag_idr_true_replicates"))

    # D3: Mixed CUT&Tag + ChIP-seq + ATAC, stage5 off — independent
    check("D3: mixed stage5 off → cuttag IDR only",
          lambda: _dry_run("D3", CFG_ON, MIXED,
                           expect_target="cuttag_idr_true_replicates",
                           expect_no_target="06_idr/"))

    # D4: cuttag_narrow false → no cuttag IDR
    check("D4: cuttag_narrow false → no cuttag IDR",
          lambda: _dry_run("D4", CFG_OFF, CT_PE_N2,
                           expect_no_target="cuttag_idr_"))

    # D5: repro disabled → no cuttag IDR
    check("D5: repro disabled → no cuttag IDR",
          lambda: _dry_run("D5", CFG_DISABLED, CT_PE_N2,
                           expect_no_target="cuttag_idr_"))

    # D6: final IDR path present
    check("D6: final IDR path",
          lambda: _dry_run("D6", CFG_ON, CT_PE_N2,
                           expect_final="replicate_validated.idr.narrowPeak"))

    # D7: consensus peak + summary still present
    check("D7: consensus peak + summary remain",
          lambda: _dry_run("D7", CFG_ON, CT_PE_N2,
                           expect_target=(
                               "cuttag.macs3.narrow.consensus.narrowPeak",
                               "cuttag.macs3.narrow.consensus.summary.tsv",
                           )))

    # D8: exact consensus final path absent for IDR experiment
    check("D8: consensus final absent for IDR experiment",
          lambda: _dry_run("D8", CFG_ON, CT_PE_N2,
                           expect_no_final=(
                               "exp_ct.cuttag.macs3.narrow."
                               "replicate_validated.consensus.narrowPeak",
                           )))

    # D9: 1 biorep → dry-run fails (validation error)
    check("D9: 1 biorep → dry-run fails",
          lambda: _dry_run("D9", CFG_ON, CT_N1,
                           expect_fail=True, expected_error="exactly 2"))

    # D10: --assay cuttag, no atac_idr_ rules
    check("D10: --assay cuttag in command, no atac_idr_",
          lambda: _dry_run("D10", CFG_ON, CT_PE_N2,
                           expect_target="--assay cuttag",
                           expect_no_target="atac_idr_"))

    # D11: no 06_idr legacy path, no SEACR IDR
    check("D11: no legacy 06_idr, no SEACR IDR",
          lambda: _dry_run("D11", CFG_ON, CT_PE_N2,
                           expect_no_target=("06_idr/", "seacr.idr",
                                             "seacr_idr", "seacr_experimental")))

    # D12: consensus.enabled false → IDR still present, consensus absent
    check("D12: consensus.enabled false → IDR still present",
          lambda: _dry_run("D12", CFG_CONS_OFF, CT_PE_N2,
                           expect_target="cuttag_idr_true_replicates",
                           expect_no_target="consensus/"))

    # D13: mixed with legacy stage5 true preserves ChIP-seq IDR
    check("D13: mixed stage5 true → legacy + cuttag IDR coexist",
          lambda: _dry_run("D13", CFG_S5, MIXED,
                           expect_target=(
                               "cuttag_idr_true_replicates",
                               "exp_cs/06_idr/",
                           )))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
