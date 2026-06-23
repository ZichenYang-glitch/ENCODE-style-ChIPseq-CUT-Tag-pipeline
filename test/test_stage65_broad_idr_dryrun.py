"""Stage 65 broad IDR dry-run tests.

Verifies DAG rules, output paths, consensus interaction, and legacy invariance.
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

    CS_B2 = ms(
        ("CB1", "chipseq", "broad", "exp_csb", "1", "PE"),
        ("CB2", "chipseq", "broad", "exp_csb", "2", "PE"),
    )
    CT_B2 = ms(
        ("TB1", "cuttag", "broad", "exp_ctb", "1", "PE"),
        ("TB2", "cuttag", "broad", "exp_ctb", "2", "PE"),
    )
    CS_B1 = ms(
        ("CB1", "chipseq", "broad", "exp_one", "1", "PE"),
    )

    BASE = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: true\n"
    CFG_CS = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
    )
    CFG_CT = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    cuttag_broad_experimental: true\n"
    )
    CFG_OFF = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: false\n"
    )
    CFG_DISABLED = BASE + (
        "reproducibility:\n  enabled: false\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
    )

    def _write_files(config_yaml, samples_tsv):
        tmp = tempfile.mkdtemp(prefix="stage65_")
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
                 expect_no_target="", expect_fail=False, expected_error=""):
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
                print("FAIL: %s\n   Expected dry-run failure" % name)
                ok = False
            if expected_error and expected_error not in output:
                print("FAIL: %s\n   Expected error '%s' not found" % name)
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
        return ok

    # D1: chipseq broad → all 7 rules + idr path
    check("D1: all 7 rules + idr paths (chipseq)",
          lambda: _dry_run("D1", CFG_CS, CS_B2,
                           expect_target=(
                               "broad_idr_macs3_biorep",
                               "broad_idr_true_replicates",
                               "06_reproducibility/idr/",
                           )))

    # D2: cuttag broad → same rules appear
    check("D2: cuttag broad → IDR rules appear",
          lambda: _dry_run("D2", CFG_CT, CT_B2,
                           expect_target="broad_idr_true_replicates"))

    # D3: chipseq_broad_experimental false → no broad IDR
    check("D3: flag false → no broad IDR",
          lambda: _dry_run("D3", CFG_OFF, CS_B2,
                           expect_no_target="broad_idr_"))

    # D4: repro disabled → no broad IDR
    check("D4: repro disabled → no broad IDR",
          lambda: _dry_run("D4", CFG_DISABLED, CS_B2,
                           expect_no_target="broad_idr_"))

    # D5: final IDR broadPeak path present
    check("D5: final IDR broadPeak path",
          lambda: _dry_run("D5", CFG_CS, CS_B2,
                           expect_target="replicate_validated.idr.broadPeak"))

    # D6: consensus peak + summary remain
    check("D6: consensus peak + summary remain",
          lambda: _dry_run("D6", CFG_CS, CS_B2,
                           expect_target=(
                               "chipseq.macs3.broad.consensus.broadPeak",
                               "consensus.summary.tsv",
                           )))

    # D7: consensus final absent for broad IDR experiment
    check("D7: broad consensus final absent",
          lambda: _dry_run("D7", CFG_CS, CS_B2,
                           expect_no_target=(
                               "exp_csb.chipseq.macs3.broad."
                               "replicate_validated.consensus.broadPeak",
                           )))

    # D8: 1 biorep → dry-run fails
    check("D8: 1 biorep → dry-run fails",
          lambda: _dry_run("D8", CFG_CS, CS_B1,
                           expect_fail=True, expected_error="exactly 2"))

    # D9: --input-file-type broadPeak in command
    check("D9: --input-file-type broadPeak",
          lambda: _dry_run("D9", CFG_CS, CS_B2,
                           expect_target="broadPeak"))

    # D10: --broad --broad-cutoff in MACS3 args
    check("D10: --broad --broad-cutoff in MACS3",
          lambda: _dry_run("D10", CFG_CS, CS_B2,
                           expect_target="broad-cutoff"))

    # D11: no SEACR IDR, no legacy 06_idr
    check("D11: no SEACR IDR, no legacy 06_idr",
          lambda: _dry_run("D11", CFG_CS, CS_B2,
                           expect_no_target=(
                               "seacr.idr", "seacr_idr",
                               "06_idr/",
                           )))

    # D12: --assay chipseq in summary command for chipseq broad
    check("D12: --assay chipseq in summary",
          lambda: _dry_run("D12", CFG_CS, CS_B2,
                           expect_target="--assay chipseq"))

    # D13: --assay cuttag in summary command for cuttag broad
    check("D13: --assay cuttag in summary",
          lambda: _dry_run("D13", CFG_CT, CT_B2,
                           expect_target="--assay cuttag"))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
