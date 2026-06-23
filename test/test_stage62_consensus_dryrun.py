"""Stage 62 consensus DAG dry-run tests.

Verifies config gating, per-mode targets, final semantics, and legacy
invariance for consensus peak DAG integration.
"""

import os
import subprocess
import sys
import tempfile


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

    # Check if snakemake is available. Prefer explicit env var because some
    # systems have a non-executable repo-local "snakemake" path earlier in PATH.
    snakemake = None
    candidates = []
    if os.environ.get("SNAKEMAKE"):
        candidates.append(os.environ["SNAKEMAKE"])
    default_snakemake = os.path.join(
        os.path.expanduser("~"), "miniconda3", "envs", "chipseq", "bin",
        "snakemake",
    )
    candidates.extend([
        default_snakemake,
        "snakemake",
    ])
    for candidate in candidates:
        try:
            r = subprocess.run(
                [candidate, "--version"], capture_output=True, text=True)
            if r.returncode == 0:
                snakemake = candidate
                break
        except OSError:
            pass

    if snakemake is None:
        print("SKIP: snakemake not available (all tests)")
        return 0

    def _write_files(config_yaml, samples_tsv):
        tmp = tempfile.mkdtemp(prefix="stage62_")
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
                 expect_no_target="", expect_final_target="",
                 expect_no_final_target=""):
        tmp, config_path = _write_files(config_yaml, samples_tsv)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = os.path.join(tmp, ".cache")
        result = subprocess.run(
            [snakemake, "-s", "workflow/Snakefile",
             "--configfile", config_path, "-n"],
            capture_output=True, text=True, timeout=120, cwd=os.getcwd(),
            env=env,
        )
        output = result.stdout + result.stderr
        ok = True
        if result.returncode != 0:
            print("FAIL: %s\n   Snakemake dry-run failed:\n%s" % (name, output[-4000:]))
            ok = False
        if expect_target and expect_target not in output:
            print("FAIL: %s\n   Expected target '%s' not found" % (name, expect_target))
            ok = False
        if expect_no_target and expect_no_target in output:
            print("FAIL: %s\n   Unexpected target '%s' found" % (name, expect_no_target))
            ok = False
        if expect_final_target and expect_final_target not in output:
            print("FAIL: %s\n   Expected final target '%s' not found" % (name, expect_final_target))
            ok = False
        if expect_no_final_target and expect_no_final_target in output:
            print("FAIL: %s\n   Unexpected final target '%s' found" % (name, expect_no_final_target))
            ok = False
        return ok

    HEADER = (
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
        "\tbowtie2_index\texperiment\tbiological_replicate\n"
    )

    def ms(*rows):
        return HEADER + "\n".join(
            "{id}\t{r1}\t{r2}\tPE\t{assay}\tT\t{pm}\ths\tidx\t{exp}\t{br}".format(
                id=id, r1="{r1}", r2="{r2}",
                assay=assay, pm=pm, exp=exp, br=br,
            )
            for id, assay, pm, exp, br in rows
        )

    CS_N2 = ms(
        ("C1", "chipseq", "narrow", "exp_cs", "1"),
        ("C2", "chipseq", "narrow", "exp_cs", "2"),
    )
    CS_B2 = ms(
        ("CB1", "chipseq", "broad", "exp_cb", "1"),
        ("CB2", "chipseq", "broad", "exp_cb", "2"),
    )
    CT_N2 = ms(
        ("T1", "cuttag", "narrow", "exp_ct", "1"),
        ("T2", "cuttag", "narrow", "exp_ct", "2"),
    )
    CT_B2 = ms(
        ("TB1", "cuttag", "broad", "exp_ctb", "1"),
        ("TB2", "cuttag", "broad", "exp_ctb", "2"),
    )
    AT_N2 = ms(
        ("A1", "atac", "narrow", "exp_at", "1"),
        ("A2", "atac", "narrow", "exp_at", "2"),
    )
    ONE_BR = ms(
        ("S1", "chipseq", "narrow", "exp_one", "1"),
    )

    BASE = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: true\n"

    CFG_REPRO_OFF = BASE + "reproducibility:\n  enabled: false\n  consensus:\n    enabled: true\n"
    CFG_CONS_OFF = BASE + "reproducibility:\n  enabled: true\n  consensus:\n    enabled: false\n"
    CFG_ON = BASE + "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
    CFG_S5 = CFG_ON + "stage5: true\n"

    # D1: reproducibility disabled → no consensus targets
    check("D1: reproducibility disabled → no consensus",
          lambda: _dry_run("D1", CFG_REPRO_OFF, CS_N2,
                           expect_no_target="06_reproducibility/consensus/"))

    # D2: consensus.enabled false → no consensus targets
    check("D2: consensus.enabled false → no consensus",
          lambda: _dry_run("D2", CFG_CONS_OFF, CS_N2,
                           expect_no_target="06_reproducibility/consensus/"))

    # D3: consensus enabled + 2 chipseq narrow bioreps
    check("D3: chipseq narrow consensus peak + summary",
          lambda: _dry_run("D3", CFG_ON, CS_N2,
                           expect_target="06_reproducibility/consensus/exp_cs.chipseq.macs3.narrow.consensus.narrowPeak"))

    # D4: consensus enabled + 2 chipseq broad → final
    check("D4: chipseq broad consensus final",
          lambda: _dry_run("D4", CFG_ON, CS_B2,
                           expect_final_target="replicate_validated.consensus.broadPeak"))

    # D5: consensus enabled + 2 cuttag narrow → final (IDR not yet implemented)
    check("D5: cuttag narrow consensus final",
          lambda: _dry_run("D5", CFG_ON, CT_N2,
                           expect_final_target="replicate_validated.consensus.narrowPeak"))

    # D6: consensus enabled + 2 cuttag broad → final
    check("D6: cuttag broad consensus final",
          lambda: _dry_run("D6", CFG_ON, CT_B2,
                           expect_final_target="replicate_validated.consensus.broadPeak"))

    # D7: consensus enabled + 2 atac narrow → NO final consensus
    check("D7: atac narrow NO final consensus",
          lambda: _dry_run("D7", CFG_ON, AT_N2,
                           expect_no_final_target="replicate_validated.consensus"))

    # D8: chipseq narrow + stage5: true → legacy IDR unchanged, NO consensus final
    check("D8: chipseq narrow stage5 → IDR final, no consensus final",
          lambda: _dry_run("D8", CFG_S5, CS_N2,
                           expect_target="06_idr/final/conservative.narrowPeak",
                           expect_no_final_target="06_reproducibility/final/exp_cs.chipseq.macs3.narrow.replicate_validated.consensus"))

    # D9: 1 biorep → no consensus targets
    check("D9: 1 biorep → no consensus",
          lambda: _dry_run("D9", CFG_ON, ONE_BR,
                           expect_no_target="06_reproducibility/consensus/"))

    # D10: No SEACR consensus targets
    check("D10: no SEACR consensus targets",
          lambda: _dry_run("D10", CFG_ON, CT_N2,
                           expect_no_target="seacr"))

    # D11: Consensus inputs not from 06_idr/ relaxed peaks
    check("D11: consensus inputs from 06_reproducibility/consensus/biorep_peaks/",
          lambda: _dry_run("D11", CFG_ON, CS_N2,
                           expect_no_target="06_idr/"))

    # D12: CUT&Tag narrow + cuttag_narrow true → consensus final still exists
    CFG_CT_N_IDR = CFG_ON + "  idr:\n    cuttag_narrow: true\n"
    check("D12: cuttag narrow + cuttag_narrow=true → consensus final (IDR not implemented)",
          lambda: _dry_run("D12", CFG_CT_N_IDR, CT_N2,
                           expect_final_target="replicate_validated.consensus.narrowPeak"))

    # D13: CUT&Tag narrow + cuttag_narrow false → consensus final
    CFG_CT_N_NOIDR = CFG_ON + "  idr:\n    cuttag_narrow: false\n"
    check("D13: cuttag narrow + cuttag_narrow=false → consensus final",
          lambda: _dry_run("D13", CFG_CT_N_NOIDR, CT_N2,
                           expect_final_target="replicate_validated.consensus.narrowPeak"))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
