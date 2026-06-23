"""Stage 64 CUT&Tag narrow IDR config validation tests.

Verifies that cuttag_narrow=true validates correctly per experiment.
"""

import os
import subprocess
import sys
import tempfile

VALIDATOR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts", "validate_samples.py",
)


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
    CT_N3 = ms(
        ("T1", "cuttag", "narrow", "exp_tri", "1", "PE"),
        ("T2", "cuttag", "narrow", "exp_tri", "2", "PE"),
        ("T3", "cuttag", "narrow", "exp_tri", "3", "PE"),
    )
    CT_B2 = ms(
        ("TB1", "cuttag", "broad", "exp_ctb", "1", "PE"),
        ("TB2", "cuttag", "broad", "exp_ctb", "2", "PE"),
    )
    CS_N2 = ms(
        ("C1", "chipseq", "narrow", "exp_cs", "1", "PE"),
        ("C2", "chipseq", "narrow", "exp_cs", "2", "PE"),
    )
    AT_N2 = ms(
        ("A1", "atac", "narrow", "exp_at", "1", "PE"),
        ("A2", "atac", "narrow", "exp_at", "2", "PE"),
    )

    BASE = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: true\n"
    CFG_ON = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
    )
    CFG_NO_STAGE4B = BASE.replace("stage4b: true", "stage4b: false") + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
    )
    CFG_DISABLED = BASE.replace("stage4b: true", "stage4b: false") + (
        "reproducibility:\n  enabled: false\n  idr:\n    cuttag_narrow: true\n"
    )
    CFG_SEACR = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
        "cuttag:\n  seacr:\n    enabled: true\n    mode: stringent\n"
    )
    CFG_BAD_RANK = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
        "idr:\n  threshold: 0.05\n  rank: invalid_rank\n  seed: 42\n"
    )
    CFG_BAD_THRESH = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
        "idr:\n  threshold: 1.5\n  rank: p.value\n  seed: 42\n"
    )

    def join(*samples):
        merged = []
        for idx, s in enumerate(samples):
            lines = s.rstrip("\n").split("\n")
            if idx == 0:
                merged.extend(lines)
            else:
                merged.extend(lines[1:])
        return "\n".join(merged) + "\n"

    def run_validator(name, config_yaml, samples_tsv,
                      expect_fail=False, expected_error=""):
        with tempfile.TemporaryDirectory(prefix="stage64_") as tmp:
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
            result = subprocess.run(
                [sys.executable, VALIDATOR, "--config", config_path],
                capture_output=True, text=True,
            )
        if expect_fail:
            if result.returncode == 0:
                print("FAIL: %s\n   Expected failure, got pass." % name)
                return False
            if (expected_error and
                expected_error not in result.stderr and
                expected_error not in result.stdout):
                print("FAIL: %s\n   Expected error '%s' not found.\n   stderr: %s"
                      % (name, expected_error, result.stderr.strip()[:200]))
                return False
            return True
        else:
            if result.returncode != 0:
                print("FAIL: %s\n   Expected pass, got:\n   %s"
                      % (name, result.stderr.strip()[:200]))
                return False
            return True

    # V1: CUT&Tag PE narrow × 2 → pass
    check("V1: PE narrow × 2 → pass",
          lambda: run_validator("V1", CFG_ON, CT_PE_N2))

    # V2: CUT&Tag SE narrow × 2 → pass (SE supported)
    check("V2: SE narrow × 2 → pass",
          lambda: run_validator("V2", CFG_ON, CT_SE_N2))

    # V3: stage4b false → fail
    check("V3: stage4b false → fail",
          lambda: run_validator("V3", CFG_NO_STAGE4B, CT_PE_N2,
                                expect_fail=True, expected_error="stage4b"))

    # V4: 1 biorep → fail
    check("V4: 1 biorep → fail",
          lambda: run_validator("V4", CFG_ON, CT_N1,
                                expect_fail=True, expected_error="exactly 2"))

    # V5: 3 bioreps → fail
    check("V5: 3 bioreps → fail",
          lambda: run_validator("V5", CFG_ON, CT_N3,
                                expect_fail=True, expected_error="exactly 2"))

    # V6: CUT&Tag broad only → fail
    check("V6: broad only → fail: no eligible",
          lambda: run_validator("V6", CFG_ON, CT_B2,
                                expect_fail=True, expected_error="no eligible"))

    # V7: ChIP-seq narrow only → fail: no eligible
    check("V7: ChIP-seq only → fail: no eligible",
          lambda: run_validator("V7", CFG_ON, CS_N2,
                                expect_fail=True, expected_error="no eligible"))

    # V8: Mixed CUT&Tag narrow + ChIP-seq + ATAC → pass
    check("V8: mixed CUT&Tag narrow + ChIP-seq + ATAC → pass",
          lambda: run_validator("V8", CFG_ON, join(CT_PE_N2, CS_N2, AT_N2)))

    # V9: Mixed CUT&Tag narrow + CUT&Tag broad → pass
    check("V9: mixed CUT&Tag narrow + broad → pass",
          lambda: run_validator("V9", CFG_ON, join(CT_PE_N2, CT_B2)))

    # V10: reproducibility.enabled false → pass (validation skipped)
    check("V10: repro disabled → pass",
          lambda: run_validator("V10", CFG_DISABLED, CT_PE_N2))

    # V11: Mixed valid + invalid (3 bioreps) → fail
    CT_INVALID_N3 = ms(
        ("TX1", "cuttag", "narrow", "exp_inv", "1", "PE"),
        ("TX2", "cuttag", "narrow", "exp_inv", "2", "PE"),
        ("TX3", "cuttag", "narrow", "exp_inv", "3", "PE"),
    )
    check("V11: valid + invalid (3 bioreps) → fail",
          lambda: run_validator("V11", CFG_ON, join(CT_PE_N2, CT_INVALID_N3),
                                expect_fail=True, expected_error="exactly 2"))

    # V12: cuttag_narrow + seacr.enabled → pass
    check("V12: cuttag_narrow + seacr.enabled → pass",
          lambda: run_validator("V12", CFG_SEACR, CT_PE_N2))

    # V13: invalid idr.rank → fail
    check("V13: invalid idr.rank → fail",
          lambda: run_validator("V13", CFG_BAD_RANK, CT_PE_N2,
                                expect_fail=True, expected_error="rank"))

    # V14: invalid idr.threshold → fail
    check("V14: invalid idr.threshold → fail",
          lambda: run_validator("V14", CFG_BAD_THRESH, CT_PE_N2,
                                expect_fail=True,
                                expected_error="threshold"))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
