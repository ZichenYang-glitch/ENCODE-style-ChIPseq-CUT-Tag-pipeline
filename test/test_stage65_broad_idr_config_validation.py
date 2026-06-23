"""Stage 65 broad IDR config validation tests.

Verifies chipseq_broad_experimental and cuttag_broad_experimental validation.
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

    CS_B2 = ms(
        ("CB1", "chipseq", "broad", "exp_csb", "1", "PE"),
        ("CB2", "chipseq", "broad", "exp_csb", "2", "PE"),
    )
    CS_B1 = ms(
        ("CB1", "chipseq", "broad", "exp_one", "1", "PE"),
    )
    CS_B3 = ms(
        ("CBX1", "chipseq", "broad", "exp_tri", "1", "PE"),
        ("CBX2", "chipseq", "broad", "exp_tri", "2", "PE"),
        ("CBX3", "chipseq", "broad", "exp_tri", "3", "PE"),
    )
    CT_B2 = ms(
        ("TB1", "cuttag", "broad", "exp_ctb", "1", "PE"),
        ("TB2", "cuttag", "broad", "exp_ctb", "2", "PE"),
    )
    CS_N2 = ms(
        ("C1", "chipseq", "narrow", "exp_cs", "1", "PE"),
        ("C2", "chipseq", "narrow", "exp_cs", "2", "PE"),
    )
    CT_N2 = ms(
        ("T1", "cuttag", "narrow", "exp_ct", "1", "PE"),
        ("T2", "cuttag", "narrow", "exp_ct", "2", "PE"),
    )

    BASE = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: true\n"
    CFG_CS_BROAD = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
    )
    CFG_CT_BROAD = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    cuttag_broad_experimental: true\n"
    )
    CFG_BOTH_BROAD = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
        "    cuttag_broad_experimental: true\n"
    )
    CFG_NO_STAGE4B = BASE.replace("stage4b: true", "stage4b: false") + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
    )
    CFG_DISABLED = BASE + (
        "reproducibility:\n  enabled: false\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
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
        with tempfile.TemporaryDirectory(prefix="stage65_") as tmp:
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

    # V1: chipseq broad × 2 → pass
    check("V1: chipseq broad × 2 → pass",
          lambda: run_validator("V1", CFG_CS_BROAD, CS_B2))

    # V2: cuttag broad × 2 → pass
    check("V2: cuttag broad × 2 → pass",
          lambda: run_validator("V2", CFG_CT_BROAD, CT_B2))

    # V3: both flags + both eligible → pass
    check("V3: both flags + both eligible → pass",
          lambda: run_validator("V3", CFG_BOTH_BROAD, join(CS_B2, CT_B2)))

    # V4: chipseq broad + stage4b false → fail
    check("V4: stage4b false → fail",
          lambda: run_validator("V4", CFG_NO_STAGE4B, CS_B2,
                                expect_fail=True, expected_error="stage4b"))

    # V5: chipseq broad + 1 biorep → fail
    check("V5: 1 biorep → fail",
          lambda: run_validator("V5", CFG_CS_BROAD, CS_B1,
                                expect_fail=True, expected_error="exactly 2"))

    # V6: chipseq broad + 3 bioreps → fail
    check("V6: 3 bioreps → fail",
          lambda: run_validator("V6", CFG_CS_BROAD, CS_B3,
                                expect_fail=True, expected_error="exactly 2"))

    # V7: chipseq broad + chipseq narrow only → fail: no eligible
    check("V7: narrow only → fail: no eligible",
          lambda: run_validator("V7", CFG_CS_BROAD, CS_N2,
                                expect_fail=True, expected_error="no eligible"))

    # V8: cuttag broad + cuttag narrow only → fail: no eligible
    check("V8: cuttag broad + narrow only → fail",
          lambda: run_validator("V8", CFG_CT_BROAD, CT_N2,
                                expect_fail=True, expected_error="no eligible"))

    # V9: reproducibility disabled → pass
    check("V9: repro disabled → pass",
          lambda: run_validator("V9", CFG_DISABLED, CS_B2))

    # V10: chipseq broad × 2 + mixed (narrow experiments) → pass
    check("V10: broad + mixed narrow → pass",
          lambda: run_validator("V10", CFG_CS_BROAD, join(CS_B2, CS_N2, CT_N2)))

    # V11: chipseq broad + valid & invalid (3 bioreps) → fail
    check("V11: valid + invalid (3 bioreps) → fail",
          lambda: run_validator("V11", CFG_CS_BROAD, join(CS_B2, CS_B3),
                                expect_fail=True, expected_error="exactly 2"))

    # V12: cuttag broad × 2 + cuttag narrow × 2 → pass
    check("V12: cuttag broad + narrow coexisting → pass",
          lambda: run_validator("V12", CFG_CT_BROAD, join(CT_B2, CT_N2)))

    # V13: chipseq broad + invalid idr.rank → fail
    CFG_BAD_RANK = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
        "idr:\n  threshold: 0.05\n  rank: invalid_rank\n  seed: 42\n"
    )
    check("V13: invalid idr.rank → fail",
          lambda: run_validator("V13", CFG_BAD_RANK, CS_B2,
                                expect_fail=True, expected_error="rank"))

    # V14: chipseq broad + invalid idr.threshold → fail
    CFG_BAD_THRESH = BASE + (
        "reproducibility:\n  enabled: true\n  idr:\n"
        "    chipseq_broad_experimental: true\n"
        "idr:\n  threshold: 1.5\n  rank: p.value\n  seed: 42\n"
    )
    check("V14: invalid idr.threshold → fail",
          lambda: run_validator("V14", CFG_BAD_THRESH, CS_B2,
                                expect_fail=True,
                                expected_error="threshold"))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
