"""Stage 58 mixed experiment IDR validation tests.

Verifies that IDR validation applies per eligible experiment, not globally
across a mixed sample sheet. Non-eligible experiments (CUT&Tag, ATAC when
stage5 is enabled, ChIP-seq broad, MNase) should not block IDR validation
for eligible experiments.
"""

import os
import subprocess
import sys
import tempfile

VALIDATOR = "scripts/validate_samples.py"

BASE_CONFIG = """\
samples: "{samples}"
use_control: false
threads: 1
stage4b: true
"""


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

    def run_validator(name, config_yaml, samples_tsv,
                      expect_fail=False, expected_error=""):
        with tempfile.TemporaryDirectory(prefix="stage58_") as tmp:
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
            if expected_error and expected_error not in result.stderr and expected_error not in result.stdout:
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

    HEADER = (
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
        "\tbowtie2_index\texperiment\tbiological_replicate\n"
    )

    def make_samples(*rows):
        return HEADER + "\n".join(
            "{id}\t{r1}\t{r2}\tPE\t{assay}\tT\t{pm}\ths\tidx\t{exp}\t{br}".format(
                id=id, r1="{r1}", r2="{r2}",
                assay=assay, pm=pm, exp=exp, br=br,
            )
            for id, assay, pm, exp, br in rows
        )

    # CHIPSEQ_NARROW_2: 2 valid ChIP-seq narrow bioreps
    CS_N2 = make_samples(
        ("C1", "chipseq", "narrow", "exp_cs", "1"),
        ("C2", "chipseq", "narrow", "exp_cs", "2"),
    )
    # CHIPSEQ_NARROW_3: 3 bioreps (invalid for IDR)
    CS_N3 = make_samples(
        ("C1", "chipseq", "narrow", "exp_cs", "1"),
        ("C2", "chipseq", "narrow", "exp_cs", "2"),
        ("C3", "chipseq", "narrow", "exp_cs", "3"),
    )
    # CHIPSEQ_BROAD_2: 2 valid ChIP-seq broad bioreps
    CS_B = make_samples(
        ("CB1", "chipseq", "broad", "exp_cb", "1"),
        ("CB2", "chipseq", "broad", "exp_cb", "2"),
    )
    # CUTTAG_NARROW_2
    CT_N = make_samples(
        ("T1", "cuttag", "narrow", "exp_ct", "1"),
        ("T2", "cuttag", "narrow", "exp_ct", "2"),
    )
    # ATAC_NARROW_2
    AT_N = make_samples(
        ("A1", "atac", "narrow", "exp_at", "1"),
        ("A2", "atac", "narrow", "exp_at", "2"),
    )
    # MNase
    MN = make_samples(
        ("M1", "mnase", "nucleosome", "exp_mn", "1"),
    )

    CFG_S5 = BASE_CONFIG + "stage5: true\n"
    CFG_ATAC = BASE_CONFIG + (
        "stage5: false\n"
        "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    )
    CFG_BOTH = BASE_CONFIG + (
        "stage5: true\n"
        "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    )

    def join(*samples):
        """Join sample-table fragments while keeping exactly one header."""
        merged = []
        for idx, sample_text in enumerate(samples):
            lines = sample_text.rstrip("\n").split("\n")
            if idx == 0:
                merged.extend(lines)
            else:
                merged.extend(lines[1:])
        return "\n".join(merged) + "\n"

    # M1: stage5 + ChIP-seq narrow + CUT&Tag → passes
    check("M1: stage5 + ChIP-seq narrow + CUT&Tag",
          lambda: run_validator("M1", CFG_S5, join(CS_N2, CT_N)))

    # M2: stage5 + ChIP-seq narrow + ATAC narrow → passes
    check("M2: stage5 + ChIP-seq narrow + ATAC narrow",
          lambda: run_validator("M2", CFG_S5, join(CS_N2, AT_N)))

    # M3: stage5 + ChIP-seq narrow + MNase → passes
    check("M3: stage5 + ChIP-seq narrow + MNase",
          lambda: run_validator("M3", CFG_S5, join(CS_N2, MN)))

    # M4: stage5 + only CUT&Tag → fails (no eligible)
    check("M4: stage5 + only CUT&Tag → fails",
          lambda: run_validator("M4", CFG_S5, CT_N, expect_fail=True,
                                expected_error="no eligible"))

    # M5: stage5 + only ChIP-seq broad → fails (broad skipped)
    check("M5: stage5 + only ChIP-seq broad → fails",
          lambda: run_validator("M5", CFG_S5, CS_B, expect_fail=True,
                                expected_error="no eligible"))

    # M6: stage5 + ChIP-seq narrow × 3 bioreps → fails
    check("M6: stage5 + ChIP-seq narrow × 3 bioreps → fails",
          lambda: run_validator("M6", CFG_S5, CS_N3, expect_fail=True,
                                expected_error="exactly 2"))

    # M7: atac_narrow + ATAC narrow + ChIP-seq narrow → passes
    check("M7: atac_narrow + ATAC narrow + ChIP-seq narrow",
          lambda: run_validator("M7", CFG_ATAC, join(AT_N, CS_N2)))

    # M8: both stage5 and atac_narrow + both eligible → passes
    CS_N2_B = make_samples(
        ("C1", "chipseq", "narrow", "exp_cs", "1"),
        ("C2", "chipseq", "narrow", "exp_cs", "2"),
    )
    AT_N_B = make_samples(
        ("A1", "atac", "narrow", "exp_at", "1"),
        ("A2", "atac", "narrow", "exp_at", "2"),
    )
    check("M8: stage5 + atac_narrow + both eligible",
          lambda: run_validator("M8", CFG_BOTH, join(CS_N2_B, AT_N_B)))

    # M9: stage5 + one valid ChIP-seq narrow + one invalid (3 bioreps) → fails
    CS_INVALID_3 = make_samples(
        ("CX1", "chipseq", "narrow", "exp_inv", "1"),
        ("CX2", "chipseq", "narrow", "exp_inv", "2"),
        ("CX3", "chipseq", "narrow", "exp_inv", "3"),
    )
    check("M9: stage5 + valid CS-narrow + invalid CS-narrow × 3",
          lambda: run_validator("M9", CFG_S5, join(CS_N2, CS_INVALID_3),
                                expect_fail=True, expected_error="exactly 2"))

    # M10: stage5 + one ChIP-seq narrow + one ChIP-seq broad → passes
    check("M10: stage5 + ChIP-seq narrow + ChIP-seq broad → passes",
          lambda: run_validator("M10", CFG_S5, join(CS_N2, CS_B)))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
