"""Stage 55 config validation tests — reproducibility.idr.atac_narrow.

Tests config validation for ATAC narrow IDR: gating, assay/peak_mode
requirements, biological replicate requirements, mixed-run permissiveness,
and idr settings preservation.
"""

import subprocess
import sys
import os

VALIDATOR = "scripts/validate_samples.py"

BASE_CONFIG = """\
samples: "{samples}"
use_control: false
threads: 1
stage4b: true
"""

ATAC_SAMPLES_2BIOREP = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "A1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t1\n"
    "A2\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t2\n"
)

ATAC_SAMPLES_3BIOREP = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "A1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t1\n"
    "A2\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t2\n"
    "A3\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t3\n"
)

ATAC_SAMPLES_1BIOREP = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "A1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t1\n"
)

ATAC_SAMPLES_BROAD = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "AB1\t{r1}\t{r2}\tPE\tatac\tT\tbroad\ths\tidx\texp1\t1\n"
    "AB2\t{r1}\t{r2}\tPE\tatac\tT\tbroad\ths\tidx\texp1\t2\n"
)

MIXED_SAMPLES = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "A1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t1\n"
    "A2\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t2\n"
    "C1\t{r1}\t{r2}\tPE\tchipseq\tT\tnarrow\ths\tidx\texp2\t1\n"
    "C2\t{r1}\t{r2}\tPE\tchipseq\tT\tnarrow\ths\tidx\texp2\t2\n"
)

MIXED_CHIPSEQ_BROAD_ATAC_NARROW = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "AN1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t1\n"
    "AN2\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\texp1\t2\n"
    "CB1\t{r1}\t{r2}\tPE\tchipseq\tT\tbroad\ths\tidx\texp2\t1\n"
    "CB2\t{r1}\t{r2}\tPE\tchipseq\tT\tbroad\ths\tidx\texp2\t2\n"
)


def main():
    import tempfile

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

    def run_validator(name, config_yaml, samples_tsv, expect_fail=False, expected_error=""):
        with tempfile.TemporaryDirectory(prefix="stage55_cfg_") as tmp:
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

    # V1: atac_narrow + enabled + 2 ATAC narrow bioreps → passes
    cfg_v1 = BASE_CONFIG + "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    check("V1: atac_narrow=true + 2 ATAC narrow bioreps", lambda: run_validator("test",
        cfg_v1, ATAC_SAMPLES_2BIOREP))

    # V2: atac_narrow + enabled=false → passes config validation
    cfg_v2 = BASE_CONFIG + "reproducibility:\n  enabled: false\n  idr:\n    atac_narrow: true\n"
    check("V2: atac_narrow=true + enabled=false → passes", lambda: run_validator("test",
        cfg_v2, ATAC_SAMPLES_2BIOREP))

    # V3: atac_narrow + stage4b=false → fails
    cfg_v3 = "samples: \"{samples}\"\nuse_control: false\nthreads: 1\nstage4b: false\n"
    cfg_v3 += "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    check("V3: atac_narrow + stage4b=false → fails", lambda: run_validator("test",
        cfg_v3, ATAC_SAMPLES_2BIOREP, expect_fail=True, expected_error="stage4b"))

    # V4a: ATAC broad remains rejected by baseline assay policy
    check("V4a: ATAC broad only + atac_narrow → fails", lambda: run_validator("test",
        cfg_v1, ATAC_SAMPLES_BROAD, expect_fail=True, expected_error="atac currently supports"))

    # V4b: non-ATAC broad + ATAC narrow eligible → passes; non-ATAC is skipped
    check("V4b: ChIP-seq broad + eligible ATAC narrow → passes", lambda: run_validator("test",
        cfg_v1, MIXED_CHIPSEQ_BROAD_ATAC_NARROW))

    # V5: ATAC narrow + 3 bioreps → fails
    check("V5: ATAC narrow + 3 bioreps → fails", lambda: run_validator("test",
        cfg_v1, ATAC_SAMPLES_3BIOREP, expect_fail=True, expected_error="exactly 2"))

    # V6: ATAC narrow + 1 biorep → fails
    check("V6: ATAC narrow + 1 biorep → fails", lambda: run_validator("test",
        cfg_v1, ATAC_SAMPLES_1BIOREP, expect_fail=True, expected_error="exactly 2"))

    # V7: atac_narrow + no ATAC experiments → fails
    NO_ATAC = (
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
        "\tbowtie2_index\texperiment\tbiological_replicate\n"
        "C1\t{r1}\t{r2}\tPE\tchipseq\tT\tnarrow\ths\tidx\texp1\t1\n"
        "C2\t{r1}\t{r2}\tPE\tchipseq\tT\tnarrow\ths\tidx\texp1\t2\n"
    )
    check("V7: no ATAC experiments → fails", lambda: run_validator("test",
        cfg_v1, NO_ATAC, expect_fail=True, expected_error="no eligible"))

    # V8: stage5=false + mixed ATAC narrow + ChIP-seq narrow → passes
    cfg_v8 = BASE_CONFIG + "stage5: false\n"
    cfg_v8 += "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    check("V8: stage5=false + mixed ATAC+ChIP-seq → passes", lambda: run_validator("test",
        cfg_v8, MIXED_SAMPLES))

    # V9a: idr settings preserved when ATAC IDR enabled + stage5=false
    # (Verified by checking the validated idr dict — indirect via dry-run)
    # For config validation, we just verify that atac_narrow + stage5=false
    # doesn't reject valid idr settings
    cfg_v9a = BASE_CONFIG + "stage5: false\n"
    cfg_v9a += "idr:\n  threshold: 0.01\n  rank: p.value\n  seed: 99\n"
    cfg_v9a += "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    check("V9a: idr settings with atac_narrow + stage5=false → passes", lambda: run_validator("test",
        cfg_v9a, ATAC_SAMPLES_2BIOREP))

    # V9b: invalid idr.rank rejected
    cfg_v9b = BASE_CONFIG + "stage5: false\n"
    cfg_v9b += "idr:\n  rank: invalid_rank\n"
    cfg_v9b += "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    check("V9b: invalid idr.rank → fails", lambda: run_validator("test",
        cfg_v9b, ATAC_SAMPLES_2BIOREP, expect_fail=True, expected_error="rank"))

    # V9c: invalid idr.threshold rejected
    cfg_v9c = BASE_CONFIG + "stage5: false\n"
    cfg_v9c += "idr:\n  threshold: 2.0\n"
    cfg_v9c += "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
    check("V9c: invalid idr.threshold → fails", lambda: run_validator("test",
        cfg_v9c, ATAC_SAMPLES_2BIOREP, expect_fail=True, expected_error="threshold"))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
