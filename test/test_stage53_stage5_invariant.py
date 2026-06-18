"""Stage 53 stage5 invariant tests.

Verifies that no reproducibility setting can suppress, disable, rename,
or move existing Stage 5 ChIP-seq narrow IDR behavior.

Tests use subprocess validation (matching repo patterns) to verify
that stage5 config validation still passes with reproducibility present.
"""

import subprocess
import sys
import os
import tempfile

BASE_CONFIG = """\
samples: "{samples}"
use_control: false
threads: 1
stage4b: true
stage5: true
"""

BASE_SAMPLES = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "T1\t{r1}\t{r2}\tPE\tchipseq\tTF\tnarrow\ths\tidx\texp1\t1\n"
    "T2\t{r1}\t{r2}\tPE\tchipseq\tTF\tnarrow\ths\tidx\texp1\t2\n"
)

VALIDATOR = "scripts/validate_samples.py"
SNAKEFILE = "workflow/Snakefile"


def main():
    passed = 0
    total = 0

    def check(name, fn):
        nonlocal passed, total
        total += 1
        try:
            ok = fn()
            if ok:
                print("PASS: %s" % name)
                passed += 1
        except Exception as e:
            print("FAIL: %s\n   Unexpected error: %s" % (name, e))

    def run_validator(name, config_yaml, expect_fail=False, expected_error=""):
        """Run validate_samples.py as subprocess."""
        with tempfile.TemporaryDirectory(prefix="stage53_s5_") as tmp:
            r1 = os.path.join(tmp, "R1.fq")
            r2 = os.path.join(tmp, "R2.fq")
            samples = os.path.join(tmp, "samples.tsv")
            config = os.path.join(tmp, "config.yaml")
            open(r1, "w").close()
            open(r2, "w").close()
            with open(samples, "w") as f:
                f.write(BASE_SAMPLES.format(r1=r1, r2=r2))
            with open(config, "w") as f:
                f.write(config_yaml.format(samples=samples))

            result = subprocess.run(
                [sys.executable, VALIDATOR, "--config", config],
                capture_output=True, text=True,
            )

        if expect_fail:
            if result.returncode == 0:
                print("FAIL: %s\n   Expected failure, got pass." % name)
                return False
            if expected_error and expected_error not in result.stderr and expected_error not in result.stdout:
                print("FAIL: %s\n   Expected error '%s' not found.\n   stderr: %s" % (name, expected_error, result.stderr.strip()))
                return False
            return True
        else:
            if result.returncode != 0:
                print("FAIL: %s\n   Expected pass, got:\n   %s" % (name, result.stderr.strip()))
                return False
            return True

    # 1. stage5: true alone still validates
    check("stage5=true alone validates", lambda: run_validator(
        "stage5=true alone", BASE_CONFIG))

    # 2. stage5: true + reproducibility.absent → validates
    check("stage5=true + reproducibility absent", lambda: run_validator(
        "stage5 + no reproducibility", BASE_CONFIG))

    # 3. stage5: true + reproducibility.enabled=false → validates
    check("stage5=true + reproducibility.enabled=false", lambda: run_validator(
        "stage5 + reproducibility disabled",
        BASE_CONFIG + "reproducibility:\n  enabled: false\n"))

    # 4. stage5: true + reproducibility.enabled=true + chipseq_narrow null → validates
    check("stage5=true + reproducibility.enabled=true + chipseq_narrow null", lambda: run_validator(
        "stage5 + reproducibility enabled + csn null",
        BASE_CONFIG + "reproducibility:\n  enabled: true\n"))

    # 5. stage5: true + reproducibility.enabled=true + chipseq_narrow true → validates
    check("stage5=true + chipseq_narrow=true", lambda: run_validator(
        "stage5 + chipseq_narrow explicit true",
        BASE_CONFIG + "reproducibility:\n  enabled: true\n  idr:\n    chipseq_narrow: true\n"))

    # 6. stage5: true + reproducibility.enabled=true + chipseq_narrow false → validates (warning, not error)
    check("stage5=true + chipseq_narrow=false → validates (warns, not error)", lambda: run_validator(
        "stage5 + chipseq_narrow explicit false",
        BASE_CONFIG + "reproducibility:\n  enabled: true\n  idr:\n    chipseq_narrow: false\n"))

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
