"""Stage 53 experimental IDR warning tests.

Verifies that experimental *_experimental flags emit informational
warnings when set to true, and that atac_narrow (established, not
experimental) does not emit such warnings.

All tests use validate_config() directly with warnings capture.
"""

import sys
import os
import warnings

from encode_pipeline.config.validate import validate_config

BASE_CONFIG = {
    "samples": "config/samples.tsv",
    "outdir": "results",
    "threads": 8,
    "mapq": 30,
    "binsize": 10,
    "trim": True,
    "use_control": False,
    "multiqc": True,
    "remove_dup": "auto",
    "extend_reads": "auto",
}


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

    # 1. chipseq_broad_experimental: true → warning
    def t1():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {"chipseq_broad_experimental": True},
            }))
            exp_warnings = [
                x for x in w
                if "Experimental IDR" in str(x.message)
                and "chipseq_broad_experimental" in str(x.message)
            ]
            ok = len(exp_warnings) >= 1
            if not ok:
                print("   Expected experimental IDR warning, found %d warnings" % len(w))
            return ok
    check("chipseq_broad_experimental=true warns", t1)

    # 2. cuttag_broad_experimental: true → warning
    def t2():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {"cuttag_broad_experimental": True},
            }))
            exp_warnings = [
                x for x in w
                if "Experimental IDR" in str(x.message)
                and "cuttag_broad_experimental" in str(x.message)
            ]
            ok = len(exp_warnings) >= 1
            if not ok:
                print("   Expected experimental IDR warning, found %d warnings" % len(w))
            return ok
    check("cuttag_broad_experimental=true warns", t2)

    # 3. atac_narrow: true → NO experimental warning
    def t3():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {"atac_narrow": True},
            }))
            exp_warnings = [
                x for x in w
                if "Experimental IDR" in str(x.message)
            ]
            ok = len(exp_warnings) == 0
            if not ok:
                print("   Unexpected experimental warnings for atac_narrow: %s" %
                      [str(x.message) for x in exp_warnings])
            return ok
    check("atac_narrow=true does NOT warn (established mode)", t3)

    # 4. cuttag_narrow: false → no warning (not experimental by default)
    def t4():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {"cuttag_narrow": False},
            }))
            exp_warnings = [
                x for x in w
                if "Experimental IDR" in str(x.message)
            ]
            ok = len(exp_warnings) == 0
            if not ok:
                print("   Unexpected experimental warnings for cuttag_narrow=false")
            return ok
    check("cuttag_narrow=false no warning", t4)

    # 5. experimental flags both true → both warnings
    def t5():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {
                    "chipseq_broad_experimental": True,
                    "cuttag_broad_experimental": True,
                },
            }))
            exp_warnings = [
                x for x in w
                if "Experimental IDR" in str(x.message)
            ]
            ok = len(exp_warnings) >= 2
            if not ok:
                print("   Expected 2 experimental warnings, found %d" % len(exp_warnings))
            return ok
    check("both experimental flags true → two warnings", t5)

    # 6. experimental flags false → no warnings
    def t6():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {
                    "chipseq_broad_experimental": False,
                    "cuttag_broad_experimental": False,
                },
            }))
            exp_warnings = [
                x for x in w
                if "Experimental IDR" in str(x.message)
            ]
            ok = len(exp_warnings) == 0
            if not ok:
                print("   Unexpected experimental warnings when flags are false")
            return ok
    check("experimental flags false → no warnings", t6)

    # 7. No experimental warnings when reproducibility disabled
    def t7():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": False,
                "idr": {"chipseq_broad_experimental": True},
            }))
            all_warnings = [x for x in w]
            ok = len(all_warnings) == 0
            if not ok:
                print("   Got %d warnings with reproducibility disabled: %s" %
                      (len(all_warnings), [str(x.message) for x in all_warnings]))
            return ok
    check("no warnings when reproducibility disabled", t7)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
