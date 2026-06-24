"""Stage 53 config validation tests — reproducibility block.

Tests parsing defaults, bad values, null inference, and unknown keys.
All tests use validate_config() directly (no subprocess).
"""

import sys
import os
import warnings

from encode_pipeline.config.validate import validate_config, ValidationError

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


def _report(name, passed):
    if passed:
        print("PASS: %s" % name)
    else:
        print("FAIL: %s" % name)
    return passed


def main():
    passed = 0
    total = 0

    def check(name, fn):
        nonlocal passed, total
        total += 1
        try:
            if fn():
                passed += 1
        except Exception as e:
            print("FAIL: %s\n   Unexpected error: %s" % (name, e))

    # 1. reproducibility absent → defaults to enabled=false
    def t1():
        cfg = validate_config(dict(BASE_CONFIG))
        r = cfg.get("reproducibility", {})
        ok = r.get("enabled") == False
        if not ok:
            print("   Expected enabled=False, got %s" % r)
        return ok
    check("absent → enabled=false", t1)

    # 2. reproducibility.enabled=false → ignores sub-keys
    def t2():
        cfg = validate_config(dict(BASE_CONFIG, reproducibility={
            "enabled": False,
            "consensus": {"min_replicates": 1},  # would error if validated
        }))
        r = cfg["reproducibility"]
        return r["enabled"] == False and "consensus" not in r
    check("enabled=false ignores sub-keys", t2)

    # 3. enabled=true with defaults
    def t3():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg = validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
            }))
        r = cfg["reproducibility"]
        ok = (
            r["enabled"] == True
            and r["consensus"]["min_replicates"] == 2
            and r["consensus"]["reciprocal_overlap"] == 0.5
            and r["consensus"]["enabled"] == True
            and r["idr"]["chipseq_narrow"] is None
            and r["idr"]["atac_narrow"] == False
            and r["idr"]["cuttag_narrow"] == False
        )
        if not ok:
            print("   Got: %s" % r)
        return ok
    check("enabled=true with defaults", t3)

    # 4. min_replicates default (2)
    def t4():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg = validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {},
            }))
        return cfg["reproducibility"]["consensus"]["min_replicates"] == 2
    check("min_replicates default=2", t4)

    # 5. min_replicates: 1 → error
    def t5():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"min_replicates": 1},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "min_replicates" in str(e) and ">= 2" in str(e)
    check("min_replicates=1 rejected", t5)

    # 6. min_replicates string "3" → coerced
    def t6():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg = validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"min_replicates": "3"},
            }))
        return cfg["reproducibility"]["consensus"]["min_replicates"] == 3
    check("min_replicates='3' coerced", t6)

    # 7. min_replicates bool → error
    def t7():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"min_replicates": True},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "min_replicates" in str(e)
    check("min_replicates=true rejected", t7)

    # 8. reciprocal_overlap default (0.5)
    def t8():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg = validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {},
            }))
        return cfg["reproducibility"]["consensus"]["reciprocal_overlap"] == 0.5
    check("reciprocal_overlap default=0.5", t8)

    # 9. reciprocal_overlap: 0 → error (must be > 0)
    def t9():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"reciprocal_overlap": 0},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "reciprocal_overlap" in str(e)
    check("reciprocal_overlap=0 rejected", t9)

    # 10. reciprocal_overlap: 1.5 → error
    def t10():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"reciprocal_overlap": 1.5},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "reciprocal_overlap" in str(e)
    check("reciprocal_overlap=1.5 rejected", t10)

    # 11. reciprocal_overlap: 1.0 → allowed (must be <= 1)
    def t11():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg = validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"reciprocal_overlap": 1.0},
            }))
        return cfg["reproducibility"]["consensus"]["reciprocal_overlap"] == 1.0
    check("reciprocal_overlap=1.0 accepted", t11)

    # 12. consensus unknown key → error
    def t12():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "consensus": {"unknown_key": True},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "unknown key" in str(e)
    check("consensus unknown key rejected", t12)

    # 13. idr atac_narrow non-bool → error
    def t13():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {"atac_narrow": "yes"},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "atac_narrow" in str(e)
    check("atac_narrow='yes' rejected", t13)

    # 14. idr unknown key → error (no seacr_experimental)
    def t14():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "idr": {"seacr_experimental": True},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "unknown key" in str(e) and "seacr_experimental" in str(e)
    check("seacr_experimental rejected (unknown key)", t14)

    # 15. reproducibility unknown top-level key → error
    def t15():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility={
                "enabled": True,
                "unknown_block": {},
            }))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "unknown key" in str(e)
    check("top-level unknown key rejected", t15)

    # 16. chipseq_narrow: null + stage5: true → null (inferred at DAG level)
    def t16():
        cfg = validate_config(dict(BASE_CONFIG, stage5=True, stage4b=True, reproducibility={
            "enabled": True,
        }))
        return cfg["reproducibility"]["idr"]["chipseq_narrow"] is None
    check("chipseq_narrow null preserved (inference deferred)", t16)

    # 17. chipseq_narrow: null + stage5: false → null
    def t17():
        cfg = validate_config(dict(BASE_CONFIG, reproducibility={
            "enabled": True,
        }))
        return cfg["reproducibility"]["idr"]["chipseq_narrow"] is None
    check("chipseq_narrow null when stage5 false", t17)

    # 18. reproducibility not a mapping → error
    def t18():
        try:
            validate_config(dict(BASE_CONFIG, reproducibility="enabled"))
            print("   Expected ValidationError")
            return False
        except ValidationError as e:
            return "mapping" in str(e)
    check("reproducibility string rejected", t18)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
