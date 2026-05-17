"""Unit tests for scripts/histone_utils.py."""

import sys
sys.path.insert(0, "scripts")
from histone_utils import classify_histone_target


def main():
    print("Starting histone_utils unit tests\n")
    tests = 0
    passed = 0

    # 1. H3K27me3 + broad -> broad_like, ok
    r = classify_histone_target("H3K27me3", "broad")
    tests += 1
    if r["inferred_histone_class"] == "broad_like" and r["peak_mode_status"] == "ok":
        print("PASS: H3K27me3 + broad -> broad_like, ok")
        passed += 1
    else:
        print("FAIL: H3K27me3 + broad -> %s, %s" % (
            r["inferred_histone_class"], r["peak_mode_status"]))

    # 2. H3K27me3 + narrow -> broad_like, mismatch
    r = classify_histone_target("H3K27me3", "narrow")
    tests += 1
    if r["inferred_histone_class"] == "broad_like" and r["peak_mode_status"] == "mismatch":
        print("PASS: H3K27me3 + narrow -> broad_like, mismatch")
        passed += 1
    else:
        print("FAIL: H3K27me3 + narrow -> %s, %s" % (
            r["inferred_histone_class"], r["peak_mode_status"]))

    # 3. H3K27ac + broad -> context_dependent, ok
    r = classify_histone_target("H3K27ac", "broad")
    tests += 1
    if r["inferred_histone_class"] == "context_dependent" and r["peak_mode_status"] == "ok":
        print("PASS: H3K27ac + broad -> context_dependent, ok")
        passed += 1
    else:
        print("FAIL: H3K27ac + broad -> %s, %s" % (
            r["inferred_histone_class"], r["peak_mode_status"]))

    # 4. H3K27ac + narrow -> context_dependent, ok
    r = classify_histone_target("H3K27ac", "narrow")
    tests += 1
    if r["inferred_histone_class"] == "context_dependent" and r["peak_mode_status"] == "ok":
        print("PASS: H3K27ac + narrow -> context_dependent, ok")
        passed += 1
    else:
        print("FAIL: H3K27ac + narrow -> %s, %s" % (
            r["inferred_histone_class"], r["peak_mode_status"]))

    # 5. H3K4me3 + narrow -> narrow_like, ok
    r = classify_histone_target("H3K4me3", "narrow")
    tests += 1
    if r["inferred_histone_class"] == "narrow_like" and r["peak_mode_status"] == "ok":
        print("PASS: H3K4me3 + narrow -> narrow_like, ok")
        passed += 1
    else:
        print("FAIL: H3K4me3 + narrow -> %s, %s" % (
            r["inferred_histone_class"], r["peak_mode_status"]))

    # 6. CTCF + narrow -> unknown, unknown
    r = classify_histone_target("CTCF", "narrow")
    tests += 1
    if r["inferred_histone_class"] == "unknown" and r["peak_mode_status"] == "unknown":
        print("PASS: CTCF + narrow -> unknown, unknown")
        passed += 1
    else:
        print("FAIL: CTCF + narrow -> %s, %s" % (
            r["inferred_histone_class"], r["peak_mode_status"]))

    # 7. h3k27me3 lowercase -> broad_like
    r = classify_histone_target("h3k27me3", "broad")
    tests += 1
    if r["inferred_histone_class"] == "broad_like":
        print("PASS: h3k27me3 lowercase -> broad_like")
        passed += 1
    else:
        print("FAIL: h3k27me3 lowercase -> %s" % r["inferred_histone_class"])

    # 8. H2A.Z punctuation -> narrow_like
    r = classify_histone_target("H2A.Z", "narrow")
    tests += 1
    if r["inferred_histone_class"] == "narrow_like":
        print("PASS: H2A.Z -> narrow_like")
        passed += 1
    else:
        print("FAIL: H2A.Z -> %s" % r["inferred_histone_class"])

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
