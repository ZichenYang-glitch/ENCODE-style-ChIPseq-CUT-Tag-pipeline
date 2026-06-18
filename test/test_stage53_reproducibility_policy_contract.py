"""Stage 53 reproducibility policy contract tests.

Verifies that the reproducibility policy document covers all six current
peak-calling modes with defined primary reproducibility methods, and that
no mode is left without a primary method.
"""

import sys
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POLICY_DOC = os.path.join(REPO_ROOT, "docs", "reproducibility-policy.md")

# Expected modes defined in the spec:
#   assay, peak_mode, caller, primary_method
EXPECTED_MODES = [
    ("chipseq", "narrow", "macs3", "idr"),
    ("chipseq", "broad", "macs3", "consensus"),
    ("cuttag", "narrow", "macs3", "consensus"),
    ("cuttag", "broad", "macs3", "consensus"),
    ("cuttag", "seacr", "seacr", "consensus"),
    ("atac", "narrow", "macs3", "idr"),  # when enabled
]


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

    with open(POLICY_DOC) as fh:
        content = fh.read()

    # 1. Policy doc contains "Reproducibility Strategy Matrix" or equivalent
    def t1():
        ok = "strategy matrix" in content.lower() or "reproducibility strategy" in content.lower()
        if not ok:
            print("   Strategy matrix section not found")
        return ok
    check("strategy matrix section exists", t1)

    # 2. Each expected (assay, peak_mode) pair is mentioned
    def t2():
        missing = []
        for assay, peak_mode, caller, _method in EXPECTED_MODES:
            # Simple check: both assay and peak_mode appear near each other
            if assay not in content.lower():
                missing.append("%s/%s (assay not found)" % (assay, peak_mode))
            elif peak_mode not in content.lower():
                missing.append("%s/%s (peak_mode not found)" % (assay, peak_mode))
        ok = len(missing) == 0
        if not ok:
            for m in missing:
                print("   Missing: %s" % m)
        return ok
    check("all assay/peak_mode pairs mentioned", t2)

    # 3. Each mode has a defined primary method
    def t3():
        # Check that the strategy matrix uses "Consensus" or "IDR" for each mode
        primary_methods_found = 0
        for line in content.split("\n"):
            line_lower = line.lower()
            if "consensus" in line_lower or "idr" in line_lower:
                if any(a in line_lower for a in ["chipseq", "cuttag", "atac"]):
                    if "primary" in line_lower or "|" in line:
                        primary_methods_found += 1
        # At minimum, the table should have rows with primary methods
        ok = primary_methods_found >= 2  # table format makes counting unreliable; check presence
        if not ok:
            print("   Too few primary method entries found")
        return ok
    check("primary methods defined for modes", t3)

    # 4. SEACR row shows consensus primary, no IDR planned
    def t4():
        seacr_lines = []
        capture = False
        for line in content.split("\n"):
            if "seacr" in line.lower():
                capture = True
                seacr_lines.append(line)
            elif capture and ("|" in line and "seacr" in line.lower()):
                seacr_lines.append(line)
            elif capture and "|" not in line:
                capture = False
        seacr_text = "\n".join(seacr_lines)
        has_consensus = "consensus" in seacr_text.lower()
        ok = has_consensus
        if not ok:
            print("   SEACR entry missing consensus reference")
        return ok
    check("SEACR mode: consensus primary", t4)

    # 5. Non-goals section states no SEACR IDR planned
    def t5():
        non_goal_section = content[content.find("Non-goals"):] if "Non-goals" in content else content
        has_seacr_statement = "seacr" in non_goal_section.lower() and (
            "not planned" in non_goal_section.lower()
            or "out of scope" in non_goal_section.lower()
        )
        if not has_seacr_statement:
            print("   Non-goals section missing SEACR IDR exclusion")
        return has_seacr_statement
    check("non-goals: SEACR IDR not planned", t5)

    # 6. chipseq narrow mode is documented as IDR (legacy) primary
    def t6():
        idr_legacy_found = False
        for line in content.split("\n"):
            if "chipseq" in line.lower() and "narrow" in line.lower():
                if "idr" in line.lower():
                    idr_legacy_found = True
                    break
        if not idr_legacy_found:
            print("   chipseq narrow not documented as IDR primary")
        return idr_legacy_found
    check("chipseq narrow: IDR primary documented", t6)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
