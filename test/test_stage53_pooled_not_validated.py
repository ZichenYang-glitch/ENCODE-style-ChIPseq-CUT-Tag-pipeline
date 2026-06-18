"""Stage 53 pooled-peaks-not-validated contract tests.

Verifies that the reproducibility policy document and output-contract
artifacts (if any) correctly classify pooled peaks as aggregate-signal,
not as replicate-validated peak outputs.

Stage 53 scope: document-level only. No file-existence assertions.
"""

import sys
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POLICY_DOC = os.path.join(REPO_ROOT, "docs", "reproducibility-policy.md")


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

    # 1. Policy doc exists
    def t1():
        ok = os.path.isfile(POLICY_DOC)
        if not ok:
            print("   File not found: %s" % POLICY_DOC)
        return ok
    check("policy doc exists", t1)

    # 2. Policy doc states pooled peaks are not validated
    def t2():
        with open(POLICY_DOC) as fh:
            content = fh.read()
        # Must contain the hard statement
        key_phrase = "Pooled Peaks Are Not Validated Peaks"
        ok = key_phrase.lower() in content.lower()
        if not ok:
            print("   Missing statement: '%s'" % key_phrase)
        return ok
    check("policy doc: pooled peaks not validated", t2)

    # 3. Policy doc states pooled peaks MUST NOT be labeled as validated
    def t3():
        with open(POLICY_DOC) as fh:
            content = fh.read()
        must_not = "MUST NOT be labeled" in content or "MUST NOT be classified" in content
        if not must_not:
            print("   Missing MUST NOT language for pooled peak classification")
        return must_not
    check("policy doc: MUST NOT label pooled as validated", t3)

    # 4. Output path section does not list pooled under 06_reproducibility
    def t4():
        with open(POLICY_DOC) as fh:
            content = fh.read()
        # 06_reproducibility section should not reference 04_peaks/pooled
        # Find the reproducibility output section
        repro_section_start = content.find("06_reproducibility/")
        if repro_section_start == -1:
            print("   Could not find 06_reproducibility/ section")
            return False
        repro_section = content[repro_section_start:repro_section_start + 3000]
        # pooled should not appear as an output under reproducibility
        has_pooled_in_repro = "04_peaks/pooled" in repro_section
        if has_pooled_in_repro:
            print("   Found 04_peaks/pooled reference in reproducibility section")
            return False
        return True
    check("no pooled path in reproducibility output section", t4)

    # 5. Strategy matrix has entry for every mode (indirect check — also
    #    covered by test_stage53_reproducibility_policy_contract.py)
    def t5():
        with open(POLICY_DOC) as fh:
            content = fh.read()
        assays = ["chipseq", "cuttag", "atac"]
        found = []
        for assay in assays:
            if assay in content.lower():
                found.append(assay)
        ok = len(found) == 3
        if not ok:
            print("   Missing assay references: expected %s, found %s" % (assays, found))
        return ok
    check("policy doc references all three assays", t5)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
