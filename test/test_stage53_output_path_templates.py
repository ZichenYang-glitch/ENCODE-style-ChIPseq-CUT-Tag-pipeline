"""Stage 53 output path template contract tests.

Verifies that output path templates documented in the reproducibility
policy follow the prescribed naming conventions. All checks are
document-level only — no file-existence assertions.

Template formats:
  consensus: <exp>.<assay>.<caller>.<peak_mode>.consensus.<suffix>
  final:     <exp>.<assay>.<caller>.<peak_mode>.replicate_validated.<method>.<suffix>
"""

import sys
import os
import re

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POLICY_DOC = os.path.join(REPO_ROOT, "docs", "reproducibility-policy.md")


def extract_templates(content):
    """Extract template-like strings from the policy doc.
    Looks for <exp>.<assay>.<caller>... patterns.
    """
    templates = set()
    # Match strings containing angle-bracket placeholders
    # Pattern: anything with <exp> followed by other fields
    pattern = r"<exp>\.[\w.<>]+"
    for match in re.finditer(pattern, content):
        templates.add(match.group(0))
    return templates


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

    templates = extract_templates(content)

    # 1. At least one consensus template found
    def t1():
        consensus_templates = [t for t in templates if ".consensus." in t]
        ok = len(consensus_templates) >= 1
        if not ok:
            print("   No consensus templates found")
        return ok
    check("consensus templates exist", t1)

    # 2. Consensus templates match <exp>.<assay>.<caller>.<peak_mode>.consensus.<suffix>
    def t2():
        field = r"(?:\w+|<[^>]+>)"
        consensus_pattern = re.compile(
            r"^<exp>\.%s\.%s\.%s\.consensus\.%s(?:\.tsv)?$"
            % (field, field, field, field)
        )
        mismatches = []
        for t in templates:
            if ".consensus." in t and "replicate_validated" not in t:
                if not consensus_pattern.match(t):
                    mismatches.append(t)
        ok = len(mismatches) == 0
        if not ok:
            for m in mismatches:
                print("   Mismatch: %s" % m)
        return ok
    check("consensus templates follow <exp>.<assay>.<caller>.<peak_mode>.consensus.<suffix>", t2)

    # 3. At least one replicate_validated template found
    def t3():
        rv_templates = [t for t in templates if "replicate_validated" in t]
        ok = len(rv_templates) >= 1
        if not ok:
            print("   No replicate_validated templates found")
        return ok
    check("replicate_validated templates exist", t3)

    # 4. Final templates match <exp>.<assay>.<caller>.<peak_mode>.replicate_validated.<method>.<suffix>
    def t4():
        field = r"(?:\w+|<[^>]+>)"
        rv_pattern = re.compile(
            r"^<exp>\.%s\.%s\.%s\.replicate_validated\.%s\.%s$"
            % (field, field, field, field, field)
        )
        mismatches = []
        for t in templates:
            if "replicate_validated" in t:
                if not rv_pattern.match(t):
                    mismatches.append(t)
        ok = len(mismatches) == 0
        if not ok:
            for m in mismatches:
                print("   Mismatch: %s" % m)
        return ok
    check("final templates follow <exp>.<assay>.<caller>.<peak_mode>.replicate_validated.<method>.<suffix>", t4)

    # 5. No SEACR IDR output template exists
    def t5():
        seacr_idr = [t for t in templates if "seacr" in t.lower() and "idr" in t.lower()]
        ok = len(seacr_idr) == 0
        if not ok:
            print("   Found SEACR+IDR template: %s" % seacr_idr)
        return ok
    check("no SEACR IDR template", t5)

    # 6. Legacy 06_idr/final/ paths documented as unchanged
    def t6():
        has_legacy_reference = "06_idr/final/" in content
        if not has_legacy_reference:
            print("   Legacy 06_idr/final/ path not documented")
        return has_legacy_reference
    check("legacy 06_idr/final/ path documented", t6)

    # 7. New 06_reproducibility namespace documented
    def t7():
        has_repro = "06_reproducibility/" in content
        if not has_repro:
            print("   06_reproducibility/ namespace not documented")
        return has_repro
    check("06_reproducibility/ namespace documented", t7)

    # 8. Consensus summary TSV template documented
    def t8():
        has_summary = ".consensus.summary.tsv" in content or "consensus.summary" in content
        if not has_summary:
            print("   Consensus summary TSV template not found")
        return has_summary
    check("consensus summary TSV documented", t8)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
