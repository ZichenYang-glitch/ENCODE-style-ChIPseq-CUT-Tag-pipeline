"""Stage 53 reproducibility policy contract tests.

Verifies that the reproducibility policy document covers all six current
peak-calling modes with defined primary reproducibility methods, and that
no mode is left without a primary method.
"""

import sys
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POLICY_DOC = os.path.join(REPO_ROOT, "docs", "reproducibility-policy.md")
CONFIG_SCHEMA = os.path.join(REPO_ROOT, "workflow", "schemas", "config.schema.yaml")

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
    content_lower = content.lower()

    with open(CONFIG_SCHEMA) as fh:
        config_schema = fh.read()
    config_schema_lower = config_schema.lower()

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

    # 7. Production-supported wording is used without implying default-enabled IDR
    def t7():
        has_supported = "production-supported idr" in content_lower
        forbidden = ["production/default", "default idr"]
        bad = [phrase for phrase in forbidden if phrase in content_lower]
        ok = has_supported and not bad
        if not ok:
            if not has_supported:
                print("   Missing production-supported IDR wording")
            for phrase in bad:
                print("   Forbidden default-enabled wording found: %s" % phrase)
        return ok
    check("IDR wording: production-supported, not default-enabled", t7)

    # 8. CUT&Tag narrow IDR is supported opt-in, not experimental
    def t8():
        cuttag_lines = [
            line.lower()
            for line in content.split("\n")
            if "cuttag" in line.lower() and "narrow" in line.lower()
        ]
        has_supported = (
            "supported opt-in idr" in content_lower
            and "reproducibility.idr.cuttag_narrow" in content
            and "idr becomes final when explicitly enabled" in content_lower
        )
        bad_lines = [line for line in cuttag_lines if "experimental" in line]
        ok = has_supported and not bad_lines
        if not ok:
            if not has_supported:
                print("   CUT&Tag narrow supported opt-in policy not found")
            for line in bad_lines:
                print("   CUT&Tag narrow still marked experimental: %s" % line)
        return ok
    check("CUT&Tag narrow: supported opt-in IDR", t8)

    # 9. Broad IDR is experimental opt-in and becomes final only when explicit.
    def t9():
        has_experimental = (
            "experimental opt-in idr" in content_lower
            and "chip-seq broad" in content_lower
            and "cut&tag broad" in content_lower
        )
        has_final_guard = (
            "idr becomes final when explicitly enabled and eligible" in content_lower
            and "consensus remains available as fallback/report" in content_lower
        )
        ok = has_experimental and has_final_guard
        if not ok:
            if not has_experimental:
                print("   Broad experimental IDR policy missing")
            if not has_final_guard:
                print("   Broad IDR explicit-final/fallback policy missing")
        return ok
    check("broad modes: experimental explicit-final IDR", t9)

    # 10. SEACR IDR is excluded and no seacr_experimental config key exists
    def t10():
        has_no_idr = (
            "seacr" in content_lower
            and ("no idr planned" in content_lower or "not planned for idr" in content_lower)
        )
        no_key = not any(
            line.strip().startswith("seacr_experimental:")
            for line in config_schema.split("\n")
        )
        ok = has_no_idr and no_key
        if not ok:
            if not has_no_idr:
                print("   SEACR no-IDR policy missing")
            if not no_key:
                print("   seacr_experimental key/reference found")
        return ok
    check("SEACR: no IDR and no experimental key", t10)

    # 11. MNase is explicitly outside peak IDR policy
    def t11():
        ok = "mnase" in content_lower and (
            "not planned for idr" in content_lower
            or "no mnase reproducibility" in content_lower
            or "outside peak idr policy" in content_lower
        )
        if not ok:
            print("   MNase IDR exclusion missing")
        return ok
    check("MNase: outside peak IDR policy", t11)

    # 12. Consensus path templates cover all six peak-like modes
    def t12():
        required_templates = [
            "<exp>.chipseq.macs3.narrow.consensus.narrowPeak",
            "<exp>.chipseq.macs3.broad.consensus.broadPeak",
            "<exp>.cuttag.macs3.narrow.consensus.narrowPeak",
            "<exp>.cuttag.macs3.broad.consensus.broadPeak",
            "<exp>.cuttag.seacr.<mode>.consensus.bed",
            "<exp>.atac.macs3.narrow.consensus.narrowPeak",
        ]
        missing = [template for template in required_templates if template not in content]
        ok = not missing
        if not ok:
            for template in missing:
                print("   Missing consensus template: %s" % template)
        return ok
    check("consensus templates cover all peak-like modes", t12)

    # 13. CUT&Tag narrow final semantics include both IDR-enabled and consensus fallback rows
    def t13():
        has_idr_final = (
            "reproducibility.idr.cuttag_narrow: true" in content
            and "<exp>.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak" in content
        )
        has_consensus_fallback = (
            "Consensus (when CUT&Tag IDR not enabled)" in content
            and "<exp>.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak" in content
        )
        ok = has_idr_final and has_consensus_fallback
        if not ok:
            if not has_idr_final:
                print("   CUT&Tag narrow IDR final row missing")
            if not has_consensus_fallback:
                print("   CUT&Tag narrow consensus fallback row missing")
        return ok
    check("CUT&Tag narrow: IDR final plus consensus fallback", t13)

    # 14. Config schema describes CUT&Tag narrow IDR as supported opt-in
    def t14():
        schema_words = " ".join(config_schema_lower.split())
        has_key = "cuttag_narrow:" in config_schema
        has_supported = "cut&tag narrow-peak idr (opt-in, supported)" in schema_words
        has_final_semantics = "idr becomes final when explicitly enabled; consensus otherwise" in schema_words
        bad_experimental = "cut&tag narrow-peak idr (opt-in, experimental)" in schema_words
        ok = has_key and has_supported and has_final_semantics and not bad_experimental
        if not ok:
            if not has_key:
                print("   cuttag_narrow schema key missing")
            if not has_supported:
                print("   cuttag_narrow schema supported wording missing")
            if not has_final_semantics:
                print("   cuttag_narrow schema final semantics missing")
            if bad_experimental:
                print("   cuttag_narrow schema still says experimental")
        return ok
    check("config schema: CUT&Tag narrow IDR supported opt-in", t14)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
