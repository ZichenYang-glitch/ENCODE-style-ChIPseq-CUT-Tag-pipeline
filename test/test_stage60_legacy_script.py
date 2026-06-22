"""Stage 60 legacy script cleanup tests.

Verifies:
- scripts/chipseq.sh is a deprecation shim (exits non-zero, mentions Snakemake)
- bash -n passes (valid syntax)
- User-facing docs (README.md) do not recommend the legacy script as active
- Archive references are under docs/archive/ or marked deprecated
"""

import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO_ROOT, "scripts", "chipseq.sh")
README = os.path.join(REPO_ROOT, "README.md")
KNOWN_ISSUES = os.path.join(REPO_ROOT, "KNOWN_ISSUES.md")


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

    # T1: bash -n passes
    def t1():
        result = subprocess.run(
            ["bash", "-n", SCRIPT], capture_output=True, text=True
        )
        return result.returncode == 0

    check("T1: bash -n passes", t1)

    # T2: executing the script exits non-zero
    def t2():
        result = subprocess.run(
            ["bash", SCRIPT], capture_output=True, text=True
        )
        return result.returncode != 0

    check("T2: script exits non-zero", t2)

    # T3: stderr+stdout mentions "deprecated"
    def t3():
        result = subprocess.run(
            ["bash", SCRIPT], capture_output=True, text=True
        )
        combined = (result.stdout + result.stderr).lower()
        return "deprecated" in combined

    check("T3: output mentions deprecated", t3)

    # T4: output mentions Snakemake and workflow/Snakefile
    def t4():
        result = subprocess.run(
            ["bash", SCRIPT], capture_output=True, text=True
        )
        combined = result.stdout + result.stderr
        return "snakemake" in combined.lower() and "workflow/Snakefile" in combined

    check("T4: output mentions Snakemake + workflow/Snakefile", t4)

    # T5: stderr is the primary output channel (error goes to stderr)
    def t5():
        result = subprocess.run(
            ["bash", SCRIPT], capture_output=True, text=True
        )
        return "deprecated" in result.stderr.lower()

    check("T5: deprecation message on stderr", t5)

    # T6: README does not recommend chipseq.sh as active entrypoint
    def t6():
        with open(README) as f:
            content = f.read()

        # Should mention it's deprecated
        deprecated_ref = "deprecated" in content.lower()
        if not deprecated_ref:
            print("   README does not mention deprecation")
            return False

        # Should NOT present it as the recommended entrypoint
        # Check that active/recommended language always points to Snakemake
        lines = content.split("\n")
        for i, line in enumerate(lines):
            low = line.lower()
            if "chipseq.sh" in low and "deprecated" not in low and "archive" not in low:
                # Check surrounding context
                ctx_start = max(0, i - 2)
                ctx_end = min(len(lines), i + 3)
                ctx = "\n".join(lines[ctx_start:ctx_end]).lower()
                if "deprecat" not in ctx and "archive" not in ctx and "legacy" not in ctx:
                    print("   Line %d: chipseq.sh referenced outside deprecated/archive/legacy context" % (i + 1))
                    return False
        return True

    check("T6: README marks chipseq.sh as deprecated", t6)

    # T7: README and KNOWN_ISSUES do not describe the legacy script as supported
    def t7():
        forbidden_phrases = [
            "retained for compatibility",
            "remains available for compatibility",
            "legacy script remains available",
            "scripts/chipseq.sh is retained",
        ]
        for path in (README, KNOWN_ISSUES):
            with open(path) as f:
                content = f.read().lower()
            for phrase in forbidden_phrases:
                if phrase in content:
                    print("   %s contains forbidden phrase: %s" % (path, phrase))
                    return False
        return True

    check("T7: docs do not preserve legacy script as supported", t7)

    # T8: Archive file exists
    def t8():
        archive_path = os.path.join(
            REPO_ROOT, "docs", "archive", "scripts", "chipseq-legacy.sh"
        )
        return os.path.isfile(archive_path)

    check("T8: archive copy exists", t8)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
