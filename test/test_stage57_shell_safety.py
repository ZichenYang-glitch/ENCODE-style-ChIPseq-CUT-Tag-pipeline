"""Stage 57 shell safety test — enforce set -e -o pipefail in all shell: blocks.

Scans workflow/Snakefile and workflow/rules/*.smk for shell: blocks.
Each block must start with set -e -o pipefail or have an explicit
exemption marker: # no pipefail: <reason>

Exits 0 if all shell blocks pass, non-zero otherwise.
"""

import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RULES_DIR = os.path.join(REPO_ROOT, "workflow", "rules")
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")

GUARD_REQUIRED = "set -e -o pipefail"
EXEMPTION_PATTERN = re.compile(r"#\s*no\s+pipefail\s*:\s*(.*)")


def find_smk_files():
    """Return list of all .smk files plus Snakefile."""
    files = [SNAKEFILE]
    if os.path.isdir(RULES_DIR):
        for fname in sorted(os.listdir(RULES_DIR)):
            if fname.endswith(".smk"):
                files.append(os.path.join(RULES_DIR, fname))
    return files


def extract_shell_blocks(content):
    """Extract shell: blocks from Snakemake content.

    Returns list of (line_number, block_text) tuples.
    line_number is the physical line number of the 'shell:' keyword.

    Conservative parser: finds 'shell:' at the start of a line (after
    optional whitespace), then collects content until a non-indented
    keyword or end of file. Multiline blocks with triple-quotes,
    single-line blocks with double-quotes, and heredoc-style are all
    handled.
    """
    blocks = []
    lines = content.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]

        # Detect shell: keyword at start of line (possibly indented)
        if re.match(r"^\s+shell\s*:", line):
            shell_line_no = i + 1  # 1-indexed
            inline_content = line.split(":", 1)[1].strip()

            # Simple one-line shell blocks, e.g. shell: "cmd".
            # Multiline blocks may write shell: followed by quotes on the next
            # line, which falls through to the indented block parser below.
            if inline_content:
                blocks.append((shell_line_no, inline_content))
                i += 1
                continue

            # Collect the shell content
            i += 1
            shell_lines = []

            # Determine indentation of the shell keyword
            shell_indent = len(line) - len(line.lstrip())

            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()

                # Empty line — part of shell block
                if not next_stripped:
                    shell_lines.append(next_line)
                    i += 1
                    continue

                next_indent = len(next_line) - len(next_line.lstrip())

                # If next line is at same or lesser indent AND looks like
                # a Snakemake keyword (not shell content), stop
                if next_indent <= shell_indent:
                    # Check if it looks like a keyword line
                    if re.match(r"^\s*(rule|checkpoint|include|configfile|"
                                r"use|localrules|onsuccess|onerror|"
                                r"shell|script|wrapper|run|"
                                r"input|output|params|log|threads|"
                                r"resources|priority|wildcard_constraints|"
                                r"conda|container|message|benchmark|"
                                r"shadow|group|default_target)\b",
                                next_line):
                        break
                    # If same indent but not a keyword, could be continuation

                shell_lines.append(next_line)
                i += 1

            if shell_lines:
                block_text = "\n".join(shell_lines)
                blocks.append((shell_line_no, block_text))
        else:
            i += 1

    return blocks


def check_block(block_text):
    """Check a single shell block. Returns error string or None."""
    stripped = block_text.strip()

    if not stripped:
        return f"empty shell block"

    # Check for exemption marker anywhere in the block
    exempt_match = EXEMPTION_PATTERN.search(block_text)
    if exempt_match:
        reason = exempt_match.group(1).strip()
        if not reason:
            return (
                f"exemption marker found but no reason provided "
                f"(add reason after '# no pipefail:')"
            )
        return None  # Exempted with reason — OK

    # For single-line blocks (double-quoted)
    if stripped.startswith('"') and ';' in stripped:
        prefix = stripped.split(';')[0].strip().strip('"')
        if GUARD_REQUIRED in prefix:
            return None
    elif stripped.startswith('"') and ';' not in stripped:
        if GUARD_REQUIRED not in stripped:
            return f"missing '{GUARD_REQUIRED}' guard"
        return None

    # For multiline blocks (triple-quoted)
    lines = stripped.split("\n")
    first_content_line = None
    for line in lines:
        ls = line.strip()
        if ls and not ls.startswith('"""') and not ls.startswith("'''"):
            first_content_line = ls
            break

    if first_content_line is None:
        return f"could not parse shell block content"

    if GUARD_REQUIRED in first_content_line:
        return None

    return f"missing '{GUARD_REQUIRED}' guard"


def main():
    files = find_smk_files()
    errors = []
    total_blocks = 0

    for filepath in files:
        relpath = os.path.relpath(filepath, REPO_ROOT)
        with open(filepath) as fh:
            content = fh.read()

        blocks = extract_shell_blocks(content)
        directive_count = len(re.findall(r"(?m)^\s+shell\s*:", content))
        if len(blocks) != directive_count:
            errors.append(
                f"  {relpath}: parser extracted {len(blocks)} shell block(s), "
                f"but found {directive_count} shell directive(s)"
            )

        for line_no, block_text in blocks:
            total_blocks += 1
            error = check_block(block_text)
            if error:
                errors.append(f"  {relpath}: line {line_no}: {error}")

    if errors:
        print(f"FAIL: {len(errors)} shell block(s) missing fail-fast guard:\n")
        for err in errors:
            print(err)
        print(f"\n  Add '{GUARD_REQUIRED}' at the start of each shell block,")
        print(f"  or add '# no pipefail: <reason>' for intentional exemptions.")
        print(f"\n  Total shell blocks checked: {total_blocks}")
        return 1

    print(f"PASS: All {total_blocks} shell block(s) have '{GUARD_REQUIRED}'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
