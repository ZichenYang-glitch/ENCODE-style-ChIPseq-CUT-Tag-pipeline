"""Require guarded Snakemake shell blocks or an explicit exemption."""

import re
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
RULES_DIR = REPO_ROOT / "workflow" / "rules"
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"

GUARD_REQUIRED = "set -e -o pipefail"
EXEMPTION_PATTERN = re.compile(r"#\s*no\s+pipefail\s*:\s*(.*)")


def _smk_files():
    files = [SNAKEFILE]
    if RULES_DIR.is_dir():
        files.extend(sorted(RULES_DIR.glob("*.smk")))
    return files


def _extract_shell_blocks(content):
    """Yield (line_number, block_text) for each shell: block.

    The parser is deliberately conservative and only inspects shell blocks.
    """
    blocks = []
    lines = content.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if re.match(r"^\s+shell\s*:", line):
            shell_line_no = i + 1
            inline_content = line.split(":", 1)[1].strip()

            if inline_content:
                blocks.append((shell_line_no, inline_content))
                i += 1
                continue

            i += 1
            shell_lines = []
            shell_indent = len(line) - len(line.lstrip())

            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()

                if not next_stripped:
                    shell_lines.append(next_line)
                    i += 1
                    continue

                next_indent = len(next_line) - len(next_line.lstrip())

                if next_indent <= shell_indent:
                    if re.match(
                        r"^\s*(rule|checkpoint|include|configfile|use|localrules|"
                        r"onsuccess|onerror|shell|script|wrapper|run|input|output|"
                        r"params|log|threads|resources|priority|wildcard_constraints|"
                        r"conda|container|message|benchmark|shadow|group|default_target)\b",
                        next_line,
                    ):
                        break

                shell_lines.append(next_line)
                i += 1

            if shell_lines:
                blocks.append((shell_line_no, "\n".join(shell_lines)))
        else:
            i += 1

    return blocks


def _check_block(block_text):
    stripped = block_text.strip()

    if not stripped:
        return "empty shell block"

    exempt_match = EXEMPTION_PATTERN.search(block_text)
    if exempt_match:
        reason = exempt_match.group(1).strip()
        if not reason:
            return "exemption marker found but no reason provided"
        return None

    if stripped.startswith('"') and ";" in stripped:
        prefix = stripped.split(";")[0].strip().strip('"')
        if GUARD_REQUIRED in prefix:
            return None
    elif stripped.startswith('"'):
        if GUARD_REQUIRED not in stripped:
            return f"missing '{GUARD_REQUIRED}' guard"
        return None

    lines = stripped.split("\n")
    first_content_line = None
    for line in lines:
        ls = line.strip()
        if ls and not ls.startswith('"""') and not ls.startswith("'''"):
            first_content_line = ls
            break

    if first_content_line is None:
        return "could not parse shell block content"

    if GUARD_REQUIRED in first_content_line:
        return None

    return f"missing '{GUARD_REQUIRED}' guard"


@pytest.mark.parametrize("smk_path", _smk_files(), ids=lambda p: p.name)
def test_all_shell_blocks_have_failfast_guard(smk_path):
    """Every shell block starts with `set -e -o pipefail` or is exempted."""
    content = smk_path.read_text(encoding="utf-8")
    blocks = _extract_shell_blocks(content)
    directive_count = len(re.findall(r"(?m)^\s+shell\s*:", content))
    assert len(blocks) == directive_count, (
        f"Parser extracted {len(blocks)} shell block(s), "
        f"but found {directive_count} shell directive(s) in {smk_path.name}"
    )

    errors = []
    for line_no, block_text in blocks:
        error = _check_block(block_text)
        if error:
            errors.append(f"{smk_path.name}:{line_no}: {error}")

    assert not errors, "\n".join(errors)
