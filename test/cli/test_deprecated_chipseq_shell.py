"""User-visible contract for the deprecated shell entry point."""

import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "chipseq.sh"


def test_deprecated_shell_has_valid_bash_syntax():
    result = subprocess.run(
        ["bash", "-n", str(SCRIPT)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_deprecated_shell_fails_closed_and_redirects_to_snakemake():
    result = subprocess.run(
        ["bash", str(SCRIPT)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode != 0
    assert result.stdout == ""
    assert "deprecated" in result.stderr.lower()
    assert "snakemake" in result.stderr.lower()
    assert "workflow/Snakefile" in result.stderr
