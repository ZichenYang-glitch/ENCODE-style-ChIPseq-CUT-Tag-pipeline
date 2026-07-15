"""Require version constraints for core tools in each Conda environment."""

from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
ENVS_DIR = REPO_ROOT / "workflow" / "envs"

CORE_TOOLS = {
    "fastqc",
    "trim-galore",
    "cutadapt",
    "bowtie2",
    "samtools",
    "bedtools",
    "macs3",
    "deeptools",
    "idr",
    "multiqc",
    "picard",
    "phantompeakqualtools",
    "preseq",
    "seacr",
    "ucsc-bedgraphtobigwig",
}

VERSION_CONSTRAINT_CHARS = set("=<>~")


def _extract_package_name(dep_text):
    text = dep_text.strip()
    for sep in (" ", "\t"):
        if sep in text:
            return text.split(sep)[0].strip()
    return text


def _has_version_constraint(dep_text):
    pkg = _extract_package_name(dep_text)
    rest = dep_text[len(pkg) :].strip()
    if not rest:
        return False
    return any(c in VERSION_CONSTRAINT_CHARS for c in rest)


def _core_tool_entries():
    entries = []
    for env_path in sorted(ENVS_DIR.glob("*.yml")):
        with open(env_path, encoding="utf-8") as fh:
            spec = yaml.safe_load(fh)
        for dep in spec.get("dependencies", []):
            if isinstance(dep, str):
                pkg = _extract_package_name(dep)
                if pkg in CORE_TOOLS:
                    entries.append((env_path.name, pkg, dep))
            elif isinstance(dep, dict):
                # pip: [...] entries are not core CLI tools in this repo
                pass
    return entries


@pytest.mark.parametrize(
    "env_name,pkg_name,dep_text",
    _core_tool_entries(),
    ids=lambda t: f"{t[0]}::{t[1]}",
)
def test_core_tool_has_version_constraint(env_name, pkg_name, dep_text):
    assert _has_version_constraint(dep_text), (
        f"{env_name}: {pkg_name!r} has no version constraint ({dep_text!r})"
    )
