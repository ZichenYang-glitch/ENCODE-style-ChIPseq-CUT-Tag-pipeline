"""Safe transition from the unpublished legacy candidate distribution."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import venv


REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGING_FILES = ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in")
CONSOLE_SCRIPTS = (
    "helixweave",
    "encode-validate",
    "encode-manifest",
    "encode-dag",
    "encode-worker",
)


def _run(
    *argv: str,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    environment = {
        key: value for key, value in os.environ.items() if key != "PYTHONPATH"
    }
    environment.update(
        {
            "PYTHONNOUSERSITE": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
        }
    )
    return subprocess.run(
        argv,
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def _copy_build_source(destination: Path, *, legacy: bool) -> None:
    destination.mkdir()
    for filename in PACKAGING_FILES:
        shutil.copy2(REPO_ROOT / filename, destination / filename)
    shutil.copytree(
        REPO_ROOT / "src",
        destination / "src",
        ignore=shutil.ignore_patterns("*.egg-info", "__pycache__", "*.pyc"),
    )
    if legacy:
        pyproject = destination / "pyproject.toml"
        content = pyproject.read_text(encoding="utf-8")
        content = content.replace('name = "helixweave"', 'name = "encode-pipeline"', 1)
        content = content.replace(
            'helixweave = "encode_pipeline.cli.app:main"\n',
            "",
            1,
        )
        pyproject.write_text(content, encoding="utf-8")


def _build_wheel(source: Path, wheelhouse: Path) -> Path:
    wheelhouse.mkdir()
    completed = _run(
        sys.executable,
        "-m",
        "pip",
        "wheel",
        "--no-index",
        "--no-deps",
        "--no-build-isolation",
        "--wheel-dir",
        str(wheelhouse),
        str(source),
    )
    assert completed.returncode == 0, completed.stderr
    wheels = tuple(wheelhouse.glob("*.whl"))
    assert len(wheels) == 1
    return wheels[0]


def _create_venv(root: Path) -> Path:
    venv.EnvBuilder(with_pip=False).create(root)
    return root / "bin" / "python"


def _pip(python: Path, *argv: str) -> subprocess.CompletedProcess[str]:
    return _run(
        sys.executable,
        "-m",
        "pip",
        "--disable-pip-version-check",
        "--python",
        str(python),
        *argv,
        cwd=python.parent.parent,
    )


def test_legacy_candidate_must_be_uninstalled_before_helixweave(
    tmp_path: Path,
) -> None:
    legacy_source = tmp_path / "legacy-source"
    current_source = tmp_path / "current-source"
    _copy_build_source(legacy_source, legacy=True)
    _copy_build_source(current_source, legacy=False)

    legacy_wheel = _build_wheel(legacy_source, tmp_path / "legacy-wheelhouse")
    current_wheel = _build_wheel(current_source, tmp_path / "current-wheelhouse")
    assert legacy_wheel.name == "encode_pipeline-0.3.0-py3-none-any.whl"
    assert current_wheel.name == "helixweave-0.3.0-py3-none-any.whl"

    unsafe_python = _create_venv(tmp_path / "unsafe")
    assert (
        _pip(
            unsafe_python,
            "install",
            "--no-index",
            "--no-deps",
            str(legacy_wheel),
            str(current_wheel),
        ).returncode
        == 0
    )
    unsafe_probe = _run(
        str(unsafe_python),
        "-c",
        (
            "import json; from importlib import metadata; "
            "print(json.dumps(metadata.packages_distributions()['encode_pipeline']))"
        ),
        cwd=tmp_path,
    )
    assert unsafe_probe.returncode == 0, unsafe_probe.stderr
    assert set(json.loads(unsafe_probe.stdout)) == {"encode-pipeline", "helixweave"}

    safe_python = _create_venv(tmp_path / "safe")
    assert (
        _pip(
            safe_python,
            "install",
            "--no-index",
            "--no-deps",
            str(legacy_wheel),
        ).returncode
        == 0
    )
    uninstall = _pip(
        safe_python,
        "uninstall",
        "-y",
        "encode-pipeline",
    )
    assert uninstall.returncode == 0, uninstall.stderr
    removed_probe = _run(
        str(safe_python),
        "-c",
        """
from importlib import metadata, util

try:
    metadata.version("encode-pipeline")
except metadata.PackageNotFoundError:
    pass
else:
    raise AssertionError("legacy distribution metadata remains installed")
assert metadata.packages_distributions().get("encode_pipeline", []) == []
assert util.find_spec("encode_pipeline") is None
""",
        cwd=tmp_path,
    )
    assert removed_probe.returncode == 0, removed_probe.stderr
    scripts = safe_python.parent
    for name in CONSOLE_SCRIPTS:
        assert not (scripts / name).exists()

    assert (
        _pip(
            safe_python,
            "install",
            "--no-index",
            "--no-deps",
            str(current_wheel),
        ).returncode
        == 0
    )
    safe_probe = _run(
        str(safe_python),
        "-c",
        """
import json
from importlib import metadata, resources

assert metadata.version("helixweave") == "0.3.0"
try:
    metadata.version("encode-pipeline")
except metadata.PackageNotFoundError:
    pass
else:
    raise AssertionError("legacy distribution metadata remains installed")
assert metadata.packages_distributions()["encode_pipeline"] == ["helixweave"]
assert resources.files("encode_pipeline.artifacts").joinpath(
    "artifact-inventory.yaml"
).is_file()
print(json.dumps(metadata.packages_distributions()["encode_pipeline"]))
""",
        cwd=tmp_path,
    )
    assert safe_probe.returncode == 0, safe_probe.stderr
    assert json.loads(safe_probe.stdout) == ["helixweave"]

    for name in CONSOLE_SCRIPTS:
        assert (scripts / name).is_file()
    module_help = _run(
        str(safe_python),
        "-m",
        "encode_pipeline",
        "--help",
        cwd=tmp_path,
    )
    console_help = _run(
        str(scripts / "helixweave"),
        "--help",
        cwd=tmp_path,
    )
    assert module_help.returncode == 0, module_help.stderr
    assert console_help.returncode == 0, console_help.stderr
    assert module_help.stdout == console_help.stdout
