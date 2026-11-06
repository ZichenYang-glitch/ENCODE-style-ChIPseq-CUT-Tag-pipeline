"""The Snakemake process must load this checkout's encode_pipeline.

verify-checkout covers the interpreter that runs the bootstrap. These tests
execute Snakemake itself and read that process's verbose import trace.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap

from conftest import prepare_profile_workdir


REPO = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO / "workflow" / "Snakefile"
BOOTSTRAP = REPO / "scripts" / "checkout_bootstrap.py"
CHECKOUT_INIT = REPO / "src" / "encode_pipeline" / "__init__.py"
PROFILE = ("--profile", "none", "--workflow-profile", "none")
PROVENANCE_LINE = (
    "source provenance check failed [module_origin_mismatch]: "
    "run this workflow under the interpreter that loads this checkout"
)
INIT_FAILURE_LINE = (
    "source provenance check failed [module_initialization_failed]: "
    "this Snakemake process located encode_pipeline but could not finish loading it"
)


def _env(**overrides: str) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    env.pop("SNAKEMAKE_PROFILE", None)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.update(overrides)
    return env


def _run(
    argv: list[str],
    *,
    cwd: Path,
    env: dict[str, str],
    timeout: int = 180,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        argv,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )


def _loaded_inits(stderr: str) -> list[str]:
    """Return encode_pipeline initializer paths recorded by python -v."""

    suffix = "/encode_pipeline/__init__.py"
    found: list[str] = []
    for line in stderr.splitlines():
        origin = ""
        if " matches " in line:
            origin = line.split(" matches ", 1)[1].strip()
        elif line.startswith("# code object from "):
            origin = line.removeprefix("# code object from ").strip().strip("'\"")
        if origin.endswith(suffix):
            found.append(origin)
    return found


def _assert_no_scientific_output(workdir: Path) -> None:
    assert not (workdir / "results").exists()
    assert list(workdir.glob("**/*.bam")) == []
    assert list(workdir.glob("**/*.bai")) == []
    assert list(workdir.glob("**/*identity*")) == []
    assert list(workdir.glob("**/*effective*")) == []


def _profile_workdir() -> tuple[Path, Path]:
    workdir, config = prepare_profile_workdir(
        str(REPO / "test" / "profiles" / "chipseq_se_noctrl")
    )
    return Path(workdir), Path(config)


def _foreign_tree(root: Path) -> Path:
    package = root / "src" / "encode_pipeline"
    package.mkdir(parents=True)
    init = package / "__init__.py"
    init.write_text("FOREIGN = True\n", encoding="utf-8")
    return init


def test_matching_checkout_dry_run_uses_the_validator_source(tmp_path: Path) -> None:
    workdir, config = _profile_workdir()
    try:
        env = _env()
        env["XDG_CACHE_HOME"] = str(workdir / ".cache")
        verified = _run(
            [
                sys.executable,
                "-v",
                "-I",
                "-S",
                str(BOOTSTRAP),
                "--repository-root",
                str(REPO),
                "verify-checkout",
            ],
            cwd=REPO,
            env=env,
            timeout=60,
        )
        assert verified.returncode == 0, verified.stderr
        validated = _run(
            [
                sys.executable,
                "-v",
                "-I",
                "-S",
                str(BOOTSTRAP),
                "--repository-root",
                str(REPO),
                "validate",
                "--config",
                str(config),
            ],
            cwd=REPO,
            env=env,
            timeout=60,
        )
        assert validated.returncode == 0, validated.stderr
        assert "OK:" in validated.stdout
        assert _loaded_inits(validated.stderr)[:1] == [str(CHECKOUT_INIT)]

        dry_run = _run(
            [
                sys.executable,
                "-v",
                "-m",
                "snakemake",
                "-s",
                str(SNAKEFILE),
                "--configfile",
                str(config),
                "--directory",
                str(workdir),
                "--cores",
                "1",
                "--dry-run",
                *PROFILE,
            ],
            cwd=workdir,
            env=env,
        )
        output = dry_run.stdout + dry_run.stderr
        assert dry_run.returncode == 0, output
        assert "Building DAG of jobs" in output
        assert PROVENANCE_LINE not in output
        assert _loaded_inits(dry_run.stderr)[:1] == [str(CHECKOUT_INIT)]
        _assert_no_scientific_output(workdir)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def test_pythonpath_foreign_package_fails_before_scientific_output(
    tmp_path: Path,
) -> None:
    foreign_init = _foreign_tree(tmp_path / "old-checkout")
    workdir, config = _profile_workdir()
    try:
        env = _env(PYTHONPATH=str(foreign_init.parent.parent))
        env["XDG_CACHE_HOME"] = str(workdir / ".cache")
        # The same interpreter still passes the isolated bootstrap check.
        verified = _run(
            [
                sys.executable,
                "-I",
                "-S",
                str(BOOTSTRAP),
                "--repository-root",
                str(REPO),
                "verify-checkout",
            ],
            cwd=REPO,
            env=env,
            timeout=60,
        )
        assert verified.returncode == 0, verified.stderr

        result = _run(
            [
                sys.executable,
                "-v",
                "-m",
                "snakemake",
                "-s",
                str(SNAKEFILE),
                "--configfile",
                str(config),
                "--directory",
                str(workdir),
                "--cores",
                "1",
                "--dry-run",
                *PROFILE,
            ],
            cwd=workdir,
            env=env,
        )
        output = result.stdout + result.stderr
        assert result.returncode == 2, output
        assert PROVENANCE_LINE in result.stderr
        assert str(foreign_init) not in PROVENANCE_LINE
        assert "Building DAG of jobs" not in output
        assert _loaded_inits(result.stderr)[:1] == [str(foreign_init)]
        _assert_no_scientific_output(workdir)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def test_shebang_editable_mapping_fails_before_scientific_output() -> None:
    root = Path("/tmp") / f"pv-{os.getpid()}-shebang"
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir()
    workdir, config = _profile_workdir()
    try:
        foreign_init = _foreign_tree(root / "old-checkout")
        venv = root / "venv"
        subprocess.check_call(
            [sys.executable, "-m", "venv", "--system-site-packages", str(venv)],
            env=_env(),
        )
        site = next((venv / "lib").glob("python*/site-packages"))
        (site / "foreign.pth").write_text(
            str(foreign_init.parent.parent) + "\n",
            encoding="utf-8",
        )
        launcher = root / "snakemake"
        python = venv / "bin" / "python"
        launcher.write_text(
            textwrap.dedent(
                f"""\
                #!{python} -v
                import runpy
                raise SystemExit(runpy.run_module("snakemake", run_name="__main__"))
                """
            ),
            encoding="utf-8",
        )
        launcher.chmod(0o755)
        assert len(f"#!{python} -v") < 128

        env = _env()
        env["XDG_CACHE_HOME"] = str(workdir / ".cache")
        result = _run(
            [
                str(launcher),
                "-s",
                str(SNAKEFILE),
                "--configfile",
                str(config),
                "--directory",
                str(workdir),
                "--cores",
                "1",
                "--dry-run",
                *PROFILE,
            ],
            cwd=workdir,
            env=env,
        )
        output = result.stdout + result.stderr
        assert result.returncode == 2, output
        assert PROVENANCE_LINE in result.stderr
        assert "Building DAG of jobs" not in output
        assert _loaded_inits(result.stderr)[:1] == [str(foreign_init)]
        _assert_no_scientific_output(workdir)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
        shutil.rmtree(root, ignore_errors=True)


def _reject_broken_located_package(
    tmp_path: Path,
    init_source: str,
    *,
    hidden_text: str | None = None,
) -> None:
    foreign_root = tmp_path / "old-checkout" / "src"
    package = foreign_root / "encode_pipeline"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(init_source, encoding="utf-8")
    workdir = tmp_path / "work"
    workdir.mkdir()
    config = workdir / "config.yaml"
    config.write_text("threads: 1\n", encoding="utf-8")
    env = _env(PYTHONPATH=str(foreign_root))
    env["XDG_CACHE_HOME"] = str(workdir / ".cache")
    result = _run(
        [
            sys.executable,
            "-v",
            "-m",
            "snakemake",
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(config),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--dry-run",
            *PROFILE,
        ],
        cwd=workdir,
        env=env,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 2, output
    assert INIT_FAILURE_LINE in result.stderr
    assert PROVENANCE_LINE not in output
    if hidden_text is not None:
        assert hidden_text not in output
    assert str(CHECKOUT_INIT) not in _loaded_inits(result.stderr)
    assert "Building DAG of jobs" not in output
    _assert_no_scientific_output(workdir)


def test_foreign_missing_dependency_is_not_replaced_by_this_checkout(
    tmp_path: Path,
) -> None:
    _reject_broken_located_package(
        tmp_path,
        "import encode_pipeline_missing_dependency\n",
        hidden_text="No module named 'encode_pipeline_missing_dependency'",
    )


def test_foreign_explicit_importerror_is_not_replaced_by_this_checkout(
    tmp_path: Path,
) -> None:
    _reject_broken_located_package(
        tmp_path,
        "raise ImportError('forced initialization failure')\n",
        hidden_text="forced initialization failure",
    )


def test_missing_package_falls_back_to_this_checkout(tmp_path: Path) -> None:
    blocker = tmp_path / "sitecustomize.py"
    blocker.write_text(
        textwrap.dedent(
            """\
            import sys
            from pathlib import Path

            removed = []
            for entry in list(sys.path):
                if (Path(entry) / "encode_pipeline" / "__init__.py").is_file():
                    removed.append(entry)
            for entry in removed:
                while entry in sys.path:
                    sys.path.remove(entry)
            """
        ),
        encoding="utf-8",
    )
    workdir = tmp_path / "work"
    workdir.mkdir()
    config = workdir / "config.yaml"
    config.write_text("threads: -1\n", encoding="utf-8")
    env = _env(PYTHONPATH=str(tmp_path))
    env["XDG_CACHE_HOME"] = str(workdir / ".cache")
    result = _run(
        [
            sys.executable,
            "-v",
            "-m",
            "snakemake",
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(config),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--dry-run",
            *PROFILE,
        ],
        cwd=workdir,
        env=env,
    )
    output = result.stdout + result.stderr
    assert "source provenance check failed" not in output
    assert result.returncode != 0
    assert _loaded_inits(result.stderr)[:1] == [str(CHECKOUT_INIT)]
    assert "Building DAG of jobs" not in output
    _assert_no_scientific_output(workdir)
