"""Fail-closed contracts for the checkout-only Python bootstrap.

These tests deliberately invoke a fresh interpreter with ``-I -S``.  The
fixtures model the editable metadata emitted by pip/setuptools without relying
on an editable installation already present on the test machine.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from importlib import util as importlib_util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace
from urllib.parse import quote
import venv

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
BOOTSTRAP_SCRIPT = REPOSITORY_ROOT / "scripts" / "checkout_bootstrap.py"
PROVENANCE_SCRIPT = REPOSITORY_ROOT / "scripts" / "source_provenance.py"
_PUBLIC_REASON = re.compile(r"\[[a-z][a-z0-9_]*\]")
_BOOTSTRAP_SPEC = importlib_util.spec_from_file_location(
    "_helixweave_checkout_bootstrap_test",
    BOOTSTRAP_SCRIPT,
)
assert _BOOTSTRAP_SPEC is not None
assert _BOOTSTRAP_SPEC.loader is not None
BOOTSTRAP = importlib_util.module_from_spec(_BOOTSTRAP_SPEC)
sys.modules[_BOOTSTRAP_SPEC.name] = BOOTSTRAP
_BOOTSTRAP_SPEC.loader.exec_module(BOOTSTRAP)


@dataclass(frozen=True)
class _CheckoutCase:
    repository: Path
    python: Path
    site_packages: Path
    product_marker: Path
    startup_marker: Path
    target_output: Path


def _file_url(path: Path) -> str:
    return f"file://{quote(str(path.resolve()))}"


def _create_venv(root: Path) -> tuple[Path, Path]:
    venv.EnvBuilder(with_pip=False).create(root)
    python = root / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    completed = subprocess.run(
        [
            str(python),
            "-c",
            "import sysconfig; print(sysconfig.get_path('purelib'))",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return python, Path(completed.stdout.strip())


def _write_editable_distribution(
    site_packages: Path,
    repository: Path,
) -> Path:
    """Write the ownership and editable identity emitted for HelixWeave."""
    site_packages.mkdir(parents=True, exist_ok=True)
    pth_name = "__editable__.helixweave-0.3.0.pth"
    pth = site_packages / pth_name
    pth.write_text(f"{repository.resolve() / 'src'}\n", encoding="utf-8")

    metadata = site_packages / "helixweave-0.3.0.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        "Metadata-Version: 2.4\nName: helixweave\nVersion: 0.3.0\n",
        encoding="utf-8",
    )
    (metadata / "top_level.txt").write_text(
        "encode_pipeline\n",
        encoding="utf-8",
    )
    (metadata / "direct_url.json").write_text(
        json.dumps(
            {
                "dir_info": {"editable": True},
                "url": _file_url(repository),
            }
        ),
        encoding="utf-8",
    )
    (metadata / "RECORD").write_text(
        f"{pth_name},,\n"
        "helixweave-0.3.0.dist-info/METADATA,,\n"
        "helixweave-0.3.0.dist-info/RECORD,,\n"
        "helixweave-0.3.0.dist-info/direct_url.json,,\n"
        "helixweave-0.3.0.dist-info/top_level.txt,,\n",
        encoding="utf-8",
    )
    return pth


def _copy_bootstrap_sources(repository: Path) -> None:
    scripts = repository / "scripts"
    scripts.mkdir(parents=True)
    if BOOTSTRAP_SCRIPT.is_file():
        shutil.copy2(BOOTSTRAP_SCRIPT, scripts / BOOTSTRAP_SCRIPT.name)
    shutil.copy2(PROVENANCE_SCRIPT, scripts / PROVENANCE_SCRIPT.name)


def _write_checkout(repository: Path, product_marker: Path) -> None:
    repository.mkdir(parents=True)
    (repository / "pyproject.toml").write_text(
        '[project]\nname = "helixweave"\nversion = "0.3.0"\n',
        encoding="utf-8",
    )
    package = repository / "src" / "encode_pipeline"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        "import os\n"
        "from pathlib import Path\n"
        "marker = os.environ.get('PRODUCT_MARKER')\n"
        "if marker:\n"
        "    Path(marker).write_text('product imported', encoding='utf-8')\n"
        '__version__ = "0.3.0"\n',
        encoding="utf-8",
    )
    (repository / "config").mkdir()
    (repository / "config" / "config.yaml").write_text("{}\n", encoding="utf-8")
    (repository / "test").mkdir()
    (repository / "test" / "test_probe.py").write_text(
        "import encode_pipeline\n\n"
        "def test_current_checkout_imported():\n"
        "    assert encode_pipeline.__version__ == '0.3.0'\n",
        encoding="utf-8",
    )
    _copy_bootstrap_sources(repository)

    target_body = (
        "import os\n"
        "from pathlib import Path\n"
        "Path(os.environ['TARGET_OUTPUT']).write_text("
        "'target started', encoding='utf-8')\n"
        "import encode_pipeline\n"
    )
    scripts = repository / "scripts"
    for name in (
        "export_openapi.py",
        "validate_samples.py",
        "run_local_platform.py",
    ):
        (scripts / name).write_text(target_body, encoding="utf-8")

    assert not product_marker.exists()


def _create_case(tmp_path: Path, *, repository_name: str = "current") -> _CheckoutCase:
    repository = tmp_path / repository_name
    product_marker = tmp_path / "product-imported"
    startup_marker = tmp_path / "startup-imported"
    target_output = tmp_path / "target-output"
    _write_checkout(repository, product_marker)
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_editable_distribution(site_packages, repository)
    return _CheckoutCase(
        repository=repository,
        python=python,
        site_packages=site_packages,
        product_marker=product_marker,
        startup_marker=startup_marker,
        target_output=target_output,
    )


def _clean_environment(case: _CheckoutCase) -> dict[str, str]:
    environment = {
        key: value
        for key, value in os.environ.items()
        if key
        not in {
            "PYTHONHOME",
            "PYTHONPATH",
            "PYTEST_ADDOPTS",
            "PYTEST_PLUGINS",
            "PYTEST_DISABLE_PLUGIN_AUTOLOAD",
        }
    }
    environment.update(
        {
            "PRODUCT_MARKER": str(case.product_marker),
            "STARTUP_MARKER": str(case.startup_marker),
            "TARGET_OUTPUT": str(case.target_output),
        }
    )
    return environment


def _run_bootstrap(
    case: _CheckoutCase,
    subcommand: str,
    *arguments: str,
    repository_argument: Path | None = None,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    child_environment = _clean_environment(case)
    if environment:
        child_environment.update(environment)
    script = case.repository / "scripts" / BOOTSTRAP_SCRIPT.name
    return subprocess.run(
        [
            str(case.python),
            "-I",
            "-S",
            str(script),
            "--repository-root",
            str(repository_argument or case.repository),
            subcommand,
            *arguments,
        ],
        cwd=case.repository,
        env=child_environment,
        capture_output=True,
        text=True,
        check=False,
    )


def _assert_no_effects(case: _CheckoutCase) -> None:
    assert not case.startup_marker.exists()
    assert not case.product_marker.exists()
    assert not case.target_output.exists()
    assert not tuple(case.repository.rglob("__pycache__"))
    assert not tuple(case.repository.rglob(".pytest_cache"))


def _assert_safe_failure(
    result: subprocess.CompletedProcess[str],
    case: _CheckoutCase,
    *private_paths: Path,
) -> None:
    assert result.returncode != 0
    assert _PUBLIC_REASON.search(result.stderr), result.stderr
    assert "Traceback" not in result.stderr
    for path in (case.repository, case.site_packages, *private_paths):
        assert str(path) not in result.stderr
        assert str(path.resolve()) not in result.stderr
    _assert_no_effects(case)


def _write_stale_checkout(root: Path, marker: Path) -> Path:
    package = root / "src" / "encode_pipeline"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        "import os\n"
        "from pathlib import Path\n"
        f"Path({str(marker)!r}).write_text('stale imported', encoding='utf-8')\n",
        encoding="utf-8",
    )
    return root


def _write_sitecustomize_attack(case: _CheckoutCase, stale: Path) -> None:
    (case.site_packages / "sitecustomize.py").write_text(
        "import os\n"
        "from pathlib import Path\n"
        "import sys\n"
        f"Path({str(case.startup_marker)!r}).write_text("
        "'sitecustomize ran', encoding='utf-8')\n"
        f"sys.path.insert(0, {str(stale / 'src')!r})\n"
        "import encode_pipeline\n",
        encoding="utf-8",
    )


def _write_executable_pth_attack(case: _CheckoutCase, stale: Path) -> Path:
    attack = case.site_packages / "orphan-stale.pth"
    attack.write_text(
        "import pathlib,sys;"
        f"pathlib.Path({str(case.startup_marker)!r}).write_text("
        "'pth ran',encoding='utf-8');"
        f"sys.path.insert(0,{str(stale / 'src')!r});"
        '__import__("encode_"+"pipeline")\n',
        encoding="utf-8",
    )
    return attack


@pytest.mark.parametrize(
    ("subcommand", "arguments", "attack_kind"),
    (
        ("openapi", ("--output", "openapi.json"), "sitecustomize"),
        ("validate", ("--config", "config/config.yaml"), "executable-pth"),
        ("local-platform", ("doctor",), "executable-pth"),
    ),
)
def test_product_commands_reject_startup_hooks_before_any_import_or_effect(
    tmp_path: Path,
    subcommand: str,
    arguments: tuple[str, ...],
    attack_kind: str,
) -> None:
    case = _create_case(tmp_path)
    stale = _write_stale_checkout(tmp_path / "private-stale", case.product_marker)
    if attack_kind == "sitecustomize":
        _write_sitecustomize_attack(case, stale)
    else:
        _write_executable_pth_attack(case, stale)

    result = _run_bootstrap(case, subcommand, *arguments)

    _assert_safe_failure(result, case, stale)
    assert "startup" in result.stderr.lower() or ".pth" in result.stderr.lower()


def test_pytest_rejects_automatic_plugin_injection_before_plugin_import(
    tmp_path: Path,
) -> None:
    case = _create_case(tmp_path)
    stale = _write_stale_checkout(tmp_path / "private-stale", case.product_marker)
    plugin = case.site_packages / "stale_pytest_plugin.py"
    plugin.write_text(
        "import pathlib,sys\n"
        f"pathlib.Path({str(case.startup_marker)!r}).write_text("
        "'plugin ran', encoding='utf-8')\n"
        f"sys.path.insert(0, {str(stale / 'src')!r})\n"
        "import encode_pipeline\n",
        encoding="utf-8",
    )

    result = _run_bootstrap(
        case,
        "pytest",
        "--collect-only",
        "test/test_probe.py",
        environment={"PYTEST_PLUGINS": "stale_pytest_plugin"},
    )

    _assert_safe_failure(result, case, stale, plugin)
    assert "plugin" in result.stderr.lower()


def test_verify_checkout_accepts_a_valid_current_editable(
    tmp_path: Path,
) -> None:
    case = _create_case(tmp_path)

    result = _run_bootstrap(case, "verify-checkout")

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
    _assert_no_effects(case)


def test_verify_checkout_accepts_a_checkout_path_containing_spaces(
    tmp_path: Path,
) -> None:
    case = _create_case(tmp_path, repository_name="current checkout with spaces")

    result = _run_bootstrap(case, "verify-checkout")

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
    _assert_no_effects(case)


def test_verify_checkout_canonicalizes_a_repository_symlink_alias(
    tmp_path: Path,
) -> None:
    case = _create_case(tmp_path)
    alias = tmp_path / "checkout-alias"
    try:
        alias.symlink_to(case.repository, target_is_directory=True)
    except OSError as error:
        pytest.fail(f"test platform must support directory symlinks: {error}")

    result = _run_bootstrap(
        case,
        "verify-checkout",
        repository_argument=alias,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
    _assert_no_effects(case)


def test_startup_failure_does_not_disclose_private_paths(
    tmp_path: Path,
) -> None:
    case = _create_case(tmp_path, repository_name="private current checkout")
    stale = _write_stale_checkout(
        tmp_path / "private stale checkout",
        case.product_marker,
    )
    attack = _write_executable_pth_attack(case, stale)

    result = _run_bootstrap(case, "openapi", "--output", "openapi.json")

    _assert_safe_failure(result, case, stale, attack)


def _content_snapshot(*roots: Path) -> dict[str, str]:
    snapshot: dict[str, str] = {}
    for index, root in enumerate(roots):
        for path in sorted(root.rglob("*")):
            if path.is_file():
                key = f"{index}:{path.relative_to(root)}"
                snapshot[key] = hashlib.sha256(path.read_bytes()).hexdigest()
            elif path.is_symlink():
                key = f"{index}:{path.relative_to(root)}"
                snapshot[key] = f"symlink:{os.readlink(path)}"
    return snapshot


def test_verify_checkout_is_read_only_and_repeatable(tmp_path: Path) -> None:
    case = _create_case(tmp_path)
    before = _content_snapshot(case.repository, case.site_packages)

    first = _run_bootstrap(case, "verify-checkout")
    after_first = _content_snapshot(case.repository, case.site_packages)
    second = _run_bootstrap(case, "verify-checkout")
    after_second = _content_snapshot(case.repository, case.site_packages)

    assert first.returncode == 0, first.stderr
    assert second.returncode == 0, second.stderr
    assert first.stdout == second.stdout == ""
    assert first.stderr == second.stderr == ""
    assert before == after_first == after_second
    _assert_no_effects(case)


@pytest.mark.parametrize(
    ("command", "expected_result"),
    (
        ("verify-checkout", 0),
        ("installed-artifact", 0),
        ("pytest", 31),
        ("openapi", 37),
        ("validate", 37),
        ("local-platform", 37),
    ),
)
def test_bootstrap_main_dispatches_only_fixed_audited_commands(
    monkeypatch: pytest.MonkeyPatch,
    command: str,
    expected_result: int,
) -> None:
    calls: list[tuple[str, object]] = []

    class FakeSourceProvenanceError(Exception):
        pass

    class FakeProvenance:
        SourceProvenanceError = FakeSourceProvenanceError

        @staticmethod
        def bootstrap_checkout(repository_root: Path) -> None:
            calls.append(("checkout", repository_root))

        @staticmethod
        def bootstrap_installed_artifact() -> None:
            calls.append(("installed", None))

    fake_sys = SimpleNamespace(
        dont_write_bytecode=False,
        flags=SimpleNamespace(isolated=1, no_site=1),
        stderr=sys.stderr,
    )
    monkeypatch.setattr(BOOTSTRAP, "sys", fake_sys)
    monkeypatch.setattr(
        BOOTSTRAP,
        "_parse_arguments",
        lambda: (
            argparse.Namespace(repository_root=REPOSITORY_ROOT, command=command),
            ["payload"],
        ),
    )
    monkeypatch.setattr(BOOTSTRAP, "_load_provenance", lambda path: FakeProvenance)
    monkeypatch.setattr(
        BOOTSTRAP,
        "_reject_pytest_plugin_injection",
        lambda arguments: calls.append(("plugins", tuple(arguments))),
    )
    monkeypatch.setattr(
        BOOTSTRAP,
        "_run_pytest",
        lambda arguments: calls.append(("pytest", tuple(arguments))) or 31,
    )
    monkeypatch.setattr(
        BOOTSTRAP,
        "_run_script",
        lambda path, arguments: calls.append((path.name, tuple(arguments))) or 37,
    )

    assert BOOTSTRAP.main() == expected_result
    assert fake_sys.dont_write_bytecode is True
    if command == "installed-artifact":
        assert ("installed", None) in calls
    else:
        assert ("checkout", REPOSITORY_ROOT) in calls
    if command == "pytest":
        assert ("plugins", ("payload",)) in calls
        assert ("pytest", ("payload",)) in calls
    elif command in {"openapi", "validate", "local-platform"}:
        assert any(call[0].endswith(".py") for call in calls)


def test_bootstrap_main_redacts_provenance_failure(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    class FakeSourceProvenanceError(Exception):
        def __str__(self) -> str:
            return "source provenance check failed [test_failure]: repair environment"

    class FakeProvenance:
        SourceProvenanceError = FakeSourceProvenanceError

        @staticmethod
        def bootstrap_checkout(repository_root: Path) -> None:
            raise FakeSourceProvenanceError

    monkeypatch.setattr(
        BOOTSTRAP,
        "sys",
        SimpleNamespace(
            dont_write_bytecode=False,
            flags=SimpleNamespace(isolated=1, no_site=1),
            stderr=sys.stderr,
        ),
    )
    monkeypatch.setattr(
        BOOTSTRAP,
        "_parse_arguments",
        lambda: (
            argparse.Namespace(
                repository_root=REPOSITORY_ROOT,
                command="verify-checkout",
            ),
            [],
        ),
    )
    monkeypatch.setattr(BOOTSTRAP, "_load_provenance", lambda path: FakeProvenance)

    assert BOOTSTRAP.main() == 2
    assert capsys.readouterr().err == (
        "source provenance check failed [test_failure]: repair environment\n"
    )


def test_bootstrap_requires_isolated_no_site_startup(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        BOOTSTRAP,
        "sys",
        SimpleNamespace(
            dont_write_bytecode=False,
            flags=SimpleNamespace(isolated=0, no_site=1),
            stderr=sys.stderr,
        ),
    )

    with pytest.raises(SystemExit, match="2"):
        BOOTSTRAP.main()

    assert "[bootstrap_startup_unsafe]" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("environment_name", "arguments"),
    (
        ("PYTEST_ADDOPTS", []),
        ("PYTEST_PLUGINS", []),
        (None, ["-p", "unsafe"]),
        (None, ["--plugins=unsafe"]),
    ),
)
def test_bootstrap_rejects_every_pytest_plugin_injection_surface(
    monkeypatch: pytest.MonkeyPatch,
    environment_name: str | None,
    arguments: list[str],
) -> None:
    monkeypatch.delenv("PYTEST_ADDOPTS", raising=False)
    monkeypatch.delenv("PYTEST_PLUGINS", raising=False)
    if environment_name is not None:
        monkeypatch.setenv(environment_name, "unsafe")

    with pytest.raises(SystemExit, match="2"):
        BOOTSTRAP._reject_pytest_plugin_injection(arguments)


@pytest.mark.parametrize(
    ("exit_code", "expected"),
    (
        (None, 0),
        (17, 17),
        ("failure", 1),
    ),
)
def test_bootstrap_script_runner_normalizes_exit_status(
    tmp_path: Path,
    exit_code: object,
    expected: int,
) -> None:
    script = tmp_path / "target.py"
    script.write_text(f"raise SystemExit({exit_code!r})\n", encoding="utf-8")

    assert BOOTSTRAP._run_script(script, ["argument"]) == expected


def test_bootstrap_loads_the_standard_library_guard() -> None:
    product_module = sys.modules.get("encode_pipeline")
    module = BOOTSTRAP._load_provenance(PROVENANCE_SCRIPT)

    assert module.DISTRIBUTION_NAME == "helixweave"
    assert sys.modules.get("encode_pipeline") is product_module


def test_bootstrap_parser_preserves_payload_arguments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(BOOTSTRAP_SCRIPT),
            "--repository-root",
            str(REPOSITORY_ROOT),
            "openapi",
            "--output",
            "schema.json",
        ],
    )

    arguments, payload = BOOTSTRAP._parse_arguments()

    assert arguments.repository_root == REPOSITORY_ROOT
    assert arguments.command == "openapi"
    assert payload == ["--output", "schema.json"]


def test_bootstrap_rejects_an_unloadable_guard(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        BOOTSTRAP.importlib.util,
        "spec_from_file_location",
        lambda *args: None,
    )

    with pytest.raises(SystemExit, match="2"):
        BOOTSTRAP._load_provenance(PROVENANCE_SCRIPT)


def test_bootstrap_rejects_a_malformed_guard(tmp_path: Path) -> None:
    malformed = tmp_path / "source_provenance.py"
    malformed.write_text("this is not valid Python !\n", encoding="utf-8")

    with pytest.raises(SystemExit, match="2"):
        BOOTSTRAP._load_provenance(malformed)


def test_bootstrap_script_runner_accepts_an_implicit_success(tmp_path: Path) -> None:
    script = tmp_path / "target.py"
    script.write_text("value = 1\n", encoding="utf-8")

    assert BOOTSTRAP._run_script(script, []) == 0


@pytest.mark.parametrize("repository_kind", ("missing", "different"))
def test_bootstrap_rejects_an_invalid_repository_argument(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    repository_kind: str,
) -> None:
    repository = tmp_path / "repository"
    if repository_kind == "different":
        repository.mkdir()
    monkeypatch.setattr(
        BOOTSTRAP,
        "sys",
        SimpleNamespace(
            dont_write_bytecode=False,
            flags=SimpleNamespace(isolated=1, no_site=1),
            stderr=sys.stderr,
        ),
    )
    monkeypatch.setattr(
        BOOTSTRAP,
        "_parse_arguments",
        lambda: (
            argparse.Namespace(
                repository_root=repository,
                command="verify-checkout",
            ),
            [],
        ),
    )

    with pytest.raises(SystemExit, match="2"):
        BOOTSTRAP.main()


def test_bootstrap_pytest_runner_loads_only_the_explicit_plugin(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, object] = {}

    def fake_main(arguments: list[str], *, plugins: list[object]) -> int:
        observed["arguments"] = arguments
        observed["plugins"] = plugins
        return 23

    monkeypatch.setattr(pytest, "main", fake_main)

    assert BOOTSTRAP._run_pytest(["test/probe.py"]) == 23
    assert observed["arguments"] == [
        "-p",
        "no:cacheprovider",
        "test/probe.py",
    ]
    assert len(observed["plugins"]) == 1
    assert os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] == "1"
