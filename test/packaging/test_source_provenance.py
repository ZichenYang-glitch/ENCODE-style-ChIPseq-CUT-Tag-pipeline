"""Current-checkout and installed-artifact source provenance contracts."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from importlib import util as importlib_util
from importlib.util import find_spec
from types import SimpleNamespace
from urllib.parse import quote
import venv

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
GUARD_SCRIPT = REPOSITORY_ROOT / "scripts" / "source_provenance.py"
VALIDATOR_SCRIPT = REPOSITORY_ROOT / "scripts" / "validate_samples.py"
OPENAPI_SCRIPT = REPOSITORY_ROOT / "scripts" / "export_openapi.py"
LOCAL_PLATFORM_SCRIPT = REPOSITORY_ROOT / "scripts" / "run_local_platform.py"
_COVERAGE_SPEC = find_spec("coverage")
assert _COVERAGE_SPEC is not None
assert _COVERAGE_SPEC.origin is not None
COVERAGE_SUPPORT_ROOT = Path(_COVERAGE_SPEC.origin).resolve().parents[1]
_PROVENANCE_SPEC = importlib_util.spec_from_file_location(
    "_helixweave_source_provenance_test",
    GUARD_SCRIPT,
)
assert _PROVENANCE_SPEC is not None
assert _PROVENANCE_SPEC.loader is not None
PROVENANCE = importlib_util.module_from_spec(_PROVENANCE_SPEC)
sys.modules[_PROVENANCE_SPEC.name] = PROVENANCE
_PROVENANCE_SPEC.loader.exec_module(PROVENANCE)


class _FakeDistribution:
    def __init__(
        self,
        root: Path,
        *,
        direct_url: dict[str, object] | None,
    ) -> None:
        self.metadata = {"Name": "helixweave"}
        self.files: tuple[Path, ...] = ()
        self._root = root
        self._direct_url = direct_url

    def locate_file(self, _path: str) -> Path:
        return self._root

    def read_text(self, name: str) -> str | None:
        if name == "top_level.txt":
            return "encode_pipeline\n"
        if name == "METADATA":
            return "Metadata-Version: 2.4\nName: helixweave\n"
        if name == "direct_url.json" and self._direct_url is not None:
            return json.dumps(self._direct_url)
        return None


def _coverage_bootstrap() -> str:
    return (
        "import os, sys\n"
        "if os.environ.get('COVERAGE_PROCESS_START') or "
        "os.environ.get('COVERAGE_PROCESS_CONFIG'):\n"
        f"    sys.path.insert(0, {str(COVERAGE_SUPPORT_ROOT)!r})\n"
        "    import coverage\n"
        "    coverage.process_startup(slug='pth')\n"
        f"    sys.path.remove({str(COVERAGE_SUPPORT_ROOT)!r})\n"
    )


def _create_checkout(root: Path) -> Path:
    package = root / "src" / "encode_pipeline"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        '__version__ = "0.3.0"\n',
        encoding="utf-8",
    )
    return root


def _file_url(path: Path) -> str:
    return f"file://{quote(str(path.resolve()))}"


def _write_distribution(
    site_packages: Path,
    *,
    source_root: Path,
    distribution_name: str = "helixweave",
    direct_root: Path | None = None,
    direct_url_text: str | None = None,
) -> None:
    site_packages.mkdir(parents=True, exist_ok=True)
    normalized = distribution_name.replace("-", "_")
    metadata = site_packages / f"{normalized}-0.3.0.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        f"Metadata-Version: 2.4\nName: {distribution_name}\nVersion: 0.3.0\n",
        encoding="utf-8",
    )
    (metadata / "top_level.txt").write_text(
        "encode_pipeline\n",
        encoding="utf-8",
    )
    if direct_url_text is None:
        direct_url_text = json.dumps(
            {
                "dir_info": {"editable": True},
                "url": _file_url(direct_root or source_root),
            }
        )
    (metadata / "direct_url.json").write_text(direct_url_text, encoding="utf-8")
    (site_packages / f"__editable__.{normalized}-0.3.0.pth").write_text(
        f"{(source_root / 'src').resolve()}\n",
        encoding="utf-8",
    )


def _write_installed_distribution(site_packages: Path) -> None:
    package = site_packages / "encode_pipeline"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        '__version__ = "0.3.0"\n',
        encoding="utf-8",
    )
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
                "archive_info": {"hash": "sha256=" + ("0" * 64)},
                "url": "file:///artifact/helixweave-0.3.0-py3-none-any.whl",
            }
        ),
        encoding="utf-8",
    )


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


def _run_installed_guard(
    python: Path,
    *,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    child_environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONPATH", "PYTHONHOME"}
    }
    child_environment["PYTHONDONTWRITEBYTECODE"] = "1"
    if environment:
        child_environment.update(environment)
    bootstrap = (
        f"{_coverage_bootstrap()}"
        "import runpy\n"
        f"sys.argv = [{str(GUARD_SCRIPT)!r}, 'installed-artifact']\n"
        f"runpy.run_path({str(GUARD_SCRIPT)!r}, run_name='__main__')\n"
    )
    return subprocess.run(
        [str(python), "-c", bootstrap],
        cwd=python.parent,
        env=child_environment,
        capture_output=True,
        text=True,
        check=False,
    )


def _run_with_site(
    code: str,
    *,
    site_packages: Path,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    child_environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONPATH", "PYTHONHOME"}
    }
    child_environment["PYTHONDONTWRITEBYTECODE"] = "1"
    if environment:
        child_environment.update(environment)
    bootstrap = (
        f"{_coverage_bootstrap()}"
        "import site\n"
        f"site.addsitedir({str(site_packages)!r})\n"
        f"sys.path.insert(0, {str(REPOSITORY_ROOT / 'scripts')!r})\n"
        f"{code}\n"
    )
    return subprocess.run(
        [sys.executable, "-S", "-c", bootstrap],
        cwd=site_packages.parent,
        env=child_environment,
        capture_output=True,
        text=True,
        check=False,
    )


def _run_guard(
    site_packages: Path,
    mode: str,
    *,
    repository_root: Path | None = None,
    before_guard: str = "",
) -> subprocess.CompletedProcess[str]:
    arguments = [mode]
    if repository_root is not None:
        arguments.extend(("--repository-root", str(repository_root)))
    return _run_with_site(
        (
            f"{before_guard}\n"
            "from source_provenance import main\n"
            f"raise SystemExit(main({arguments!r}))"
        ),
        site_packages=site_packages,
    )


def test_checkout_mode_accepts_only_the_requested_checkout(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


def test_checkout_mode_rejects_a_sibling_checkout(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=stale)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[module_search_location_mismatch]" in result.stderr


def test_checkout_mode_rejects_module_metadata_disagreement(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=current,
        direct_root=stale,
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_source_mismatch]" in result.stderr


def test_checkout_mode_rejects_multiple_package_search_locations(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    outside = _create_checkout(tmp_path / "outside")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    before_guard = f"""
from importlib.machinery import ModuleSpec
class MultipleLocationFinder:
    def find_spec(self, fullname, path=None, target=None):
        if fullname != "encode_pipeline":
            return None
        spec = ModuleSpec(
            fullname,
            loader=None,
            origin={str(current / "src" / "encode_pipeline" / "__init__.py")!r},
            is_package=True,
        )
        spec.submodule_search_locations = [
            {str(current / "src" / "encode_pipeline")!r},
            {str(outside / "src" / "encode_pipeline")!r},
        ]
        return spec
sys.meta_path.insert(0, MultipleLocationFinder())
"""

    result = _run_guard(
        site_packages,
        "checkout",
        repository_root=current,
        before_guard=before_guard,
    )

    assert result.returncode == 2
    assert "[module_search_location_mismatch]" in result.stderr


def test_checkout_mode_rejects_a_symlink_escape(tmp_path: Path) -> None:
    current = tmp_path / "current"
    (current / "src").mkdir(parents=True)
    outside = _create_checkout(tmp_path / "outside")
    (current / "src" / "encode_pipeline").symlink_to(
        outside / "src" / "encode_pipeline",
        target_is_directory=True,
    )
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[source_root_invalid]" in result.stderr


def test_checkout_mode_fails_stably_for_missing_metadata(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    site_packages.mkdir()
    (site_packages / "current.pth").write_text(
        f"{current / 'src'}\n",
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_missing]" in result.stderr


def test_checkout_mode_fails_stably_for_damaged_metadata(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=current,
        direct_url_text="{not-json",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr


def test_checkout_mode_rejects_stale_known_metadata_without_inventory(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    stale_metadata = site_packages / "encode_pipeline-0.2.0.dist-info"
    stale_metadata.mkdir()
    (stale_metadata / "METADATA").write_text(
        "Metadata-Version: 2.4\nName: encode-pipeline\nVersion: 0.2.0\n",
        encoding="utf-8",
    )
    (stale_metadata / "direct_url.json").write_text(
        json.dumps(
            {
                "dir_info": {"editable": True},
                "url": _file_url(stale),
            }
        ),
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_identity_mismatch]" in result.stderr
    assert str(stale) not in result.stderr


def test_checkout_mode_rejects_known_metadata_without_name_or_inventory(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    damaged_metadata = site_packages / "encode_pipeline-0.2.0.dist-info"
    damaged_metadata.mkdir()

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert "Traceback" not in result.stderr
    assert str(damaged_metadata) not in result.stderr


def test_checkout_mode_fails_stably_for_damaged_record(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    site_packages.mkdir()
    (site_packages / "current.pth").write_text(
        f"{current / 'src'}\n",
        encoding="utf-8",
    )
    metadata = site_packages / "helixweave-0.3.0.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        "Metadata-Version: 2.4\nName: helixweave\nVersion: 0.3.0\n",
        encoding="utf-8",
    )
    (metadata / "RECORD").write_bytes(b"\xff")
    (metadata / "direct_url.json").write_text(
        json.dumps(
            {
                "dir_info": {"editable": True},
                "url": _file_url(current),
            }
        ),
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert "Traceback" not in result.stderr
    assert str(metadata) not in result.stderr


def test_checkout_mode_rejects_source_metadata_pointing_to_sibling(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    source_metadata = current / "src" / "helixweave.egg-info"
    source_metadata.mkdir()
    (source_metadata / "PKG-INFO").write_text(
        "Metadata-Version: 2.4\nName: helixweave\nVersion: 0.3.0\n",
        encoding="utf-8",
    )
    (source_metadata / "top_level.txt").write_text(
        "encode_pipeline\n",
        encoding="utf-8",
    )
    (source_metadata / "direct_url.json").write_text(
        json.dumps(
            {
                "dir_info": {"editable": True},
                "url": _file_url(stale),
            }
        ),
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_source_mismatch]" in result.stderr
    assert str(stale) not in result.stderr


def test_checkout_errors_do_not_disclose_private_paths(tmp_path: Path) -> None:
    private_root = _create_checkout(tmp_path / "private-current")
    stale = _create_checkout(tmp_path / "private-stale")
    site_packages = tmp_path / "private-site-packages"
    _write_distribution(site_packages, source_root=stale)

    result = _run_guard(
        site_packages,
        "checkout",
        repository_root=private_root,
    )

    assert result.returncode == 2
    assert str(private_root) not in result.stderr
    assert str(stale) not in result.stderr
    assert str(site_packages) not in result.stderr


def test_installed_mode_rejects_an_external_editable_path(tmp_path: Path) -> None:
    external = _create_checkout(tmp_path / "external")
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_distribution(site_packages, source_root=external)

    result = _run_installed_guard(python)

    assert result.returncode == 2
    assert "[installed_editable_forbidden]" in result.stderr
    assert str(external) not in result.stderr


def test_installed_mode_rejects_checkout_mixed_with_installed_artifact(
    tmp_path: Path,
) -> None:
    checkout = _create_checkout(tmp_path / "checkout")
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(site_packages)

    result = _run_installed_guard(
        python,
        environment={"PYTHONPATH": str(checkout / "src")},
    )

    assert result.returncode == 2
    assert "[installed_source_mismatch]" in result.stderr
    assert str(checkout) not in result.stderr


def test_installed_mode_rejects_noneditable_directory_install(
    tmp_path: Path,
    monkeypatch,
) -> None:
    source = _create_checkout(tmp_path / "source")
    site_packages = tmp_path / "venv" / "site-packages"
    package = site_packages / "encode_pipeline"
    package.mkdir(parents=True)
    initializer = package / "__init__.py"
    initializer.write_text("", encoding="utf-8")
    distribution = _FakeDistribution(
        site_packages,
        direct_url={
            "dir_info": {},
            "url": _file_url(source),
        },
    )
    spec = SimpleNamespace(
        origin=str(initializer),
        submodule_search_locations=[str(package)],
    )
    monkeypatch.setattr(
        PROVENANCE,
        "_installed_site_roots",
        lambda: (site_packages,),
    )
    monkeypatch.setattr(
        PROVENANCE,
        "_claiming_distributions",
        lambda: (distribution,),
    )
    monkeypatch.setattr(PROVENANCE, "_package_spec", lambda: spec)

    with pytest.raises(PROVENANCE.SourceProvenanceError) as error:
        PROVENANCE.verify_installed_artifact()
    assert error.value.reason_code == "installed_source_mismatch"
    assert str(source) not in str(error.value)


def test_installed_mode_accepts_one_current_site_package_in_process(
    tmp_path: Path,
    monkeypatch,
) -> None:
    site_packages = tmp_path / "venv" / "site-packages"
    package = site_packages / "encode_pipeline"
    package.mkdir(parents=True)
    initializer = package / "__init__.py"
    initializer.write_text("", encoding="utf-8")
    distribution = _FakeDistribution(
        site_packages,
        direct_url={
            "archive_info": {"hash": "sha256=" + ("0" * 64)},
            "url": "file:///artifact/helixweave.whl",
        },
    )
    spec = SimpleNamespace(
        origin=str(initializer),
        submodule_search_locations=[str(package)],
    )
    monkeypatch.setattr(
        PROVENANCE,
        "_installed_site_roots",
        lambda: (site_packages,),
    )
    monkeypatch.setattr(
        PROVENANCE,
        "_claiming_distributions",
        lambda: (distribution,),
    )
    monkeypatch.setattr(PROVENANCE, "_package_spec", lambda: spec)

    PROVENANCE.verify_installed_artifact()


def test_installed_site_roots_require_the_current_isolated_prefix(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prefix = tmp_path / "venv"
    site_packages = prefix / "lib" / "site-packages"
    site_packages.mkdir(parents=True)
    monkeypatch.setattr(PROVENANCE.sys, "prefix", str(prefix))
    monkeypatch.setattr(PROVENANCE.sys, "base_prefix", str(tmp_path / "base"))
    monkeypatch.setattr(
        PROVENANCE.sysconfig,
        "get_path",
        lambda _name: str(site_packages),
    )

    assert PROVENANCE._installed_site_roots() == (site_packages,)

    monkeypatch.setattr(PROVENANCE.sys, "base_prefix", str(prefix))
    with pytest.raises(PROVENANCE.SourceProvenanceError) as error:
        PROVENANCE._installed_site_roots()
    assert error.value.reason_code == "installed_environment_invalid"


def test_guard_cli_dispatch_and_public_failure_are_stable(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    calls: list[tuple[str, Path | None]] = []
    monkeypatch.setattr(
        PROVENANCE,
        "verify_checkout",
        lambda root: calls.append(("checkout", root)),
    )
    monkeypatch.setattr(
        PROVENANCE,
        "verify_installed_artifact",
        lambda: calls.append(("installed-artifact", None)),
    )

    assert (
        PROVENANCE.main(
            ["checkout", "--repository-root", str(tmp_path)],
        )
        == 0
    )
    assert PROVENANCE.main(["installed-artifact"]) == 0
    assert calls == [
        ("checkout", tmp_path),
        ("installed-artifact", None),
    ]

    failure = PROVENANCE.SourceProvenanceError(
        "installed_source_mismatch",
        "use a clean environment",
    )
    monkeypatch.setattr(
        PROVENANCE,
        "verify_installed_artifact",
        lambda: (_ for _ in ()).throw(failure),
    )
    assert PROVENANCE.main(["installed-artifact"]) == 2
    assert capsys.readouterr().err == f"{failure}\n"

    monkeypatch.setattr(
        PROVENANCE,
        "verify_checkout",
        lambda _root: (_ for _ in ()).throw(failure),
    )
    with pytest.raises(SystemExit, match="2"):
        PROVENANCE.require_checkout(tmp_path)
    assert capsys.readouterr().err == f"{failure}\n"


def test_pytest_fails_before_collection_imports_product(
    tmp_path: Path,
) -> None:
    fake_repository = _create_checkout(tmp_path / "private-current")
    (fake_repository / "scripts").mkdir()
    shutil.copy2(
        GUARD_SCRIPT,
        fake_repository / "scripts" / "source_provenance.py",
    )
    fake_test_root = fake_repository / "test"
    fake_test_root.mkdir()
    shutil.copy2(
        REPOSITORY_ROOT / "test" / "conftest.py",
        fake_test_root / "conftest.py",
    )
    marker = tmp_path / "collection-imported"
    (fake_test_root / "test_import_marker.py").write_text(
        "import os\n"
        "from pathlib import Path\n"
        "Path(os.environ['IMPORT_MARKER']).write_text('imported')\n",
        encoding="utf-8",
    )
    stale = _create_checkout(tmp_path / "private-stale")
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONPATH", "PYTHONHOME"}
    }
    environment.update(
        {
            "IMPORT_MARKER": str(marker),
            "PYTHONPATH": str(stale / "src"),
            "PYTHONDONTWRITEBYTECODE": "1",
        }
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "--collect-only",
            str(fake_test_root / "test_import_marker.py"),
            "-p",
            "no:cacheprovider",
        ],
        cwd=fake_repository,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 2
    assert "[module_search_location_mismatch]" in result.stderr
    assert not marker.exists()
    assert str(fake_repository) not in result.stderr
    assert str(stale) not in result.stderr


def test_validator_rejects_stale_editable_before_product_import(
    tmp_path: Path,
) -> None:
    stale = _create_checkout(tmp_path / "stale")
    validate_module = stale / "src" / "encode_pipeline" / "cli" / "validate.py"
    validate_module.parent.mkdir()
    (validate_module.parent / "__init__.py").write_text("", encoding="utf-8")
    validate_module.write_text(
        "import os\n"
        "from pathlib import Path\n"
        "def main():\n"
        "    Path(os.environ['IMPORT_MARKER']).write_text('stale')\n",
        encoding="utf-8",
    )
    marker = tmp_path / "validator-imported"
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=stale,
        distribution_name="encode-pipeline",
    )

    result = _run_with_site(
        (
            "import runpy\n"
            f"sys.argv = [{str(VALIDATOR_SCRIPT)!r}, '--help']\n"
            f"runpy.run_path({str(VALIDATOR_SCRIPT)!r}, run_name='__main__')"
        ),
        site_packages=site_packages,
        environment={"IMPORT_MARKER": str(marker)},
    )

    assert result.returncode == 2
    assert "[distribution_identity_mismatch]" in result.stderr
    assert not marker.exists()


def test_openapi_rejects_stale_editable_before_product_import(
    tmp_path: Path,
) -> None:
    stale = _create_checkout(tmp_path / "stale")
    api_main = stale / "src" / "encode_pipeline" / "api" / "main.py"
    api_main.parent.mkdir()
    (api_main.parent / "__init__.py").write_text("", encoding="utf-8")
    api_main.write_text(
        "import os\n"
        "from pathlib import Path\n"
        "Path(os.environ['IMPORT_MARKER']).write_text('stale')\n"
        "raise RuntimeError('stale product imported')\n",
        encoding="utf-8",
    )
    marker = tmp_path / "openapi-imported"
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=stale,
        distribution_name="encode-pipeline",
    )

    result = _run_with_site(
        (
            "import runpy\n"
            f"sys.argv = [{str(OPENAPI_SCRIPT)!r}]\n"
            f"runpy.run_path({str(OPENAPI_SCRIPT)!r}, run_name='__main__')"
        ),
        site_packages=site_packages,
        environment={"IMPORT_MARKER": str(marker)},
    )

    assert result.returncode == 2
    assert "[distribution_identity_mismatch]" in result.stderr
    assert not marker.exists()


def test_local_platform_rejects_stale_editable_metadata_before_product_import(
    tmp_path: Path,
) -> None:
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=stale,
        distribution_name="encode-pipeline",
    )

    result = _run_with_site(
        (
            "import runpy\n"
            f"sys.argv = [{str(LOCAL_PLATFORM_SCRIPT)!r}, '--help']\n"
            f"runpy.run_path({str(LOCAL_PLATFORM_SCRIPT)!r}, run_name='__main__')"
        ),
        site_packages=site_packages,
    )

    assert result.returncode == 2
    assert "[distribution_identity_mismatch]" in result.stderr
