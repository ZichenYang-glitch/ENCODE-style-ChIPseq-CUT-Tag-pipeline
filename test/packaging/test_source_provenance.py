"""Current-checkout and installed-artifact source provenance contracts."""

from __future__ import annotations

import base64
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from importlib import machinery
from importlib import util as importlib_util
from importlib.util import find_spec
from urllib.parse import quote
import venv

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
GUARD_SCRIPT = REPOSITORY_ROOT / "scripts" / "source_provenance.py"
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
_AUDITED_HOOKS = {
    "a1_coverage.pth": (
        b'import sys; exec(\'import os\\n\\nif os.getenv("COVERAGE_PROCESS_START") '
        b'or os.getenv("COVERAGE_PROCESS_CONFIG"):\\n try:\\n  import coverage\\n '
        b'except:\\n  pass\\n else:\\n  coverage.process_startup(slug="pth")\')\n',
        "coverage",
        "7.15.1",
    ),
    "coloredlogs.pth": (
        b'import os; exec(\'try: __import__("coloredlogs").auto_install() if '
        b'os.environ.get("COLOREDLOGS_AUTO_INSTALL") else None\\nexcept '
        b"ImportError: pass')\n",
        "coloredlogs",
        "15.0.1",
    ),
    "distutils-precedence.pth": (
        b"import os; var = 'SETUPTOOLS_USE_DISTUTILS'; enabled = "
        b"os.environ.get(var, 'local') == 'local'; enabled and "
        b"__import__('_distutils_hack').add_shim(); \n",
        "setuptools",
        "82.0.1",
    ),
    "sphinxcontrib_jsmath-1.0.1-py3.9-nspkg.pth": (
        b"import sys, types, os;p = os.path.join(sys._getframe(1).f_locals"
        b"['sitedir'], *('sphinxcontrib',));importlib = __import__('importlib.util');"
        b"__import__('importlib.machinery');m = sys.modules.setdefault("
        b"'sphinxcontrib', importlib.util.module_from_spec(importlib.machinery."
        b"PathFinder.find_spec('sphinxcontrib', [os.path.dirname(p)])));m = m or "
        b"sys.modules.setdefault('sphinxcontrib', types.ModuleType('sphinxcontrib'));"
        b"mp = (m or []) and m.__dict__.setdefault('__path__',[]);(p not in mp) "
        b"and mp.append(p)\n",
        "sphinxcontrib-jsmath",
        "1.0.1",
    ),
}


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
    root.mkdir(parents=True)
    (root / "pyproject.toml").write_text(
        '[project]\nname = "helixweave"\nversion = "0.3.0"\n',
        encoding="utf-8",
    )
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
    version: str = "0.3.0",
    direct_root: Path | None = None,
    direct_url_text: str | None = None,
    metadata_directory: str | None = None,
    include_inventory: bool = True,
) -> None:
    site_packages.mkdir(parents=True, exist_ok=True)
    normalized = distribution_name.replace("-", "_")
    metadata = site_packages / (
        metadata_directory or f"{normalized}-{version}.dist-info"
    )
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        f"Metadata-Version: 2.4\nName: {distribution_name}\nVersion: {version}\n",
        encoding="utf-8",
    )
    if include_inventory:
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
    pth_name = f"__editable__.{normalized}-{version}.pth"
    (site_packages / pth_name).write_text(
        f"{(source_root / 'src').resolve()}\n",
        encoding="utf-8",
    )
    if include_inventory:
        (metadata / "RECORD").write_text(
            "\n".join(
                (
                    f"{pth_name},,",
                    f"{metadata.name}/METADATA,,",
                    f"{metadata.name}/top_level.txt,,",
                    f"{metadata.name}/direct_url.json,,",
                    f"{metadata.name}/RECORD,,",
                )
            )
            + "\n",
            encoding="utf-8",
        )


def _write_audited_hook_owner(
    site_packages: Path,
    hook_name: str,
    *,
    distribution_name: str | None = None,
    version: str | None = None,
    raw: bytes | None = None,
) -> Path:
    expected_raw, expected_name, expected_version = _AUDITED_HOOKS[hook_name]
    hook_raw = expected_raw if raw is None else raw
    hook = site_packages / hook_name
    hook.write_bytes(hook_raw)
    owner_name = distribution_name or expected_name
    owner_version = version or expected_version
    metadata_name = owner_name.replace("-", "_")
    metadata = site_packages / f"{metadata_name}-{owner_version}.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        f"Metadata-Version: 2.4\nName: {owner_name}\nVersion: {owner_version}\n",
        encoding="utf-8",
    )
    digest = base64.urlsafe_b64encode(hashlib.sha256(hook_raw).digest()).decode("ascii")
    (metadata / "RECORD").write_text(
        f"{hook_name},sha256={digest.rstrip('=')},{len(hook_raw)}\n",
        encoding="utf-8",
    )
    return hook


def _write_conda_hook_inventory(
    prefix: Path,
    site_packages: Path,
    hook_name: str,
    *,
    version: str,
    build: str,
    raw: bytes,
) -> Path:
    metadata_root = prefix / "conda-meta"
    metadata_root.mkdir(exist_ok=True)
    digest = hashlib.sha256(raw).hexdigest()
    metadata = metadata_root / f"setuptools-{version}-{build}.json"
    metadata.write_text(
        json.dumps(
            {
                "name": "setuptools",
                "version": version,
                "build": build,
                "paths_data": {
                    "paths": [
                        {
                            "_path": (site_packages / hook_name)
                            .relative_to(prefix)
                            .as_posix(),
                            "path_type": "hardlink",
                            "sha256": digest,
                            "sha256_in_prefix": digest,
                            "size_in_bytes": len(raw),
                        }
                    ]
                },
            }
        ),
        encoding="utf-8",
    )
    return metadata


def _write_installed_distribution(
    site_packages: Path,
    *,
    include_inventory: bool = True,
    metadata_directory: str = "helixweave-0.3.0.dist-info",
) -> None:
    package = site_packages / "encode_pipeline"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        '__version__ = "0.3.0"\n',
        encoding="utf-8",
    )
    metadata = site_packages / metadata_directory
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
    if include_inventory:
        (metadata / "RECORD").write_text(
            "\n".join(
                (
                    "encode_pipeline/__init__.py,,",
                    "helixweave-0.3.0.dist-info/METADATA,,",
                    "helixweave-0.3.0.dist-info/top_level.txt,,",
                    "helixweave-0.3.0.dist-info/direct_url.json,,",
                    "helixweave-0.3.0.dist-info/RECORD,,",
                )
            )
            + "\n",
            encoding="utf-8",
        )


def _write_unknown_distribution_metadata(
    site_packages: Path,
    *,
    inventory: str,
    metadata_content: str | bytes | None = None,
) -> Path:
    metadata = site_packages / "opaque_owner-9.9.dist-info"
    metadata.mkdir()
    if metadata_content is not None:
        if isinstance(metadata_content, bytes):
            (metadata / "METADATA").write_bytes(metadata_content)
        else:
            (metadata / "METADATA").write_text(
                metadata_content,
                encoding="utf-8",
            )
    if inventory in {"top-level", "both"}:
        (metadata / "top_level.txt").write_text(
            "encode_pipeline\n",
            encoding="utf-8",
        )
    elif inventory == "unrelated":
        (metadata / "top_level.txt").write_text(
            "unrelated_package\n",
            encoding="utf-8",
        )
    if inventory in {"record", "both"}:
        (metadata / "RECORD").write_text(
            "encode_pipeline/__init__.py,,\n",
            encoding="utf-8",
        )
    elif inventory == "editable-mapping":
        (metadata / "RECORD").write_text(
            "opaque-owner.pth,,\n",
            encoding="utf-8",
        )
    elif inventory == "unrelated":
        (metadata / "RECORD").write_text(
            "unrelated_package/__init__.py,,\n",
            encoding="utf-8",
        )
    return metadata


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
    assert mode == "checkout"
    assert repository_root is not None
    return _run_with_site(
        (
            f"{before_guard}\n"
            "from pathlib import Path\n"
            "from source_provenance import SourceProvenanceError, verify_checkout\n"
            "try:\n"
            f"    verify_checkout(Path({str(repository_root)!r}), "
            f"site_roots=(Path({str(site_packages)!r}),))\n"
            "except SourceProvenanceError as error:\n"
            "    print(error, file=sys.stderr)\n"
            "    raise SystemExit(2) from None\n"
        ),
        site_packages=site_packages,
    )


def _run_with_unknown_metadata(
    tmp_path: Path,
    *,
    mode: str,
    inventory: str,
    metadata_content: str | bytes | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    if mode == "checkout":
        current = _create_checkout(tmp_path / "current")
        site_packages = tmp_path / "site-packages"
        _write_distribution(site_packages, source_root=current)
        metadata = _write_unknown_distribution_metadata(
            site_packages,
            inventory=inventory,
            metadata_content=metadata_content,
        )
        result = _run_guard(
            site_packages,
            "checkout",
            repository_root=current,
        )
        return result, metadata
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(site_packages)
    metadata = _write_unknown_distribution_metadata(
        site_packages,
        inventory=inventory,
        metadata_content=metadata_content,
    )
    return _run_installed_guard(python), metadata


def test_checkout_mode_accepts_only_the_requested_checkout(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


@pytest.mark.parametrize("mode", ("checkout", "installed-artifact"))
@pytest.mark.parametrize("inventory", ("top-level", "record", "both"))
def test_namespace_inventory_cannot_hide_an_unknown_claimant_without_metadata(
    tmp_path: Path,
    mode: str,
    inventory: str,
) -> None:
    result, metadata = _run_with_unknown_metadata(
        tmp_path,
        mode=mode,
        inventory=inventory,
    )

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert str(metadata) not in result.stderr
    assert "encode_pipeline/__init__.py" not in result.stderr
    assert "Traceback" not in result.stderr


def test_audited_editable_mapping_cannot_hide_unknown_claimant_metadata(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    metadata = _write_unknown_distribution_metadata(
        site_packages,
        inventory="editable-mapping",
    )
    (site_packages / "opaque-owner.pth").write_text(
        f"{stale / 'src'}\n",
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert str(metadata) not in result.stderr
    assert str(stale) not in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("mode", ("checkout", "installed-artifact"))
@pytest.mark.parametrize(
    "metadata_content",
    (
        b"\xffprivate metadata must not be disclosed",
        (
            "Metadata-Version: 2.4\n"
            "Name: helixweave\n"
            "Name: conflicting-owner\n"
            "Version: 0.3.0\n"
        ),
    ),
)
def test_namespace_inventory_rejects_malformed_or_ambiguous_metadata(
    tmp_path: Path,
    mode: str,
    metadata_content: str | bytes,
) -> None:
    result, metadata = _run_with_unknown_metadata(
        tmp_path,
        mode=mode,
        inventory="both",
        metadata_content=metadata_content,
    )

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert str(metadata) not in result.stderr
    assert "private metadata" not in result.stderr
    assert "conflicting-owner" not in result.stderr
    assert "UnicodeDecodeError" not in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("mode", ("checkout", "installed-artifact"))
def test_valid_claimant_conflicts_with_an_unknown_namespace_owner(
    tmp_path: Path,
    mode: str,
) -> None:
    result, metadata = _run_with_unknown_metadata(
        tmp_path,
        mode=mode,
        inventory="both",
        metadata_content=("Metadata-Version: 2.4\nName: opaque-owner\nVersion: 9.9\n"),
    )

    assert result.returncode == 2
    assert "[distribution_claimant_conflict]" in result.stderr
    assert str(metadata) not in result.stderr
    assert "opaque-owner" not in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("mode", ("checkout", "installed-artifact"))
def test_unrelated_unknown_distribution_without_namespace_inventory_is_ignored(
    tmp_path: Path,
    mode: str,
) -> None:
    result, _metadata = _run_with_unknown_metadata(
        tmp_path,
        mode=mode,
        inventory="unrelated",
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


def test_checkout_accepts_valid_claimant_with_an_unknown_metadata_directory(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=current,
        metadata_directory="opaque_owner-9.9.dist-info",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


def test_installed_artifact_accepts_valid_claimant_with_unknown_directory(
    tmp_path: Path,
) -> None:
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(
        site_packages,
        metadata_directory="opaque_owner-9.9.dist-info",
    )

    result = _run_installed_guard(python)

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


@pytest.mark.parametrize("root_kind", ("missing", "file"))
def test_checkout_mode_rejects_an_invalid_repository_root(
    tmp_path: Path,
    root_kind: str,
) -> None:
    repository_root = tmp_path / "invalid-repository"
    if root_kind == "file":
        repository_root.write_text("not a checkout\n", encoding="utf-8")
    site_packages = tmp_path / "site-packages"
    site_packages.mkdir()

    result = _run_guard(
        site_packages,
        "checkout",
        repository_root=repository_root,
    )

    assert result.returncode == 2
    assert "[repository_root_invalid]" in result.stderr
    assert str(repository_root) not in result.stderr


def test_checkout_mode_rejects_two_physical_claimants_with_the_same_version(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    _write_distribution(
        site_packages,
        source_root=current,
        distribution_name="HelixWeave",
        metadata_directory="HELIXWEAVE-copy-0.3.0.dist-info",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_claimant_conflict]" in result.stderr


def test_checkout_mode_rejects_a_claimant_version_mismatch(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=current,
        version="0.2.0",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_version_mismatch]" in result.stderr


def test_checkout_mode_rejects_missing_namespace_ownership_inventory(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(
        site_packages,
        source_root=current,
        include_inventory=False,
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[namespace_ownership_unproven]" in result.stderr


def test_checkout_mode_rejects_an_orphan_stale_pth(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    (site_packages / "zzz-orphan-stale.pth").write_text(
        f"{stale / 'src'}\n",
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[namespace_mapping_conflict]" in result.stderr
    assert str(stale) not in result.stderr


def test_checkout_mode_accepts_exact_inventoried_nonproduct_hooks(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    for hook_name in _AUDITED_HOOKS:
        _write_audited_hook_owner(site_packages, hook_name)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    ("owner_name", "owner_version", "raw_suffix"),
    (
        ("opaque-owner", None, b""),
        (None, "0.0", b""),
        (None, None, b"# extra startup metadata\n"),
    ),
)
def test_checkout_mode_rejects_a_hook_with_drifted_shape_or_owner(
    tmp_path: Path,
    owner_name: str | None,
    owner_version: str | None,
    raw_suffix: bytes,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "coloredlogs.pth"
    expected_raw = _AUDITED_HOOKS[hook_name][0]
    _write_audited_hook_owner(
        site_packages,
        hook_name,
        distribution_name=owner_name,
        version=owner_version,
        raw=expected_raw + raw_suffix,
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr
    assert str(site_packages) not in result.stderr


def test_checkout_mode_rejects_an_exact_hook_without_owner_inventory(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "coloredlogs.pth"
    (site_packages / hook_name).write_bytes(_AUDITED_HOOKS[hook_name][0])

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr


@pytest.mark.parametrize(
    ("build", "expected_returncode"),
    (
        ("pyh332efcf_0", 0),
        ("unreviewed_0", 2),
    ),
)
def test_checkout_mode_reconciles_exact_conda_hook_inventory(
    tmp_path: Path,
    build: str,
    expected_returncode: int,
) -> None:
    current = _create_checkout(tmp_path / "current")
    prefix = tmp_path / "environment"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "distutils-precedence.pth"
    raw, _, version = _AUDITED_HOOKS[hook_name]
    (site_packages / hook_name).write_bytes(raw)
    _write_conda_hook_inventory(
        prefix,
        site_packages,
        hook_name,
        version=version,
        build=build,
        raw=raw,
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == expected_returncode, result.stderr
    if expected_returncode:
        assert "[pth_mapping_unsafe]" in result.stderr


def test_checkout_mode_uses_the_nearest_conda_environment_inventory(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    root_prefix = tmp_path / "micromamba"
    prefix = root_prefix / "envs" / "ci-fast"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "distutils-precedence.pth"
    raw, _, version = _AUDITED_HOOKS[hook_name]
    (site_packages / hook_name).write_bytes(raw)
    _write_conda_hook_inventory(
        prefix,
        site_packages,
        hook_name,
        version=version,
        build="pyh332efcf_0",
        raw=raw,
    )
    (root_prefix / "conda-meta").mkdir()

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 0, result.stderr


def test_checkout_mode_rejects_duplicate_record_hook_owners(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "coloredlogs.pth"
    _write_audited_hook_owner(site_packages, hook_name)
    original = site_packages / "coloredlogs-15.0.1.dist-info"
    duplicate = site_packages / "coloredlogs-copy.dist-info"
    duplicate.mkdir()
    for filename in ("METADATA", "RECORD"):
        (duplicate / filename).write_bytes((original / filename).read_bytes())

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr


def test_checkout_mode_rejects_cross_inventory_owner_version_conflict(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    prefix = tmp_path / "environment"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "distutils-precedence.pth"
    raw, _, _ = _AUDITED_HOOKS[hook_name]
    _write_audited_hook_owner(
        site_packages,
        hook_name,
        version="83.0.0",
        raw=raw,
    )
    _write_conda_hook_inventory(
        prefix,
        site_packages,
        hook_name,
        version="82.0.1",
        build="pyh332efcf_0",
        raw=raw,
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr


def test_checkout_mode_requires_conda_inventory_when_prefix_is_present(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    prefix = tmp_path / "environment"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    _write_audited_hook_owner(site_packages, "coloredlogs.pth")
    (prefix / "conda-meta").mkdir()

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr


@pytest.mark.parametrize("target_location", ("inside", "outside"))
def test_checkout_mode_rejects_a_conda_metadata_root_symlink(
    tmp_path: Path,
    target_location: str,
) -> None:
    current = _create_checkout(tmp_path / "current")
    prefix = tmp_path / "environment"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    _write_audited_hook_owner(site_packages, "coloredlogs.pth")
    target = (
        prefix / "metadata-real"
        if target_location == "inside"
        else tmp_path / "external-metadata"
    )
    target.mkdir()
    (prefix / "conda-meta").symlink_to(target, target_is_directory=True)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr
    assert str(target) not in result.stderr


def test_checkout_mode_rejects_a_nondirectory_conda_metadata_root(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    prefix = tmp_path / "environment"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    _write_audited_hook_owner(site_packages, "coloredlogs.pth")
    (prefix / "conda-meta").write_text("not a metadata directory\n", encoding="utf-8")

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr


def test_checkout_mode_rejects_broken_record_symlink_with_valid_conda_owner(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    prefix = tmp_path / "environment"
    site_packages = prefix / "lib" / "python3.12" / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook_name = "distutils-precedence.pth"
    raw, _, version = _AUDITED_HOOKS[hook_name]
    _write_audited_hook_owner(site_packages, hook_name, raw=raw)
    owner_record = site_packages / f"setuptools-{version}.dist-info" / "RECORD"
    owner_record.unlink()
    owner_record.symlink_to(tmp_path / "missing-record")
    _write_conda_hook_inventory(
        prefix,
        site_packages,
        hook_name,
        version=version,
        build="pyh332efcf_0",
        raw=raw,
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr


def test_checkout_mode_rejects_an_obfuscated_executable_pth(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    (site_packages / "a1_coverage.pth").write_text(
        'import sys; __import__("encode_" + "pipeline") '
        'if "never" in sys.modules else None\n',
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr
    assert str(site_packages) not in result.stderr


def test_checkout_mode_rejects_a_pth_file_symlink_escape(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    pth = site_packages / "__editable__.helixweave-0.3.0.pth"
    external_pth = tmp_path / "external-editable.pth"
    pth.replace(external_pth)
    pth.symlink_to(external_pth)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr
    assert str(external_pth) not in result.stderr


@pytest.mark.parametrize("target_kind", ("directory", "archive"))
def test_checkout_mode_rejects_an_extra_existing_pth_target(
    tmp_path: Path,
    target_kind: str,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    target = tmp_path / f"unrelated-{target_kind}"
    if target_kind == "directory":
        target.mkdir()
    else:
        target.write_bytes(b"not a trusted import archive")
    (site_packages / "orphan-path.pth").write_text(
        f"{target}\n",
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[namespace_mapping_conflict]" in result.stderr
    assert str(target) not in result.stderr


@pytest.mark.parametrize(
    "hook_name",
    (
        "sitecustomize",
        "sitecustomize.pyc",
        f"usercustomize{machinery.EXTENSION_SUFFIXES[-1]}",
    ),
)
def test_checkout_mode_rejects_every_direct_startup_hook_shape(
    tmp_path: Path,
    hook_name: str,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    hook = site_packages / hook_name
    if "." in hook_name:
        hook.write_bytes(b"must never load")
    else:
        hook.mkdir()
        (hook / "__init__.py").write_text(
            "raise AssertionError('must never import')\n",
            encoding="utf-8",
        )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[startup_hook_unsafe]" in result.stderr
    assert str(hook) not in result.stderr


def test_checkout_mode_rejects_a_stale_editable_finder(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    (site_packages / "__editable___stale_finder.py").write_text(
        "raise AssertionError('finder content must never execute')\n",
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[pth_mapping_unsafe]" in result.stderr
    assert "AssertionError" not in result.stderr


def test_checkout_mode_rejects_a_second_symlink_alias_mapping(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    source_alias = tmp_path / "source-alias"
    source_alias.symlink_to(current / "src", target_is_directory=True)
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    (site_packages / "zzz-source-alias.pth").write_text(
        f"{source_alias}\n",
        encoding="utf-8",
    )

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[namespace_mapping_conflict]" in result.stderr


def test_checkout_mode_rejects_a_sibling_checkout(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    stale = _create_checkout(tmp_path / "stale")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=stale)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_source_mismatch]" in result.stderr


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
    assert "[distribution_claimant_conflict]" in result.stderr
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


@pytest.mark.parametrize(
    "metadata_name",
    ("METADATA", "top_level.txt", "RECORD", "direct_url.json"),
)
def test_checkout_mode_rejects_metadata_file_symlink_escape(
    tmp_path: Path,
    metadata_name: str,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    metadata = site_packages / "helixweave-0.3.0.dist-info"
    external_file = tmp_path / f"external-{metadata_name}"
    (metadata / metadata_name).replace(external_file)
    (metadata / metadata_name).symlink_to(external_file)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert str(external_file) not in result.stderr


def test_checkout_mode_rejects_metadata_directory_symlink_escape(
    tmp_path: Path,
) -> None:
    current = _create_checkout(tmp_path / "current")
    site_packages = tmp_path / "site-packages"
    _write_distribution(site_packages, source_root=current)
    metadata = site_packages / "helixweave-0.3.0.dist-info"
    external_metadata = tmp_path / "external-metadata"
    metadata.replace(external_metadata)
    metadata.symlink_to(external_metadata, target_is_directory=True)

    result = _run_guard(site_packages, "checkout", repository_root=current)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert str(external_metadata) not in result.stderr


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


def test_installed_mode_rejects_an_artifact_without_ownership_inventory(
    tmp_path: Path,
) -> None:
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(site_packages, include_inventory=False)

    result = _run_installed_guard(python)

    assert result.returncode == 2
    assert "[namespace_ownership_unproven]" in result.stderr


def test_installed_mode_rejects_metadata_directory_symlink_escape(
    tmp_path: Path,
) -> None:
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(site_packages)
    metadata = site_packages / "helixweave-0.3.0.dist-info"
    external_metadata = tmp_path / "external-metadata"
    metadata.replace(external_metadata)
    metadata.symlink_to(external_metadata, target_is_directory=True)

    result = _run_installed_guard(python)

    assert result.returncode == 2
    assert "[distribution_metadata_invalid]" in result.stderr
    assert str(external_metadata) not in result.stderr


def test_installed_mode_accepts_a_record_owned_artifact(tmp_path: Path) -> None:
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(site_packages)

    result = _run_installed_guard(python)

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


def test_installed_mode_rejects_noneditable_directory_install(
    tmp_path: Path,
) -> None:
    source = _create_checkout(tmp_path / "source")
    python, site_packages = _create_venv(tmp_path / "venv")
    _write_installed_distribution(site_packages)
    (site_packages / "helixweave-0.3.0.dist-info" / "direct_url.json").write_text(
        json.dumps(
            {
                "dir_info": {},
                "url": _file_url(source),
            }
        ),
        encoding="utf-8",
    )

    result = _run_installed_guard(python)

    assert result.returncode == 2
    assert "[installed_source_mismatch]" in result.stderr
    assert str(source) not in result.stderr


def test_clean_bootstrap_discovers_a_venv_site_root_under_no_site(
    tmp_path: Path,
) -> None:
    python, site_packages = _create_venv(tmp_path / "venv")
    code = (
        "import importlib.util,sys\n"
        f"spec=importlib.util.spec_from_file_location('provenance',{str(GUARD_SCRIPT)!r})\n"
        "module=importlib.util.module_from_spec(spec)\n"
        "sys.modules[spec.name]=module\n"
        "spec.loader.exec_module(module)\n"
        "print(module.environment_site_roots()[0])\n"
    )
    result = subprocess.run(
        [str(python), "-I", "-S", "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert Path(result.stdout.strip()) == site_packages


def test_checkout_layout_loads_without_python_311_tomllib(tmp_path: Path) -> None:
    current = _create_checkout(tmp_path / "current")
    code = f"""
import importlib.abc
import importlib.util
from pathlib import Path
import sys

class NoTomllib(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "tomllib":
            raise ModuleNotFoundError("simulated Python 3.10 stdlib")
        return None

sys.meta_path.insert(0, NoTomllib())
spec = importlib.util.spec_from_file_location("provenance_310", {str(GUARD_SCRIPT)!r})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(module._checkout_layout(Path({str(current)!r})).version)
"""
    result = subprocess.run(
        [sys.executable, "-I", "-S", "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "0.3.0"


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
