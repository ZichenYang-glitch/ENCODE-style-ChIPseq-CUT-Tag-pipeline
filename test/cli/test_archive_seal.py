"""Operator CLI contract for sealing an exact track version."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import pytest

from encode_pipeline.cli import app
from encode_pipeline.cli import archive as archive_cli
from encode_pipeline.persistence.runtime import DATABASE_URL_ENV
from encode_pipeline.services.archive_sealing import ArchiveSealError
from encode_pipeline.workers.settings import (
    ARCHIVE_ROOT_ENV,
    WorkerSettings,
    load_worker_settings,
)


PROJECT = "prj_" + "a" * 32
EXPERIMENT = "exp_" + "b" * 32
REVISION = "expr_" + "c" * 32
RESULT = "rreg_" + "1" * 32
RESULT_REVISION = "rrev_" + "2" * 32


def _settings(tmp_path: Path) -> WorkerSettings:
    return load_worker_settings(
        {
            DATABASE_URL_ENV: f"sqlite:///{tmp_path / 'platform.db'}",
            "ENCODE_PIPELINE_WORKSPACE_ROOT": str(tmp_path / "workspaces"),
            ARCHIVE_ROOT_ENV: str(tmp_path / "archive"),
        }
    )


def _args(tmp_path: Path, *extra: str) -> list[str]:
    return [
        "seal",
        "registered",
        "--database-url",
        f"sqlite:///{tmp_path / 'platform.db'}",
        "--archive-root",
        str(tmp_path / "archive"),
        "--project-id",
        PROJECT,
        "--experiment-id",
        EXPERIMENT,
        "--experiment-revision-id",
        REVISION,
        "--track-key",
        "signal",
        "--display-name",
        "signal.bw",
        "--assembly",
        "hg38",
        "--track-type",
        "bigwig",
        "--result-id",
        RESULT,
        "--result-revision-id",
        RESULT_REVISION,
        *extra,
    ]


def test_help_states_that_a_filename_is_not_accepted(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as captured:
        archive_cli.main(["seal", "--help"])

    assert captured.value.code == 0
    output = capsys.readouterr().out
    assert "does not accept a filename" in output
    assert "registered" in output
    assert "run" in output
    assert "ENCODE_PIPELINE_ARCHIVE_ROOT" in output


def test_filename_and_latest_flags_are_rejected(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as captured:
        archive_cli.main(_args(tmp_path, "--filename", "signal.bw", "--latest"))

    assert captured.value.code == 2
    error = capsys.readouterr().err
    assert "unrecognized arguments" in error
    assert "--filename" in error
    assert "--latest" in error


def test_registered_stdout_identifies_the_revision_without_a_path(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    secret = "/secret/pool/signal.bw"
    seen: dict[str, WorkerSettings] = {}

    class _Sealer:
        def seal_registered(self, **kwargs: object) -> dict[str, object]:
            assert "path" not in kwargs
            assert kwargs["result_revision_id"] == RESULT_REVISION
            assert kwargs["display_name"] == "signal.bw"
            return {
                "archive_version_id": "arc-1",
                "assembly": "hg38",
                "display_name": "signal.bw",
                "experiment_id": EXPERIMENT,
                "experiment_revision_id": REVISION,
                "has_index": False,
                "members": {
                    "data": {"sha256": "ab" * 32, "size_bytes": 9},
                },
                "path": secret,
                "project_id": PROJECT,
                "source": {
                    "kind": "registered",
                    "result_id": RESULT,
                    "result_revision_id": RESULT_REVISION,
                },
                "track_key": "signal",
                "track_type": "bigwig",
                "version_number": 1,
            }

    def factory(settings: WorkerSettings) -> _Sealer:
        seen["settings"] = settings
        return _Sealer()

    assert archive_cli.main(_args(tmp_path), sealer_factory=factory) == 1
    captured = capsys.readouterr()
    assert secret not in captured.out
    assert secret not in captured.err
    assert json.loads(captured.err)["error"]["code"] == "archive_integrity"


def test_stdout_lists_the_sealed_identity(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    class _Sealer:
        def seal_registered(self, **kwargs: object) -> dict[str, object]:
            return {
                "archive_version_id": "arc-1",
                "assembly": "hg38",
                "display_name": kwargs["display_name"],
                "experiment_id": kwargs["experiment_id"],
                "experiment_revision_id": kwargs["experiment_revision_id"],
                "has_index": False,
                "members": {"data": {"sha256": "ab" * 32, "size_bytes": 9}},
                "project_id": kwargs["project_id"],
                "source": {
                    "kind": "registered",
                    "result_id": kwargs["result_id"],
                    "result_revision_id": kwargs["result_revision_id"],
                },
                "track_key": kwargs["track_key"],
                "track_type": kwargs["track_type"],
                "version_number": 1,
            }

    assert (
        archive_cli.main(_args(tmp_path), sealer_factory=lambda _settings: _Sealer())
        == 0
    )
    receipt = json.loads(capsys.readouterr().out)
    assert receipt["archive_version_id"] == "arc-1"
    assert receipt["version_number"] == 1
    assert receipt["source"]["result_revision_id"] == RESULT_REVISION
    assert receipt["has_index"] is False
    assert set(receipt) == {
        "archive_version_id",
        "assembly",
        "display_name",
        "experiment_id",
        "experiment_revision_id",
        "has_index",
        "members",
        "project_id",
        "source",
        "track_key",
        "track_type",
        "version_number",
    }


def test_archive_root_flag_overrides_the_environment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    selected = tmp_path / "selected-archive"
    ignored = tmp_path / "ignored-archive"
    monkeypatch.setenv(ARCHIVE_ROOT_ENV, str(ignored))
    seen: dict[str, Path | None] = {}
    arguments = _args(tmp_path)
    arguments[arguments.index("--archive-root") + 1] = str(selected)

    def factory(settings: WorkerSettings) -> SimpleNamespace:
        seen["root"] = settings.archive_root
        raise ArchiveSealError("archive_invalid", "stop after settings")

    assert archive_cli.main(arguments, sealer_factory=factory) == 1
    assert seen["root"] == selected


def test_primary_cli_dispatches_archive(monkeypatch: pytest.MonkeyPatch) -> None:
    observed: list[list[str]] = []
    monkeypatch.setattr(
        archive_cli,
        "main",
        lambda arguments, **kwargs: observed.append(list(arguments)) or 11,
    )

    assert app.main(["archive", "seal", "--help"]) == 11
    assert observed == [["seal", "--help"]]


def test_compose_passes_worker_settings_archive_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive = tmp_path / "archive"
    archive.mkdir()
    workspace = tmp_path / "workspaces"
    workspace.mkdir()
    settings = _settings(tmp_path)
    assert settings.archive_root == archive
    seen: dict[str, object] = {}

    class _ArchiveObjects:
        def __init__(self, root: Path) -> None:
            seen["root"] = root

    class _ArchiveRepository:
        def __init__(self, session_factory: object) -> None:
            seen["session_factory"] = session_factory

    objects_module = ModuleType("encode_pipeline.services.archive_objects")
    objects_module.ArchiveObjects = _ArchiveObjects  # type: ignore[attr-defined]
    tracks_module = ModuleType("encode_pipeline.persistence.archived_tracks")
    tracks_module.ArchiveRepository = _ArchiveRepository  # type: ignore[attr-defined]
    monkeypatch.setitem(
        sys.modules,
        "encode_pipeline.services.archive_objects",
        objects_module,
    )
    monkeypatch.setitem(
        sys.modules,
        "encode_pipeline.persistence.archived_tracks",
        tracks_module,
    )

    class _Persistence:
        def __init__(self) -> None:
            self.engine = object()
            self.repository = object()
            self.catalog_repository = object()
            self.input_registry_repository = object()
            self.closed = False

        def close(self) -> None:
            self.closed = True

    persistence = _Persistence()
    monkeypatch.setattr(
        "encode_pipeline.persistence.runtime.open_existing_run_persistence",
        lambda url: persistence,
    )
    sentinel = object()
    monkeypatch.setattr(
        "encode_pipeline.persistence.database.create_session_factory",
        lambda engine: sentinel,
    )

    class _RunService:
        def __init__(self, registry: object, repository: object = None, **kwargs: object) -> None:
            seen["run_repository"] = repository

    class _Downloads:
        def __init__(
            self,
            *,
            run_service: object,
            workspace_root: Path,
            relocation_root: Path | None = None,
        ) -> None:
            seen["workspace"] = workspace_root
            seen["relocation"] = relocation_root

    monkeypatch.setattr("encode_pipeline.services.runs.RunService", _RunService)
    monkeypatch.setattr(
        "encode_pipeline.services.artifact_downloads.ArtifactDownloadService",
        _Downloads,
    )
    monkeypatch.setattr(
        "encode_pipeline.services.defaults.create_default_workflow_registry",
        lambda: object(),
    )

    _sealer, opened = archive_cli.compose_sealer(settings)

    assert seen["root"] is settings.archive_root
    assert seen["workspace"] == workspace
    assert seen["relocation"] is None
    assert seen["session_factory"] is sentinel
    assert seen["run_repository"] is persistence.repository
    assert opened is persistence
    opened.close()
