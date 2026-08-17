"""Machine-readable local doctor recovery visibility contracts."""

from __future__ import annotations

from contextlib import contextmanager
import json
from types import SimpleNamespace

import pytest

from encode_pipeline.cli import admin
from encode_pipeline.cli import local_platform
from encode_pipeline.cli.local_platform import (
    EnvironmentCheck,
    RecoveryDoctorCheck,
    RecoveryDoctorCounts,
    RecoveryDoctorStatus,
    WorkflowCheck,
    run_recovery_doctor,
)


def test_recovery_doctor_missing_database_is_not_configured_and_read_only(
    tmp_path,
) -> None:
    database_path = tmp_path / "missing" / "platform.db"

    check = run_recovery_doctor(
        database_url=f"sqlite:///{database_path}",
        redis_url="redis://127.0.0.1:6379/0",
        queue_name="default",
    )

    assert check == RecoveryDoctorCheck(
        RecoveryDoctorStatus.NOT_CONFIGURED,
        "RUN_RECOVERY_DATABASE_NOT_CONFIGURED",
        RecoveryDoctorCounts(),
    )
    assert not database_path.parent.exists()


@pytest.mark.parametrize(
    ("summary", "expected_status", "expected_reason"),
    [
        (
            SimpleNamespace(
                runs_examined=2,
                database_gaps=0,
                queue_gaps=0,
                callback_gaps=0,
                result_indexing_gaps=0,
                cleanup_gaps=0,
                queue_unavailable=False,
            ),
            RecoveryDoctorStatus.READY,
            "RUN_RECOVERY_READY",
        ),
        (
            SimpleNamespace(
                runs_examined=2,
                database_gaps=0,
                queue_gaps=1,
                callback_gaps=0,
                result_indexing_gaps=0,
                cleanup_gaps=0,
                queue_unavailable=False,
            ),
            RecoveryDoctorStatus.ATTENTION_REQUIRED,
            "RUN_RECOVERY_ATTENTION_REQUIRED",
        ),
        (
            SimpleNamespace(
                runs_examined=2,
                database_gaps=0,
                queue_gaps=2,
                callback_gaps=0,
                result_indexing_gaps=0,
                cleanup_gaps=0,
                queue_unavailable=True,
            ),
            RecoveryDoctorStatus.UNAVAILABLE,
            "RUN_RECOVERY_QUEUE_UNAVAILABLE",
        ),
    ],
)
def test_recovery_doctor_maps_typed_summary_to_stable_status(
    tmp_path,
    monkeypatch,
    summary,
    expected_status,
    expected_reason,
) -> None:
    database_path = tmp_path / "platform.db"
    database_path.touch()
    observed = []

    class Recovery:
        def summarize(self):
            return summary

    @contextmanager
    def factory(database_url, *, read_only, redis_url, queue_name):
        observed.append((database_url, read_only, redis_url, queue_name))
        yield Recovery()

    monkeypatch.setattr(admin, "_open_run_recovery", factory)

    check = run_recovery_doctor(
        database_url=f"sqlite:///{database_path}",
        redis_url="redis://127.0.0.1:6379/0",
        queue_name="default",
    )

    assert check.status is expected_status
    assert check.reason_code == expected_reason
    assert check.counts.queue_gaps == summary.queue_gaps
    assert observed == [
        (
            f"sqlite:///{database_path}",
            True,
            "redis://127.0.0.1:6379/0",
            "default",
        )
    ]


def test_json_doctor_adds_safe_recovery_counts_without_runtime_mutation(
    tmp_path, monkeypatch, capsys
) -> None:
    runtime_root = tmp_path / "runtime"
    calls: list[object] = []
    monkeypatch.setattr(
        local_platform,
        "run_environment_doctor",
        lambda _root, *, redis_url: (
            calls.append(("environment", redis_url))
            or (EnvironmentCheck("Python", "3.12.13"),)
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "run_workflow_doctor",
        lambda _root, *, environ: (
            calls.append("workflows")
            or (
                WorkflowCheck(
                    "bulk-rnaseq",
                    "Bulk RNA-seq",
                    "available",
                    "not_configured",
                    "WORKFLOW_EXECUTION_NOT_CONFIGURED",
                ),
            )
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "run_recovery_doctor",
        lambda *, database_url, redis_url, queue_name: (
            calls.append(("recovery", database_url, redis_url, queue_name))
            or SimpleNamespace(
                status="attention_required",
                reason_code="RUN_RECOVERY_ATTENTION_REQUIRED",
                counts=SimpleNamespace(
                    runs_examined=7,
                    database_gaps=1,
                    queue_gaps=2,
                    callback_gaps=3,
                    result_indexing_gaps=4,
                    cleanup_gaps=5,
                ),
            )
        ),
        raising=False,
    )
    monkeypatch.setattr(
        local_platform,
        "_port_available",
        lambda *_args: calls.append("port") or True,
    )

    assert (
        local_platform.main(["--doctor", "--json", "--runtime-root", str(runtime_root)])
        == 0
    )

    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == (
        '{"environment":[{"name":"Python","version":"3.12.13"}],'
        '"run_recovery":{"counts":{"callback_gaps":3,"cleanup_gaps":5,'
        '"database_gaps":1,"queue_gaps":2,"result_indexing_gaps":4,'
        '"runs_examined":7},"reason_code":"RUN_RECOVERY_ATTENTION_REQUIRED",'
        '"status":"attention_required"},"schema_version":"1.0.0",'
        '"workflows":[{"authoring":"available","execution":"not_configured",'
        '"name":"Bulk RNA-seq","reason_code":'
        '"WORKFLOW_EXECUTION_NOT_CONFIGURED","workflow_id":"bulk-rnaseq"}]}\n'
    )
    assert calls == [
        ("environment", "redis://127.0.0.1:6379/0"),
        "workflows",
        (
            "recovery",
            f"sqlite:///{runtime_root / 'platform.db'}",
            "redis://127.0.0.1:6379/0",
            "encode-pipeline-demo",
        ),
    ]
    assert not runtime_root.exists()


def test_json_doctor_collects_independent_failures_without_disclosure(
    tmp_path, monkeypatch, capsys
) -> None:
    private_value = str(tmp_path / "private" / "runtime")
    monkeypatch.setattr(
        local_platform,
        "run_environment_doctor",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError(private_value)),
    )
    monkeypatch.setattr(
        local_platform,
        "run_workflow_doctor",
        lambda *_args, **_kwargs: (
            WorkflowCheck(
                "bulk-rnaseq",
                "Bulk RNA-seq",
                "available",
                "available",
                "WORKFLOW_EXECUTION_READY",
            ),
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "run_recovery_doctor",
        lambda **_kwargs: RecoveryDoctorCheck(
            RecoveryDoctorStatus.NOT_CONFIGURED,
            "RUN_RECOVERY_DATABASE_NOT_CONFIGURED",
            RecoveryDoctorCounts(),
        ),
    )

    assert local_platform.main(["--doctor", "--json"]) == 1

    captured = capsys.readouterr()
    assert captured.err == ""
    assert private_value not in captured.out
    assert json.loads(captured.out) == {
        "schema_version": "1.0.0",
        "environment": [],
        "workflows": [
            {
                "workflow_id": "bulk-rnaseq",
                "name": "Bulk RNA-seq",
                "authoring": "available",
                "execution": "available",
                "reason_code": "WORKFLOW_EXECUTION_READY",
            }
        ],
        "run_recovery": {
            "status": "not_configured",
            "reason_code": "RUN_RECOVERY_DATABASE_NOT_CONFIGURED",
            "counts": {
                "runs_examined": 0,
                "database_gaps": 0,
                "queue_gaps": 0,
                "callback_gaps": 0,
                "result_indexing_gaps": 0,
                "cleanup_gaps": 0,
            },
        },
        "errors": [
            {
                "component": "environment",
                "reason_code": "ENVIRONMENT_CHECK_FAILED",
            }
        ],
    }


def test_json_doctor_sanitizes_unexpected_recovery_failure_and_exits_one(
    monkeypatch, capsys
) -> None:
    monkeypatch.setattr(
        local_platform,
        "run_environment_doctor",
        lambda *_args, **_kwargs: (EnvironmentCheck("Python", "3.12.13"),),
    )
    monkeypatch.setattr(
        local_platform,
        "run_workflow_doctor",
        lambda *_args, **_kwargs: (),
    )
    monkeypatch.setattr(
        local_platform,
        "run_recovery_doctor",
        lambda **_kwargs: (_ for _ in ()).throw(
            RuntimeError("redis://user:secret@private-host/0")
        ),
    )

    assert local_platform.main(["--doctor", "--json"]) == 1

    captured = capsys.readouterr()
    assert captured.err == ""
    assert "secret" not in captured.out
    recovery = json.loads(captured.out)["run_recovery"]
    assert recovery == {
        "status": "unavailable",
        "reason_code": "RUN_RECOVERY_INTERNAL_ERROR",
        "counts": {
            "runs_examined": 0,
            "database_gaps": 0,
            "queue_gaps": 0,
            "callback_gaps": 0,
            "result_indexing_gaps": 0,
            "cleanup_gaps": 0,
        },
    }


def test_prose_doctor_does_not_call_recovery_or_change_existing_output(
    tmp_path, monkeypatch, capsys
) -> None:
    monkeypatch.setattr(
        local_platform,
        "run_environment_doctor",
        lambda _root, *, redis_url: (EnvironmentCheck("Python", "3.12.13"),),
    )
    monkeypatch.setattr(
        local_platform,
        "run_workflow_doctor",
        lambda _root, *, environ: (
            WorkflowCheck(
                "bulk-rnaseq",
                "Bulk RNA-seq",
                "available",
                "available",
                "WORKFLOW_EXECUTION_READY",
            ),
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "run_recovery_doctor",
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("prose doctor must remain unchanged")
        ),
        raising=False,
    )

    assert local_platform.main(["--doctor"]) == 0

    assert capsys.readouterr().out == (
        "HelixWeave environment check\n"
        "Python: 3.12.13\n"
        "HelixWeave workflow availability\n"
        "bulk-rnaseq: authoring=available execution=available "
        "reason=WORKFLOW_EXECUTION_READY\n"
    )
