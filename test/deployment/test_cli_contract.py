from __future__ import annotations

from dataclasses import dataclass
import io
import json
from pathlib import Path

import pytest

from encode_pipeline.deployment.cli import (
    DOCTOR_RESULT_SCHEMA,
    EXIT_DATA,
    EXIT_INCOMPATIBLE,
    EXIT_OK,
    EXIT_OPERATION,
    EXIT_UNAVAILABLE,
    EXIT_USAGE,
    OPERATION_RESULT_SCHEMA,
    STATUS_RESULT_SCHEMA,
    VERIFY_RESULT_SCHEMA,
    CommandRequest,
    PublicCommandResult,
    main,
    parse_command,
)
from encode_pipeline.deployment.errors import DeploymentError, fail


IDENTITY = f"sha256-{'a' * 64}"


def _components() -> dict[str, object]:
    return {
        component: {"active": None, "previous": None, "staged": None}
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }


def _services() -> dict[str, object]:
    return {
        service: {"state": "stopped", "identity": None}
        for service in ("api", "worker", "redis", "bulk_rnaseq_docker")
    }


def _status_result() -> dict[str, object]:
    return {
        "schema_version": STATUS_RESULT_SCHEMA,
        "state_identity": IDENTITY,
        "generation": 1,
        "components": _components(),
        "interrupted": False,
        "partial_staging_count": 0,
        "pending_transaction_count": 0,
        "orphaned_deployment_count": 0,
        "database_schema_identity": None,
        "services": _services(),
    }


def _operation_result(command: str, component: str) -> dict[str, object]:
    states = {"install": "staged", "upgrade": "activated", "rollback": "rolled-back"}
    return {
        "schema_version": OPERATION_RESULT_SCHEMA,
        "operation": command,
        "state": states[command],
        "component": component,
        "deployment_identity": IDENTITY,
        "state_identity": IDENTITY,
    }


def _doctor_result(*, healthy: bool) -> dict[str, object]:
    checks = (
        ("deployment-state", "deployment"),
        ("configuration", "host"),
        ("permissions", "host"),
        ("database", "database"),
        ("redis", "queue"),
        ("worker", "queue"),
        ("frontend", "platform"),
        ("encode-runtime", "runtime"),
        ("bulk-rnaseq-runtime", "runtime"),
        ("references", "reference"),
    )
    return {
        "schema_version": DOCTOR_RESULT_SCHEMA,
        "status": "healthy" if healthy else "unhealthy",
        "ready": healthy,
        "checks": [
            {
                "check_id": check_id,
                "category": category,
                "state": "pass" if healthy else "fail",
                "reason_code": "READY" if healthy else "DOCTOR_CHECK_FAILED",
            }
            for check_id, category in checks
        ],
    }


@dataclass
class RecordingBackend:
    result: PublicCommandResult | None = None
    called: tuple[object, ...] | None = None
    error: Exception | None = None

    def _call(self, *values: object) -> PublicCommandResult:
        self.called = values
        if self.error is not None:
            raise self.error
        if self.result is not None:
            return self.result
        command = values[0]
        if command == "status":
            return PublicCommandResult(_status_result())
        if command == "doctor":
            return PublicCommandResult(_doctor_result(healthy=True))
        if command == "verify":
            return PublicCommandResult(
                {
                    "schema_version": VERIFY_RESULT_SCHEMA,
                    "verified": True,
                    "deployment": _status_result(),
                    "frontend_identity": IDENTITY,
                    "database_schema_identity": None,
                }
            )
        return PublicCommandResult(_operation_result(command, values[1]))

    def install(self, component: str, bundle_path: Path) -> PublicCommandResult:
        return self._call("install", component, bundle_path)

    def status(self) -> PublicCommandResult:
        return self._call("status")

    def doctor(self) -> PublicCommandResult:
        return self._call("doctor")

    def verify(self) -> PublicCommandResult:
        return self._call("verify")

    def upgrade(self, component: str, bundle_path: Path) -> PublicCommandResult:
        return self._call("upgrade", component, bundle_path)

    def rollback(self, component: str, deployment_identity: str) -> PublicCommandResult:
        return self._call("rollback", component, deployment_identity)


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (("status",), CommandRequest("status")),
        (("doctor",), CommandRequest("doctor")),
        (("verify",), CommandRequest("verify")),
        (
            (
                "install",
                "--component",
                "platform",
                "--bundle",
                "/srv/releases/platform.tar",
            ),
            CommandRequest("install", "platform", Path("/srv/releases/platform.tar")),
        ),
        (
            (
                "upgrade",
                "--component",
                "encode-runtime",
                "--bundle",
                "/srv/releases/encode.tar",
            ),
            CommandRequest(
                "upgrade", "encode-runtime", Path("/srv/releases/encode.tar")
            ),
        ),
        (
            (
                "rollback",
                "--component",
                "bulk-rnaseq-runtime",
                "--identity",
                IDENTITY,
            ),
            CommandRequest(
                "rollback",
                "bulk-rnaseq-runtime",
                deployment_identity=IDENTITY,
            ),
        ),
    ],
)
def test_exact_command_grammar(argv: tuple[str, ...], expected: CommandRequest) -> None:
    assert parse_command(argv) == expected


@pytest.mark.parametrize(
    "argv",
    [
        (),
        ("install",),
        ("install", "--component", "platform", "--bundle", "relative.tar"),
        ("install", "--component", "platform", "--bundle", "/tmp/../secret"),
        ("install", "--bundle", "/tmp/a", "--component", "platform"),
        ("upgrade", "--component", "platform", "--bundle", "/tmp/a", "--also-runtime"),
        ("rollback", "--component", "unknown", "--identity", IDENTITY),
        (
            "rollback",
            "--component",
            "platform",
            "--identity",
            "sha256-not-an-identity",
        ),
        ("status", "--verbose"),
        ("doctor\ninstall",),
        ("exec", "/bin/sh", "-c", "id"),
    ],
)
def test_command_grammar_rejects_injection_and_ambiguous_scope(
    argv: tuple[str, ...],
) -> None:
    with pytest.raises(DeploymentError) as caught:
        parse_command(argv)

    assert caught.value.issue.code == "DEPLOYMENT_COMMAND_INVALID"
    assert "/tmp" not in str(caught.value)


def test_success_is_one_canonical_json_document() -> None:
    output = io.StringIO()
    errors = io.StringIO()
    backend = RecordingBackend()

    exit_code = main(
        ("status",),
        backend=backend,
        stdout=output,
        stderr=errors,
    )

    assert exit_code == EXIT_OK
    assert errors.getvalue() == ""
    assert output.getvalue().count("\n") == 1
    receipt = json.loads(output.getvalue())
    assert receipt["command"] == "status"
    assert receipt["result"] == _status_result()
    assert (
        output.getvalue()
        == json.dumps(
            receipt,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    )
    assert backend.called == ("status",)


def test_doctor_can_report_not_ready_with_a_stable_nonzero_exit() -> None:
    report = _doctor_result(healthy=False)
    backend = RecordingBackend(PublicCommandResult(report, exit_code=EXIT_UNAVAILABLE))
    output = io.StringIO()

    exit_code = main(("doctor",), backend=backend, stdout=output)

    assert exit_code == EXIT_UNAVAILABLE
    assert json.loads(output.getvalue())["status"] == "not-ready"
    assert json.loads(output.getvalue())["result"] == report


def test_doctor_result_rejects_unregistered_reason_codes() -> None:
    report = _doctor_result(healthy=True)
    report["checks"][0]["reason_code"] = "SYNTACTICALLY_VALID_BUT_UNREGISTERED"

    with pytest.raises(DeploymentError) as caught:
        PublicCommandResult(report)

    assert caught.value.issue.code == "DEPLOYMENT_RESULT_INVALID"


@pytest.mark.parametrize(
    ("command", "result"),
    [
        ("status", PublicCommandResult(_status_result(), EXIT_UNAVAILABLE)),
        ("doctor", PublicCommandResult(_doctor_result(healthy=True), EXIT_UNAVAILABLE)),
        (
            "verify",
            PublicCommandResult(
                {
                    "schema_version": VERIFY_RESULT_SCHEMA,
                    "verified": True,
                    "deployment": _status_result(),
                    "frontend_identity": IDENTITY,
                    "database_schema_identity": None,
                },
                EXIT_DATA,
            ),
        ),
    ],
)
def test_result_state_and_exit_code_must_agree(
    command: str,
    result: PublicCommandResult,
) -> None:
    errors = io.StringIO()

    exit_code = main(
        (command,),
        backend=RecordingBackend(result=result),
        stderr=errors,
    )

    assert exit_code == EXIT_OPERATION
    assert json.loads(errors.getvalue())["issue"]["code"] == (
        "DEPLOYMENT_RESULT_INVALID"
    )


@pytest.mark.parametrize(
    ("error", "expected_exit"),
    [
        (RuntimeError("/private/reference.fa token=secret"), EXIT_OPERATION),
        (
            fail(
                "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
                "/private/database is incompatible",
            ),
            EXIT_INCOMPATIBLE,
        ),
        (
            fail(
                "DEPLOYMENT_SCHEMA_OBSERVATION_DEFERRED",
                "/private/database observer unavailable",
                recoverable=True,
            ),
            EXIT_UNAVAILABLE,
        ),
    ],
)
def test_errors_are_redacted_and_have_stable_exit_classes(
    error: Exception,
    expected_exit: int,
) -> None:
    output = io.StringIO()
    errors = io.StringIO()
    backend = RecordingBackend(error=error)

    exit_code = main(
        ("verify",),
        backend=backend,
        stdout=output,
        stderr=errors,
    )

    assert exit_code == expected_exit
    assert output.getvalue() == ""
    receipt = json.loads(errors.getvalue())
    assert receipt["status"] == "error"
    assert "private" not in errors.getvalue()
    assert "secret" not in errors.getvalue()
    assert "reference.fa" not in errors.getvalue()


def test_backend_cannot_publish_paths_secrets_or_arbitrary_objects() -> None:
    with pytest.raises(DeploymentError) as path_error:
        PublicCommandResult({"private_path": "/private/reference.fa"})
    assert path_error.value.issue.code == "DEPLOYMENT_RESULT_INVALID"

    with pytest.raises(DeploymentError):
        PublicCommandResult(
            {"schema_version": STATUS_RESULT_SCHEMA, "result": object()}
        )

    for unsafe in (
        {**_status_result(), "message": "permission denied"},
        {**_status_result(), "credential": "abc"},
    ):
        with pytest.raises(DeploymentError) as caught:
            PublicCommandResult(unsafe)
        assert caught.value.issue.code == "DEPLOYMENT_RESULT_INVALID"


def test_phase_a_default_backend_fails_closed_without_host_action() -> None:
    errors = io.StringIO()

    exit_code = main(("status",), stderr=errors)

    assert exit_code == EXIT_UNAVAILABLE
    assert json.loads(errors.getvalue())["issue"]["code"] == (
        "DEPLOYMENT_INTEGRATION_DEFERRED"
    )


def test_invalid_usage_never_echoes_the_rejected_argument() -> None:
    errors = io.StringIO()

    exit_code = main(("install", "--bundle", "/private/secret"), stderr=errors)

    assert exit_code == EXIT_USAGE
    assert "/private" not in errors.getvalue()
    assert "secret" not in errors.getvalue()
