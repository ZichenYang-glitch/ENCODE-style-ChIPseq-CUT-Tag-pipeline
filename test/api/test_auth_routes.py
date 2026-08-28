"""LAN authentication route and fail-closed authorization boundary tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

fastapi = pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import TEST_AUTH_USER_ID, seed_test_authentication  # noqa: E402
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.services.authentication import digest_session_token  # noqa: E402


@pytest.fixture
def app(tmp_path):
    application = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    yield application
    application.state.persistence.close()


@pytest.fixture
def seeded_app(app):
    seed_test_authentication(app)
    return app


def _login(
    client: ApiTestClient,
    username: str = "root-admin",
    password: str = "correct horse battery staple",
) -> dict[str, str]:
    response = client.post(
        "/api/v1/auth/login",
        json={"username": username, "password": password},
    )
    assert response.status_code == 200
    cookies = {}
    for header in response.headers.get_list("set-cookie"):
        name, rest = header.split("=", 1)
        cookies[name] = rest.split(";", 1)[0]
    return cookies


def _bootstrap_admin(app, username: str = "root-admin") -> None:
    password_manager = app.state.account_administration_service._password_manager
    from datetime import datetime, timezone

    from encode_pipeline.platform.authentication import (
        UserAccount,
        UserRole,
        UserStatus,
    )
    from encode_pipeline.services.authentication import new_user_id

    now = datetime.now(timezone.utc)
    account = UserAccount(
        user_id=new_user_id(),
        username=username,
        role=UserRole.ADMINISTRATOR,
        status=UserStatus.ENABLED,
        password_hash=password_manager.hash_password("correct horse battery staple"),
        created_at=now,
        updated_at=now,
        password_changed_at=now,
    )
    app.state.authentication_repository.create_account(account)


def _cookie_header(cookies: dict[str, str]) -> str:
    return "; ".join(f"{name}={value}" for name, value in cookies.items())


def _member_session_headers(app, username: str) -> dict[str, str]:
    from datetime import datetime, timedelta, timezone

    from encode_pipeline.services.authentication import (
        new_session_record,
        new_session_secrets,
    )

    member = app.state.authentication_repository.get_account_by_username(username)
    secrets = new_session_secrets()
    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=member.user_id,
            secrets=secrets,
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(hours=1),
        )
    )
    policy = app.state.auth_cookie_policy
    return {
        "Cookie": (
            f"{policy.session_cookie.name}={secrets.session_token}; "
            f"{policy.csrf_cookie.name}={secrets.csrf_token}"
        ),
        "X-CSRF-Token": secrets.csrf_token,
    }


def test_setup_required_until_the_first_administrator_exists(app) -> None:
    with ApiTestClient(app) as client:
        state = client.get("/api/v1/auth/session")
        assert state.status_code == 200
        assert state.json()["setup_required"] is True
        assert state.json()["authenticated"] is False

        login = client.post(
            "/api/v1/auth/login",
            json={"username": "root-admin", "password": "x" * 15},
        )
        assert login.status_code == 503
        assert login.json()["issues"][0]["code"] == "SETUP_REQUIRED"

        protected = client.get("/api/v1/workflows")
        assert protected.status_code == 503
        assert protected.json()["issues"][0]["code"] == "SETUP_REQUIRED"


def test_anonymous_requests_fail_closed_on_business_and_admin_routes(
    seeded_app,
) -> None:
    with ApiTestClient(seeded_app) as client:
        client.app.state.test_auth_tokens = None
        for method, url in (
            ("get", "/api/v1/workflows"),
            ("get", "/api/v1/runs"),
            ("get", "/api/v1/auth/accounts"),
            ("get", "/api/v1/auth/preferences/terminal-email"),
        ):
            response = getattr(client, method)(url)
            assert response.status_code == 401, url
            assert response.json()["issues"][0]["code"] == "AUTHENTICATION_REQUIRED"


def test_login_issues_cookies_and_rejects_bad_credentials(app) -> None:
    _bootstrap_admin(app)
    with ApiTestClient(app) as client:
        bad = client.post(
            "/api/v1/auth/login",
            json={"username": "root-admin", "password": "wrong password here"},
        )
        assert bad.status_code == 401
        assert bad.json()["issues"][0]["code"] == "INVALID_CREDENTIALS"
        assert "root-admin" not in bad.text

        cookies = _login(client)
        assert "helixweave_session" in cookies
        assert "helixweave_csrf" in cookies

        policy = app.state.auth_cookie_policy
        session_token = cookies[policy.session_cookie.name]
        csrf_token = cookies[policy.csrf_cookie.name]
        headers = {
            "Cookie": _cookie_header(cookies),
            "X-CSRF-Token": csrf_token,
        }
        state = client.get("/api/v1/auth/session", headers=headers)
        assert state.status_code == 200
        body = state.json()
        assert body["authenticated"] is True
        assert body["principal"]["username"] == "root-admin"

        record = app.state.authentication_repository.get_session(
            digest_session_token(session_token)
        )
        assert record.csrf_digest is not None


def test_mutations_require_csrf_and_logout_revokes(app) -> None:
    _bootstrap_admin(app)
    with ApiTestClient(app) as client:
        cookies = _login(client)
        cookie_header = _cookie_header(cookies)

        missing_csrf = client.post(
            "/api/v1/auth/accounts",
            headers={"Cookie": cookie_header},
            json={"username": "lab-member", "password": "member password phrase"},
        )
        assert missing_csrf.status_code == 403
        assert missing_csrf.json()["issues"][0]["code"] == "CSRF_INVALID"

        headers = {"Cookie": cookie_header, "X-CSRF-Token": cookies["helixweave_csrf"]}
        created = client.post(
            "/api/v1/auth/accounts",
            headers=headers,
            json={"username": "lab-member", "password": "member password phrase"},
        )
        assert created.status_code == 200
        assert created.json()["account"]["role"] == "member"

        policy = app.state.auth_cookie_policy
        session_token = cookies[policy.session_cookie.name]
        logout = client.post("/api/v1/auth/logout", headers=headers)
        assert logout.status_code == 200
        record = app.state.authentication_repository.get_session(
            digest_session_token(session_token)
        )
        assert record.revoked_at is not None

        after = client.get("/api/v1/workflows", headers=headers)
        assert after.status_code == 401


def test_login_same_origin_json_contract(app) -> None:
    """v1 supports a same-origin JSON login only.

    A same-origin JSON login creates a session; form-urlencoded and text/plain
    submissions cannot create one; responses carry no CORS reflection.
    """
    _bootstrap_admin(app)
    policy = app.state.auth_cookie_policy
    session_cookie_prefix = f"{policy.session_cookie.name}="
    with ApiTestClient(app) as client:
        same_origin = client.post(
            "/api/v1/auth/login",
            json={
                "username": "root-admin",
                "password": "correct horse battery staple",
            },
            headers={"Origin": "http://testserver"},
        )
        assert same_origin.status_code == 200
        assert same_origin.json()["ok"] is True
        assert any(
            header.startswith(session_cookie_prefix)
            for header in same_origin.headers.get_list("set-cookie")
        )
        assert "access-control-allow-origin" not in same_origin.headers

        for kwargs in (
            {
                "data": {
                    "username": "root-admin",
                    "password": "correct horse battery staple",
                }
            },
            {
                "content": (
                    '{"username": "root-admin", '
                    '"password": "correct horse battery staple"}'
                ),
                "headers": {"Content-Type": "text/plain"},
            },
        ):
            rejected = client.post("/api/v1/auth/login", **kwargs)
            assert rejected.status_code == 400
            assert not any(
                header.startswith(session_cookie_prefix)
                for header in rejected.headers.get_list("set-cookie")
            )


def test_cross_origin_preflight_gets_no_permissive_cors(app) -> None:
    """CORS is not configured: a cross-origin preflight must not receive a
    permissive Access-Control-Allow-Origin."""
    _bootstrap_admin(app)
    with ApiTestClient(app) as client:
        preflight = client.request(
            "OPTIONS",
            "/api/v1/auth/login",
            headers={
                "Origin": "https://attacker.example",
                "Access-Control-Request-Method": "POST",
                "Access-Control-Request-Headers": "content-type",
            },
        )
        allow_origin = preflight.headers.get("access-control-allow-origin")
        assert allow_origin is None or allow_origin not in (
            "*",
            "https://attacker.example",
        )


def test_member_cannot_reach_admin_account_routes(app) -> None:
    _bootstrap_admin(app)
    with ApiTestClient(app) as client:
        cookies = _login(client)
        headers = {
            "Cookie": _cookie_header(cookies),
            "X-CSRF-Token": cookies["helixweave_csrf"],
        }
        client.post(
            "/api/v1/auth/accounts",
            headers=headers,
            json={"username": "lab-member", "password": "member password phrase"},
        )
    member = app.state.authentication_repository.get_account_by_username("lab-member")
    member_headers = _member_session_headers(app, "lab-member")
    with ApiTestClient(app) as client:
        allowed = client.get("/api/v1/workflows/", headers=member_headers)
        assert allowed.status_code == 200
        for method, url in (
            ("get", "/api/v1/auth/accounts"),
            ("post", "/api/v1/auth/accounts"),
            ("post", f"/api/v1/auth/accounts/{member.user_id}/status"),
            ("post", f"/api/v1/auth/accounts/{member.user_id}/password"),
            ("post", f"/api/v1/auth/accounts/{member.user_id}/sessions/revoke"),
        ):
            kwargs = {"headers": member_headers}
            if method == "post":
                kwargs["json"] = {"enabled": False, "password": "x" * 15}
            response = getattr(client, method)(url, **kwargs)
            assert response.status_code == 403, url
            assert response.json()["issues"][0]["code"] == "ADMINISTRATOR_REQUIRED"


def test_member_terminal_email_preference_is_address_free_and_csrf_protected(
    app,
) -> None:
    _bootstrap_admin(app)
    with ApiTestClient(app) as client:
        admin_cookies = _login(client)
        admin_headers = {
            "Cookie": _cookie_header(admin_cookies),
            "X-CSRF-Token": admin_cookies["helixweave_csrf"],
        }
        created = client.post(
            "/api/v1/auth/accounts",
            headers=admin_headers,
            json={"username": "lab-member", "password": "member password phrase"},
        )
        assert created.status_code == 200

        admin_preference = client.get(
            "/api/v1/auth/preferences/terminal-email",
            headers=admin_headers,
        )
        assert admin_preference.status_code == 403
        assert admin_preference.json()["issues"][0]["code"] == "MEMBER_REQUIRED"

    app.state.account_administration_service.set_notification_email_for_username(
        "lab-member",
        "member@example.test",
    )
    member_headers = _member_session_headers(app, "lab-member")
    cookie_only = {"Cookie": member_headers["Cookie"]}
    with ApiTestClient(app) as client:
        current = client.get(
            "/api/v1/auth/preferences/terminal-email",
            headers=cookie_only,
        )
        assert current.status_code == 200
        assert current.json() == {
            "terminal_email_enabled": True,
            "address_configured": True,
        }
        assert "member@example.test" not in current.text

        rejected = client.request(
            "PATCH",
            "/api/v1/auth/preferences/terminal-email",
            headers=cookie_only,
            json={"terminal_email_enabled": False},
        )
        assert rejected.status_code == 403
        assert rejected.json()["issues"][0]["code"] == "CSRF_INVALID"

        updated = client.request(
            "PATCH",
            "/api/v1/auth/preferences/terminal-email",
            headers=member_headers,
            json={"terminal_email_enabled": False},
        )
        assert updated.status_code == 200
        assert updated.json() == {
            "terminal_email_enabled": False,
            "address_configured": True,
        }

        extra = client.request(
            "PATCH",
            "/api/v1/auth/preferences/terminal-email",
            headers=member_headers,
            json={
                "terminal_email_enabled": True,
                "notification_email": "attacker@example.test",
            },
        )
        assert extra.status_code == 400
        assert "attacker@example.test" not in extra.text


def test_admin_account_lifecycle_routes(app) -> None:
    _bootstrap_admin(app)
    with ApiTestClient(app) as client:
        cookies = _login(client)
        headers = {
            "Cookie": _cookie_header(cookies),
            "X-CSRF-Token": cookies["helixweave_csrf"],
        }

        listed = client.get("/api/v1/auth/accounts", headers=headers)
        assert listed.status_code == 200
        assert [a["username"] for a in listed.json()["accounts"]] == ["root-admin"]

        created = client.post(
            "/api/v1/auth/accounts",
            headers=headers,
            json={"username": "lab-member", "password": "member password phrase"},
        )
        member_id = created.json()["account"]["user_id"]

        conflict = client.post(
            "/api/v1/auth/accounts",
            headers=headers,
            json={"username": "lab-member", "password": "member password phrase"},
        )
        assert conflict.status_code == 409

        disabled = client.post(
            f"/api/v1/auth/accounts/{member_id}/status",
            headers=headers,
            json={"enabled": False},
        )
        assert disabled.status_code == 200
        assert disabled.json()["account"]["status"] == "disabled"

        reset = client.post(
            f"/api/v1/auth/accounts/{member_id}/password",
            headers=headers,
            json={"password": "a new member password phrase"},
        )
        assert reset.status_code == 200

        revoked = client.post(
            f"/api/v1/auth/accounts/{member_id}/sessions/revoke",
            headers=headers,
        )
        assert revoked.status_code == 200
        assert revoked.json()["revoked_count"] == 0

        missing = client.post(
            "/api/v1/auth/accounts/usr_ffffffffffffffffffffffffffffffff/status",
            headers=headers,
            json={"enabled": False},
        )
        assert missing.status_code == 404


def test_seeded_member_session_works_with_test_client(seeded_app) -> None:
    with ApiTestClient(seeded_app) as client:
        response = client.get("/api/v1/workflows")
        assert response.status_code == 200
        accounts = client.get("/api/v1/auth/accounts")
        assert accounts.status_code == 200
        assert any(
            account["user_id"] == TEST_AUTH_USER_ID
            for account in accounts.json()["accounts"]
        )


def test_artifact_download_is_audited_with_actor_and_target(seeded_app) -> None:
    from types import SimpleNamespace

    from encode_pipeline.platform.results import Result
    from encode_pipeline.platform.security_audit import AuditAction, AuditResourceKind

    plan = SimpleNamespace(
        iter_bytes=lambda: iter((b"artifact-bytes",)),
        close=lambda: None,
        media_type="application/octet-stream",
        content_disposition='attachment; filename="artifact.bin"',
        size_bytes=15,
    )

    class StubDownloadService:
        def prepare(self, run_id, artifact_id, **kwargs):
            return Result.success(plan)

    seeded_app.state.artifact_download_service = StubDownloadService()
    with ApiTestClient(seeded_app) as client:
        response = client.get(
            "/api/v1/runs/run-1/artifacts/artifact-1/download"
            "?generation="
            + "artifactgen-"
            + "a" * 64
            + "&revision="
            + "artifactrev-"
            + "b" * 64
        )
        assert response.status_code == 200
        assert response.content == b"artifact-bytes"

    (event,) = seeded_app.state.authentication_repository.list_security_audit_events()
    assert event.action is AuditAction.ARTIFACT_DOWNLOAD
    assert event.actor_user_id == TEST_AUTH_USER_ID
    assert event.resource is not None
    assert event.resource.kind is AuditResourceKind.ARTIFACT
    assert "run-1" not in event.resource.resource_id
