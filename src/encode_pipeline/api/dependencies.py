"""FastAPI dependencies shared by API routes."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastapi import Depends, Request

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.authentication import (
    CSRF_HEADER_NAME,
    request_requires_csrf,
)
from encode_pipeline.services.validation import ValidationService

if TYPE_CHECKING:
    from encode_pipeline.platform.authentication import AuthenticatedPrincipal
    from encode_pipeline.services.authentication import BrowserSessionCookiePolicy
    from encode_pipeline.services.authentication_service import (
        AccountAdministrationService,
        AuthenticationService,
    )
    from encode_pipeline.services.agent import AgentService
    from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
    from encode_pipeline.services.artifact_publications import (
        ArtifactPublicationQueryService,
    )
    from encode_pipeline.services.preflight import LocalPreflightService
    from encode_pipeline.services.run_submission import RunSubmissionService
    from encode_pipeline.services.run_cancellation import RunCancellationService
    from encode_pipeline.services.runs import RunService
    from encode_pipeline.services.reference_profiles import ReferenceProfileService
    from encode_pipeline.services.validated_inputs import (
        ValidatedInputService,
        ValidatedRunCreationService,
    )


async def get_registry(request: Request) -> WorkflowRegistry:
    """Return the app registry."""
    return request.app.state.registry


async def get_validation_service(request: Request) -> ValidationService:
    """Return the app validation service."""
    return request.app.state.validation_service


async def get_validated_input_service(request: Request) -> "ValidatedInputService":
    """Return successful-validation snapshot orchestration."""
    return request.app.state.validated_input_service


async def get_reference_profile_service(
    request: Request,
) -> "ReferenceProfileService":
    """Return the workflow-neutral Reference Profile catalog and resolver."""
    return request.app.state.reference_profile_service


async def get_validated_run_creation_service(
    request: Request,
) -> "ValidatedRunCreationService":
    """Return snapshot-only run creation orchestration."""
    return request.app.state.validated_run_creation_service


async def get_agent_service(request: Request) -> "AgentService":
    """Return the app agent service."""
    return request.app.state.agent_service


async def get_run_service(request: Request) -> "RunService":
    """Return the app run service."""
    return request.app.state.run_service


async def get_artifact_download_service(
    request: Request,
) -> "ArtifactDownloadService":
    """Return the descriptor-safe artifact download service."""
    return request.app.state.artifact_download_service


async def get_artifact_publication_service(
    request: Request,
) -> "ArtifactPublicationQueryService":
    """Return the disclosure-safe publication query service."""
    return request.app.state.artifact_publication_service


async def get_run_submission_service(request: Request) -> "RunSubmissionService":
    """Return the app durable run submission service."""
    return request.app.state.run_submission_service


async def get_run_cancellation_service(request: Request) -> "RunCancellationService":
    """Return the app durable cancellation orchestration service."""
    return request.app.state.run_cancellation_service


async def get_preflight_service(request: Request) -> "LocalPreflightService":
    """Return the app local preflight service."""
    return request.app.state.preflight_service


class AuthApiError(Exception):
    """Authentication or authorization failure mapped to a stable envelope."""

    def __init__(self, status_code: int, code: str, message: str) -> None:
        self.status_code = status_code
        self.code = code
        self.message = message
        super().__init__(message)


async def get_authentication_service(request: Request) -> "AuthenticationService":
    """Return the login/session lifecycle service."""
    return request.app.state.authentication_service


async def get_account_administration_service(
    request: Request,
) -> "AccountAdministrationService":
    """Return the administrator/local-operator account service."""
    return request.app.state.account_administration_service


async def get_auth_cookie_policy(
    request: Request,
) -> "BrowserSessionCookiePolicy":
    """Return the deployment-configured browser cookie policy."""
    return request.app.state.auth_cookie_policy


def _session_token(
    request: Request,
    policy: "BrowserSessionCookiePolicy",
) -> str | None:
    return request.cookies.get(policy.session_cookie.name)


async def get_optional_principal(
    request: Request,
    authentication: "AuthenticationService" = Depends(get_authentication_service),
    policy: "BrowserSessionCookiePolicy" = Depends(get_auth_cookie_policy),
) -> "AuthenticatedPrincipal | None":
    """Resolve the current browser session, or None when unauthenticated."""
    token = _session_token(request, policy)
    if token is None:
        return None
    return authentication.resolve_session(token)


async def require_principal(
    request: Request,
    authentication: "AuthenticationService" = Depends(get_authentication_service),
    policy: "BrowserSessionCookiePolicy" = Depends(get_auth_cookie_policy),
) -> "AuthenticatedPrincipal":
    """Fail closed unless the request carries one valid enabled-member session."""
    if not authentication.setup_complete():
        raise AuthApiError(
            503,
            "SETUP_REQUIRED",
            "The deployment requires initial administrator setup.",
        )
    principal = await get_optional_principal(request, authentication, policy)
    if principal is None:
        raise AuthApiError(
            401,
            "AUTHENTICATION_REQUIRED",
            "Authentication is required.",
        )
    return principal


async def require_administrator(
    principal: "AuthenticatedPrincipal" = Depends(require_principal),
) -> "AuthenticatedPrincipal":
    """Require the resolved principal to hold the administrator role."""
    if not principal.is_administrator:
        raise AuthApiError(
            403,
            "ADMINISTRATOR_REQUIRED",
            "The operation requires an administrator.",
        )
    return principal


async def enforce_csrf(
    request: Request,
    authentication: "AuthenticationService" = Depends(get_authentication_service),
    policy: "BrowserSessionCookiePolicy" = Depends(get_auth_cookie_policy),
) -> None:
    """Require cookie, header, and durable digest to agree on unsafe methods."""
    if not request_requires_csrf(request.method):
        return
    session_token = _session_token(request, policy)
    cookie_token = request.cookies.get(policy.csrf_cookie.name)
    header_token = request.headers.get(CSRF_HEADER_NAME)
    if session_token is None or not authentication.session_csrf_valid(
        session_token,
        cookie_token,
        header_token,
    ):
        raise AuthApiError(403, "CSRF_INVALID", "The CSRF token is invalid.")
