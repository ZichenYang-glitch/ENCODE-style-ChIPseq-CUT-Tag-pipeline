"""LAN authentication session and account administration API routes."""

from __future__ import annotations

from fastapi import APIRouter, Depends, Request, Response
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import (
    AuthApiError,
    enforce_csrf,
    get_account_administration_service,
    get_auth_cookie_policy,
    get_authentication_service,
    require_administrator,
    require_principal,
)
from encode_pipeline.api.models import (
    AccountCreateRequest,
    AccountListResponse,
    AccountMutationResponse,
    AccountPasswordResetRequest,
    AccountSessionsRevokeResponse,
    AccountStatusRequest,
    AccountSummaryResponse,
    AuthErrorResponse,
    IssueResponse,
    LoginRequest,
    LoginResponse,
    PrincipalResponse,
    SessionStateResponse,
    TerminalEmailPreferenceRequest,
    TerminalEmailPreferenceResponse,
)
from encode_pipeline.platform.authentication import (
    AuthenticatedPrincipal,
    UserAccount,
    UserRole,
    validate_user_id,
)
from encode_pipeline.services.authentication import BrowserSessionCookiePolicy
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
    AuthenticationActor,
    AuthenticationError,
    AuthenticationService,
)

router = APIRouter(prefix="/auth", tags=["auth"])

_ADMIN_DEPENDENCIES = [Depends(require_administrator), Depends(enforce_csrf)]


def _issue(code: str, message: str) -> IssueResponse:
    return IssueResponse(code=code, message=message, source="api")


def _auth_error(code: str, message: str, status_code: int) -> JSONResponse:
    body = AuthErrorResponse(issues=[_issue(code, message)])
    return JSONResponse(status_code=status_code, content=body.model_dump(mode="json"))


def _client_identity(request: Request) -> str:
    client = request.client
    if client is None or not client.host:
        return "unknown"
    return client.host


def _principal_response(principal: AuthenticatedPrincipal) -> PrincipalResponse:
    return PrincipalResponse(
        user_id=principal.user_id,
        username=principal.username,
        role=principal.role.value,
    )


def _account_summary(account: UserAccount) -> AccountSummaryResponse:
    summary = account.to_public_summary()
    return AccountSummaryResponse(**summary)


def _require_member(principal: AuthenticatedPrincipal) -> None:
    if principal.role is not UserRole.MEMBER:
        raise AuthApiError(
            403,
            "MEMBER_REQUIRED",
            "Terminal email preference is available to members.",
        )


def _service_error_body(error: AuthenticationError) -> JSONResponse:
    status = {
        "INVALID_CREDENTIALS": 401,
        "LOGIN_RATE_LIMITED": 401,
        "ADMINISTRATOR_REQUIRED": 403,
        "OPERATION_CONFLICT": 409,
        "RESOURCE_NOT_FOUND": 404,
    }[error.reason_code]
    body = AuthErrorResponse(issues=[_issue(error.reason_code, str(error))])
    return JSONResponse(status_code=status, content=body.model_dump(mode="json"))


@router.post(
    "/login",
    operation_id="login",
    response_model=LoginResponse,
    responses={401: {"model": AuthErrorResponse}},
)
def login(
    payload: LoginRequest,
    request: Request,
    response: Response,
    authentication: AuthenticationService = Depends(get_authentication_service),
    policy: BrowserSessionCookiePolicy = Depends(get_auth_cookie_policy),
) -> LoginResponse | JSONResponse:
    """Authenticate one bounded login and issue browser session cookies."""
    if not authentication.setup_complete():
        return _auth_error(
            "SETUP_REQUIRED",
            "The deployment requires initial administrator setup.",
            503,
        )
    try:
        session = authentication.login(
            payload.username,
            payload.password,
            client_identity=_client_identity(request),
        )
    except AuthenticationError:
        return _auth_error(
            "INVALID_CREDENTIALS",
            "The username or password is invalid.",
            401,
        )
    response.set_cookie(
        value=session.secrets.session_token,
        **policy.session_cookie.to_response_kwargs(),
    )
    response.set_cookie(
        value=session.secrets.csrf_token,
        **policy.csrf_cookie.to_response_kwargs(),
    )
    return LoginResponse(
        ok=True,
        principal=_principal_response(session.principal),
        issues=[],
    )


@router.post(
    "/logout",
    operation_id="logout",
    response_model=SessionStateResponse,
    dependencies=[Depends(enforce_csrf)],
)
def logout(
    request: Request,
    response: Response,
    authentication: AuthenticationService = Depends(get_authentication_service),
    policy: BrowserSessionCookiePolicy = Depends(get_auth_cookie_policy),
) -> SessionStateResponse:
    """Revoke the current session and clear both cookies."""
    token = request.cookies.get(policy.session_cookie.name)
    if token is not None:
        authentication.logout(token)
    for directive in (policy.session_cookie, policy.csrf_cookie):
        response.delete_cookie(key=directive.name, path=directive.path)
    return SessionStateResponse(
        ok=True,
        setup_required=not authentication.setup_complete(),
        authenticated=False,
        principal=None,
        issues=[],
    )


@router.get(
    "/session", response_model=SessionStateResponse, operation_id="session_state"
)
def session_state(
    request: Request,
    authentication: AuthenticationService = Depends(get_authentication_service),
    policy: BrowserSessionCookiePolicy = Depends(get_auth_cookie_policy),
) -> SessionStateResponse:
    """Report non-sensitive setup and current-session state."""
    setup_complete = authentication.setup_complete()
    principal = None
    token = request.cookies.get(policy.session_cookie.name)
    if token is not None:
        principal = authentication.resolve_session(token)
    return SessionStateResponse(
        ok=True,
        setup_required=not setup_complete,
        authenticated=principal is not None,
        principal=None if principal is None else _principal_response(principal),
        issues=[],
    )


@router.get(
    "/preferences/terminal-email",
    operation_id="get_terminal_email_preference",
    response_model=TerminalEmailPreferenceResponse,
    responses={401: {"model": AuthErrorResponse}, 403: {"model": AuthErrorResponse}},
)
def get_terminal_email_preference(
    principal: AuthenticatedPrincipal = Depends(require_principal),
    authentication: AuthenticationService = Depends(get_authentication_service),
) -> TerminalEmailPreferenceResponse | JSONResponse:
    """Return the member's address-free terminal-email preference."""
    _require_member(principal)
    try:
        preference = authentication.get_terminal_email_preference(principal)
    except AuthenticationError as error:
        return _service_error_body(error)
    return TerminalEmailPreferenceResponse(
        terminal_email_enabled=preference.terminal_email_enabled,
        address_configured=preference.address_configured,
    )


@router.patch(
    "/preferences/terminal-email",
    operation_id="set_terminal_email_preference",
    response_model=TerminalEmailPreferenceResponse,
    dependencies=[Depends(enforce_csrf)],
    responses={401: {"model": AuthErrorResponse}, 403: {"model": AuthErrorResponse}},
)
def set_terminal_email_preference(
    payload: TerminalEmailPreferenceRequest,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    authentication: AuthenticationService = Depends(get_authentication_service),
) -> TerminalEmailPreferenceResponse | JSONResponse:
    """Update only the member's terminal-email opt-out flag."""
    _require_member(principal)
    try:
        preference = authentication.set_terminal_email_enabled(
            principal,
            payload.terminal_email_enabled,
        )
    except AuthenticationError as error:
        return _service_error_body(error)
    return TerminalEmailPreferenceResponse(
        terminal_email_enabled=preference.terminal_email_enabled,
        address_configured=preference.address_configured,
    )


@router.get(
    "/accounts",
    operation_id="list_accounts",
    response_model=AccountListResponse,
    dependencies=_ADMIN_DEPENDENCIES,
)
def list_accounts(
    administration: AccountAdministrationService = Depends(
        get_account_administration_service
    ),
) -> AccountListResponse:
    """List safe account summaries for administrators."""
    accounts = [_account_summary(account) for account in administration.list_accounts()]
    return AccountListResponse(ok=True, accounts=accounts, issues=[])


@router.post(
    "/accounts",
    operation_id="create_member_account",
    response_model=AccountMutationResponse,
    dependencies=_ADMIN_DEPENDENCIES,
)
def create_member_account(
    payload: AccountCreateRequest,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    administration: AccountAdministrationService = Depends(
        get_account_administration_service
    ),
) -> AccountMutationResponse | JSONResponse:
    """Create one member account; administrators only."""
    try:
        account = administration.create_member(
            AuthenticationActor.for_principal(principal),
            payload.username,
            payload.password,
        )
    except AuthenticationError as error:
        return _service_error_body(error)
    return AccountMutationResponse(
        ok=True,
        account=_account_summary(account),
        issues=[],
    )


@router.post(
    "/accounts/{user_id}/status",
    operation_id="set_account_status",
    response_model=AccountMutationResponse,
    dependencies=_ADMIN_DEPENDENCIES,
)
def set_account_status(
    user_id: str,
    payload: AccountStatusRequest,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    administration: AccountAdministrationService = Depends(
        get_account_administration_service
    ),
) -> AccountMutationResponse | JSONResponse:
    """Enable or disable one account; administrators only."""
    try:
        validate_user_id(user_id)
        account = administration.set_account_status(
            AuthenticationActor.for_principal(principal),
            user_id,
            enabled=payload.enabled,
        )
    except ValueError:
        return _auth_error("RESOURCE_NOT_FOUND", "The account does not exist.", 404)
    except AuthenticationError as error:
        return _service_error_body(error)
    return AccountMutationResponse(
        ok=True,
        account=_account_summary(account),
        issues=[],
    )


@router.post(
    "/accounts/{user_id}/password",
    operation_id="reset_account_password",
    response_model=AccountMutationResponse,
    dependencies=_ADMIN_DEPENDENCIES,
)
def reset_account_password(
    user_id: str,
    payload: AccountPasswordResetRequest,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    administration: AccountAdministrationService = Depends(
        get_account_administration_service
    ),
) -> AccountMutationResponse | JSONResponse:
    """Reset one account password and revoke its sessions; administrators only."""
    try:
        validate_user_id(user_id)
        account = administration.reset_password(
            AuthenticationActor.for_principal(principal),
            user_id,
            payload.password,
        )
    except ValueError:
        return _auth_error("RESOURCE_NOT_FOUND", "The account does not exist.", 404)
    except AuthenticationError as error:
        return _service_error_body(error)
    return AccountMutationResponse(
        ok=True,
        account=_account_summary(account),
        issues=[],
    )


@router.post(
    "/accounts/{user_id}/sessions/revoke",
    operation_id="revoke_account_sessions",
    response_model=AccountSessionsRevokeResponse,
    dependencies=_ADMIN_DEPENDENCIES,
)
def revoke_account_sessions(
    user_id: str,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    administration: AccountAdministrationService = Depends(
        get_account_administration_service
    ),
) -> AccountSessionsRevokeResponse | JSONResponse:
    """Revoke every session of one account; administrators only."""
    try:
        validate_user_id(user_id)
        revoked = administration.revoke_sessions(
            AuthenticationActor.for_principal(principal),
            user_id,
        )
    except ValueError:
        return _auth_error("RESOURCE_NOT_FOUND", "The account does not exist.", 404)
    except AuthenticationError as error:
        return _service_error_body(error)
    return AccountSessionsRevokeResponse(ok=True, revoked_count=revoked, issues=[])
