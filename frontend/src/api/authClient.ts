import { ApiError } from './fetcher';

export interface AuthPrincipal {
  user_id: string;
  username: string;
  role: string;
}

export interface SessionState {
  ok: boolean;
  setup_required: boolean;
  authenticated: boolean;
  principal: AuthPrincipal | null;
}

const CSRF_COOKIE_NAME = 'helixweave_csrf';

function readCookie(name: string): string | null {
  if (typeof document === 'undefined') {
    return null;
  }
  const prefix = `${name}=`;
  for (const part of document.cookie.split(';')) {
    const trimmed = part.trim();
    if (trimmed.startsWith(prefix)) {
      return decodeURIComponent(trimmed.slice(prefix.length));
    }
  }
  return null;
}

export function csrfHeader(): Record<string, string> {
  const token = readCookie(CSRF_COOKIE_NAME);
  return token === null ? {} : { 'X-CSRF-Token': token };
}

async function parseEnvelope(response: Response): Promise<unknown> {
  try {
    return await response.json();
  } catch {
    return null;
  }
}

async function parseSessionState(response: Response): Promise<SessionState> {
  const body = (await parseEnvelope(response)) as SessionState | null;
  if (body == null || typeof body.ok !== 'boolean') {
    throw new ApiError(response.status, 'AUTH_SESSION_INVALID', 'Session state is invalid.');
  }
  return body;
}

export async function fetchSessionState(): Promise<SessionState> {
  const response = await fetch('/api/v1/auth/session', {
    method: 'GET',
    credentials: 'same-origin',
  });
  return parseSessionState(response);
}

export async function login(
  username: string,
  password: string,
): Promise<AuthPrincipal> {
  const response = await fetch('/api/v1/auth/login', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ username, password }),
    credentials: 'same-origin',
  });
  if (!response.ok) {
    const body = (await parseEnvelope(response)) as {
      issues?: { code?: string; message?: string }[];
    } | null;
    const issue = body?.issues?.[0];
    throw new ApiError(
      response.status,
      issue?.code ?? 'LOGIN_FAILED',
      issue?.message ?? 'Login failed.',
    );
  }
  const body = (await parseEnvelope(response)) as {
    ok: boolean;
    principal: AuthPrincipal | null;
  } | null;
  if (body == null || !body.ok || body.principal == null) {
    throw new ApiError(response.status, 'LOGIN_FAILED', 'Login failed.');
  }
  return body.principal;
}

export async function logout(): Promise<void> {
  const response = await fetch('/api/v1/auth/logout', {
    method: 'POST',
    headers: { ...csrfHeader() },
    credentials: 'same-origin',
  });
  if (!response.ok) {
    throw new ApiError(response.status, 'LOGOUT_FAILED', 'Logout failed.');
  }
}
