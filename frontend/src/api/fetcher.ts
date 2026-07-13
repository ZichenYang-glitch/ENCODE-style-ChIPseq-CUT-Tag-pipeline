export interface PublicApiIssue {
  code: string;
  message: string;
  severity?: 'error' | 'warning' | 'info';
  path?: string | null;
  source?: string | null;
  hint?: string | null;
}

export class ApiError extends Error {
  status: number;
  code: string;
  issues: PublicApiIssue[];

  constructor(
    status: number,
    code: string,
    message: string,
    issues: PublicApiIssue[] = [],
  ) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.code = code;
    this.issues = issues;
  }
}

function getApiBaseUrl(): string {
  const envBase = import.meta.env?.VITE_API_BASE_URL;
  if (typeof envBase !== 'string' || envBase.trim() === '') {
    return '';
  }
  return envBase.replace(/\/$/, '');
}

function publicIssue(value: unknown): PublicApiIssue | null {
  if (value == null || typeof value !== 'object') return null;
  const issue = value as Record<string, unknown>;
  if (typeof issue.code !== 'string' || typeof issue.message !== 'string') {
    return null;
  }

  const result: PublicApiIssue = {
    code: issue.code,
    message: issue.message,
  };
  if (
    issue.severity === 'error' ||
    issue.severity === 'warning' ||
    issue.severity === 'info'
  ) {
    result.severity = issue.severity;
  }
  for (const field of ['path', 'source', 'hint'] as const) {
    const value = issue[field];
    if (typeof value === 'string' || value === null) result[field] = value;
  }
  return result;
}

function normalizeError(
  response: Response,
  body: unknown,
): { code: string; message: string; issues: PublicApiIssue[] } {
  if (
    body != null &&
    typeof body === 'object' &&
    'issues' in body &&
    Array.isArray(body.issues) &&
    body.issues.length > 0 &&
    typeof body.issues[0] === 'object' &&
    body.issues[0] != null
  ) {
    const issues = body.issues
      .map(publicIssue)
      .filter((issue): issue is PublicApiIssue => issue !== null);
    const issue = issues[0];
    return {
      code: issue?.code ?? 'API_ERROR',
      message:
        issue && issue.message.length > 0
          ? issue.message
          : `Request failed with status ${response.status}`,
      issues,
    };
  }

  if (
    body != null &&
    typeof body === 'object' &&
    'detail' in body &&
    typeof body.detail === 'string'
  ) {
    return { code: 'API_ERROR', message: body.detail, issues: [] };
  }

  return {
    code: 'API_ERROR',
    message: `Request failed with status ${response.status}`,
    issues: [],
  };
}

async function request<T>(
  url: string,
  init: RequestInit,
  decode: (response: Response) => Promise<T>,
): Promise<T> {
  const baseUrl = getApiBaseUrl();
  const fullUrl = baseUrl ? `${baseUrl}${url}` : url;

  const response = await fetch(fullUrl, {
    ...init,
    credentials: 'omit',
  });

  if (!response.ok) {
    let body: unknown;
    try {
      body = await response.json();
    } catch {
      body = null;
    }
    const { code, message, issues } = normalizeError(response, body);
    throw new ApiError(response.status, code, message, issues);
  }

  return decode(response);
}

export const fetcher = async <T>(
  url: string,
  init: RequestInit = {},
): Promise<T> => request(url, init, (response) => response.json() as Promise<T>);

export const blobFetcher = async <T extends Blob = Blob>(
  url: string,
  init: RequestInit = {},
): Promise<T> => request(url, init, (response) => response.blob() as Promise<T>);
