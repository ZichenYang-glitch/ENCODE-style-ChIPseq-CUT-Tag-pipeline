export class ApiError extends Error {
  status: number;
  code: string;

  constructor(status: number, code: string, message: string) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.code = code;
  }
}

function getApiBaseUrl(): string {
  const envBase = import.meta.env?.VITE_API_BASE_URL;
  if (typeof envBase !== 'string' || envBase.trim() === '') {
    return '';
  }
  return envBase.replace(/\/$/, '');
}

function normalizeError(
  response: Response,
  body: unknown,
): { code: string; message: string } {
  if (
    body != null &&
    typeof body === 'object' &&
    'issues' in body &&
    Array.isArray(body.issues) &&
    body.issues.length > 0 &&
    typeof body.issues[0] === 'object' &&
    body.issues[0] != null
  ) {
    const issue = body.issues[0];
    return {
      code: typeof issue.code === 'string' ? issue.code : 'API_ERROR',
      message:
        typeof issue.message === 'string' && issue.message.length > 0
          ? issue.message
          : `Request failed with status ${response.status}`,
    };
  }

  if (
    body != null &&
    typeof body === 'object' &&
    'detail' in body &&
    typeof body.detail === 'string'
  ) {
    return { code: 'API_ERROR', message: body.detail };
  }

  return {
    code: 'API_ERROR',
    message: `Request failed with status ${response.status}`,
  };
}

export const fetcher = async <T>(url: string, init: RequestInit = {}): Promise<T> => {
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
    const { code, message } = normalizeError(response, body);
    throw new ApiError(response.status, code, message);
  }

  return response.json() as Promise<T>;
};
