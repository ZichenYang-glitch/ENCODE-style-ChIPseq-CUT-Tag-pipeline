import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { fetcher, blobFetcher, ApiError } from './fetcher';

describe('fetcher', () => {
  let originalFetch: typeof globalThis.fetch;

  beforeEach(() => {
    originalFetch = globalThis.fetch;
  });

  afterEach(() => {
    globalThis.fetch = originalFetch;
    vi.unstubAllEnvs();
  });

  it('prefixes VITE_API_BASE_URL when set', async () => {
    vi.stubEnv('VITE_API_BASE_URL', 'http://localhost:8000');
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => ({ ok: true, workflows: [], issues: [] }),
    } as Response);
    globalThis.fetch = fetchMock;

    await fetcher('/api/v1/workflows/', {});

    expect(fetchMock).toHaveBeenCalledWith(
      'http://localhost:8000/api/v1/workflows/',
      expect.objectContaining({ credentials: 'omit' }),
    );
  });

  it('uses a relative URL when VITE_API_BASE_URL is empty', async () => {
    vi.stubEnv('VITE_API_BASE_URL', '');
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => ({ ok: true, workflows: [], issues: [] }),
    } as Response);
    globalThis.fetch = fetchMock;

    await fetcher('/api/v1/workflows/', {});

    expect(fetchMock).toHaveBeenCalledWith(
      '/api/v1/workflows/',
      expect.objectContaining({ credentials: 'omit' }),
    );
  });

  it('strips a trailing slash from VITE_API_BASE_URL', async () => {
    vi.stubEnv('VITE_API_BASE_URL', 'http://localhost:8000/');
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => ({ ok: true, workflows: [], issues: [] }),
    } as Response);
    globalThis.fetch = fetchMock;

    await fetcher('/api/v1/workflows/', {});

    expect(fetchMock).toHaveBeenCalledWith(
      'http://localhost:8000/api/v1/workflows/',
      expect.objectContaining({ credentials: 'omit' }),
    );
  });

  it('normalizes backend issue envelopes into ApiError', async () => {
    vi.stubEnv('VITE_API_BASE_URL', '');
    const fetchMock = vi.fn().mockResolvedValue({
      ok: false,
      status: 404,
      json: async () => ({
        ok: false,
        run: null,
        issues: [
          {
            code: 'RUN_NOT_FOUND',
            message: 'Run was not found.',
            severity: 'error',
            path: 'run_id',
            source: 'repository',
            hint: 'Refresh the run list.',
            technical_message: '/private/path/database.sqlite failed',
            context: { database_path: '/private/path/database.sqlite' },
          },
        ],
      }),
    } as Response);
    globalThis.fetch = fetchMock;

    await expect(fetcher('/api/v1/runs/run-123', {})).rejects.toMatchObject({
      status: 404,
      code: 'RUN_NOT_FOUND',
      message: 'Run was not found.',
      issues: [
        {
          code: 'RUN_NOT_FOUND',
          message: 'Run was not found.',
          severity: 'error',
          path: 'run_id',
          source: 'repository',
          hint: 'Refresh the run list.',
        },
      ],
    });

    try {
      await fetcher('/api/v1/runs/run-123', {});
    } catch (error) {
      expect(JSON.stringify(error)).not.toContain('/private/path');
      expect(JSON.stringify(error)).not.toContain('database_path');
    }
  });

  it('does not leak raw response bodies when backend returns non-JSON', async () => {
    vi.stubEnv('VITE_API_BASE_URL', '');
    const fetchMock = vi.fn().mockResolvedValue({
      ok: false,
      status: 500,
      json: async () => {
        throw new Error('invalid json');
      },
    } as unknown as Response);
    globalThis.fetch = fetchMock;

    await expect(fetcher('/api/v1/runs/run-123', {})).rejects.toBeInstanceOf(
      ApiError,
    );
  });

  it('returns successful binary responses as Blob values', async () => {
    vi.stubEnv('VITE_API_BASE_URL', '');
    const source = new Blob(['durable artifact'], { type: 'text/plain' });
    globalThis.fetch = vi.fn().mockResolvedValue({
      ok: true,
      blob: async () => source,
    } as Response);

    const result = await blobFetcher('/api/v1/runs/run-1/artifacts/a/download');

    expect(result).toBe(source);
  });

  it('uses the same redacted Issue handling for binary request failures', async () => {
    vi.stubEnv('VITE_API_BASE_URL', '');
    globalThis.fetch = vi.fn().mockResolvedValue({
      ok: false,
      status: 409,
      json: async () => ({
        ok: false,
        issues: [
          {
            code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
            message: 'Artifact content is no longer available as indexed.',
            technical_message: '/private/workspace/source.tsv changed',
            context: { path: '/private/workspace/source.tsv' },
          },
        ],
      }),
    } as Response);

    await expect(
      blobFetcher('/api/v1/runs/run-1/artifacts/a/download'),
    ).rejects.toMatchObject({
      status: 409,
      code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
      issues: [
        {
          code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
          message: 'Artifact content is no longer available as indexed.',
        },
      ],
    });
  });
});
