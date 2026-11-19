import { createHash } from 'node:crypto';

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';

import { catalogDownloadRegisteredResultRevision } from './generated/registered-results/registered-results';
import { ApiError } from './fetcher';

const base = process.env.HELIXWEAVE_REGISTERED_RESULT_BASE;
const projectId = process.env.HELIXWEAVE_REGISTERED_RESULT_PROJECT;
const revisionId = process.env.HELIXWEAVE_REGISTERED_RESULT_REVISION;
const cookie = process.env.HELIXWEAVE_REGISTERED_RESULT_COOKIE;
const expectedSha = process.env.HELIXWEAVE_REGISTERED_RESULT_SHA256;
const enabled = Boolean(base && projectId && revisionId && cookie && expectedSha);

describe.skipIf(!enabled)('generated registered result download', () => {
  const realFetch = globalThis.fetch.bind(globalThis);

  beforeEach(() => {
    vi.stubEnv('VITE_API_BASE_URL', base ?? '');
    globalThis.fetch = (input: RequestInfo | URL, init?: RequestInit) => {
      const headers = new Headers(init?.headers);
      headers.set('cookie', cookie ?? '');
      return realFetch(input, { ...init, headers });
    };
  });

  afterEach(() => {
    globalThis.fetch = realFetch;
    vi.unstubAllEnvs();
  });

  it('returns the exact revision bytes from catalogDownloadRegisteredResultRevision', async () => {
    const blob = await catalogDownloadRegisteredResultRevision(
      projectId ?? '',
      revisionId ?? '',
    );
    const bytes = new Uint8Array(await blob.arrayBuffer());

    expect(blob).toBeInstanceOf(Blob);
    expect(bytes.byteLength).toBe(50);
    expect(createHash('sha256').update(bytes).digest('hex')).toBe(expectedSha);
  });

  it('decodes a missing revision as the JSON error contract', async () => {
    const missing = `rrev_${'0'.repeat(32)}`;

    try {
      await catalogDownloadRegisteredResultRevision(projectId ?? '', missing);
    } catch (error) {
      expect(error).toBeInstanceOf(ApiError);
      expect(error).toMatchObject({
        status: 404,
        code: 'catalog_not_found',
      });
      return;
    }
    throw new Error('missing revision returned a body');
  });
});
