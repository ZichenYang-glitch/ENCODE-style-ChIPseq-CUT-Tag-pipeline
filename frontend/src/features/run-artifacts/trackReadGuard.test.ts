import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  ensureTrackReadGuard,
  hasTrackReadFailure,
  resetTrackReadGuardForTests,
  subscribeTrackReadFailure,
} from './trackReadGuard';

const trackUrl =
  'https://helix.example/api/v1/runs/run-1/browser-tracks/artifact-pet/bytes/generation/revision';

afterEach(() => {
  resetTrackReadGuardForTests();
  vi.unstubAllGlobals();
});

describe('trackReadGuard', () => {
  it('reports a failed same-origin track read and still returns the response', async () => {
    const seen: string[] = [];
    const calls: string[] = [];
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL) => {
        calls.push(String(input));
        return new Response('no', { status: 409 });
      }),
    );
    ensureTrackReadGuard();
    const unsubscribe = subscribeTrackReadFailure((artifactId) => {
      seen.push(artifactId);
    });

    const response = await fetch(trackUrl);
    expect(response.status).toBe(409);
    expect(calls).toEqual([trackUrl]);
    expect(seen).toEqual(['artifact-pet']);
    expect(hasTrackReadFailure('artifact-pet')).toBe(true);
    unsubscribe();
  });

  it('keeps a successful range read and ignores unrelated requests', async () => {
    const seen: string[] = [];
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response(new Uint8Array([1]), { status: 206 })),
    );
    ensureTrackReadGuard();
    subscribeTrackReadFailure((artifactId) => {
      seen.push(artifactId);
    });

    expect((await fetch(trackUrl)).status).toBe(206);
    expect((await fetch('/api/v1/auth/session')).status).toBe(206);
    expect(seen).toEqual([]);
    expect(hasTrackReadFailure('artifact-pet')).toBe(false);
  });

  it('reports the embedded track when the proxy fetch is rejected', async () => {
    const seen: string[] = [];
    const proxy =
      'https://epigenome.wustl.edu/cors/127.0.0.1:14174/api/v1/runs/run-1/browser-tracks/artifact-peak/bytes/generation/revision.tbi';
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => {
        throw new TypeError('Failed to fetch');
      }),
    );
    ensureTrackReadGuard();
    subscribeTrackReadFailure((artifactId) => {
      seen.push(artifactId);
    });

    await expect(fetch(proxy)).rejects.toThrow('Failed to fetch');
    expect(seen).toEqual(['artifact-peak']);
  });
});
