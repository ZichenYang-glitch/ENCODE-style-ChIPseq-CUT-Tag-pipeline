import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  beginTrackRead,
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

  it('does not record a caller abort, including one that later returns an error', async () => {
    const seen: string[] = [];
    const controller = new AbortController();
    vi.stubGlobal(
      'fetch',
      vi.fn(
        async (_input: RequestInfo | URL, init?: RequestInit) =>
          await new Promise<Response>((_resolve, reject) => {
            init?.signal?.addEventListener(
              'abort',
              () => reject(new DOMException('cancelled', 'AbortError')),
              { once: true },
            );
          }),
      ),
    );
    ensureTrackReadGuard();
    subscribeTrackReadFailure((artifactId) => {
      seen.push(artifactId);
    });
    beginTrackRead('artifact-pet', controller.signal);
    const pending = fetch(trackUrl, { signal: controller.signal });
    controller.abort();
    await expect(pending).rejects.toMatchObject({ name: 'AbortError' });
    expect(seen).toEqual([]);
    expect(hasTrackReadFailure('artifact-pet')).toBe(false);
  });

  it('ignores a superseded generation after a newer read has started', async () => {
    const seen: string[] = [];
    let release: (response: Response) => void = () => {};
    let calls = 0;
    vi.stubGlobal(
      'fetch',
      vi.fn(() => {
        calls += 1;
        if (calls === 1) {
          return new Promise<Response>((resolve) => {
            release = resolve;
          });
        }
        return Promise.resolve(new Response('changed', { status: 409 }));
      }),
    );
    ensureTrackReadGuard();
    subscribeTrackReadFailure((artifactId) => {
      seen.push(artifactId);
    });
    const first = new AbortController();
    beginTrackRead('artifact-pet', first.signal);
    const firstRead = fetch(trackUrl, { signal: first.signal });
    const second = new AbortController();
    beginTrackRead('artifact-pet', second.signal);
    release(new Response('changed', { status: 409 }));
    expect((await firstRead).status).toBe(409);
    expect(seen).toEqual([]);
    expect(hasTrackReadFailure('artifact-pet')).toBe(false);

    const current = await fetch(trackUrl, { signal: second.signal });
    expect(current.status).toBe(409);
    expect(seen).toEqual(['artifact-pet']);
  });

  it('keeps a rejected read on one run from failing the same artifact id on another run', async () => {
    const seen: Array<{ artifactId: string; runId: string | null | undefined }> = [];
    const run1 =
      'https://helix.example/api/v1/runs/run-1/browser-tracks/artifact-pet/bytes/generation/revision';
    const run2 =
      'https://helix.example/api/v1/runs/run-2/browser-tracks/artifact-pet/bytes/generation/revision';
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL) =>
        String(input).includes('/runs/run-1/')
          ? new Response('no', { status: 409 })
          : new Response(new Uint8Array([1]), { status: 206 }),
      ),
    );
    ensureTrackReadGuard();
    subscribeTrackReadFailure((artifactId, detail) => {
      seen.push({ artifactId, runId: detail?.runId });
    });
    const first = new AbortController();
    beginTrackRead('artifact-pet', first.signal, 'run-1');
    const second = new AbortController();
    beginTrackRead('artifact-pet', second.signal, 'run-2');

    expect((await fetch(run1, { signal: first.signal })).status).toBe(409);
    expect((await fetch(run2, { signal: second.signal })).status).toBe(206);
    expect(seen).toEqual([{ artifactId: 'artifact-pet', runId: 'run-1' }]);
    expect(hasTrackReadFailure('artifact-pet')).toBe(true);
  });
});
