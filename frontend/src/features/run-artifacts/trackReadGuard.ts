const TRACK_BYTES = /\/browser-tracks\/([^/?#]+)\/bytes\//;
const GUARD = Symbol.for('helixweave.track-read-guard');

type Listener = (artifactId: string) => void;

const listeners = new Set<Listener>();
const failed = new Set<string>();
let originalFetch: typeof fetch | null = null;
let violationListener: ((event: SecurityPolicyViolationEvent) => void) | null = null;

function artifactIdFromUrl(value: string): string | null {
  try {
    return new URL(value, window.location.origin).href.match(TRACK_BYTES)?.[1] ?? null;
  } catch {
    return null;
  }
}

function notify(artifactId: string): void {
  failed.add(artifactId);
  for (const listener of listeners) listener(artifactId);
}

function requestUrl(input: RequestInfo | URL): string {
  if (typeof input === 'string') return input;
  if (input instanceof URL) return input.toString();
  return input.url;
}

export function hasTrackReadFailure(artifactId: string): boolean {
  return failed.has(artifactId);
}

export function subscribeTrackReadFailure(listener: Listener): () => void {
  listeners.add(listener);
  for (const artifactId of failed) listener(artifactId);
  return () => {
    listeners.delete(listener);
  };
}

export function clearTrackReadFailures(): void {
  failed.clear();
}

export function ensureTrackReadGuard(): void {
  const current = globalThis.fetch;
  if (typeof current !== 'function' || GUARD in current) return;
  const bound = current.bind(globalThis);
  originalFetch = bound;
  const guarded: typeof fetch = async (input, init) => {
    const artifactId = artifactIdFromUrl(requestUrl(input));
    try {
      const response = await bound(input, init);
      if (artifactId && !response.ok) notify(artifactId);
      return response;
    } catch (error) {
      if (artifactId) notify(artifactId);
      throw error;
    }
  };
  Object.defineProperty(guarded, GUARD, { value: true });
  globalThis.fetch = guarded;
  if (violationListener === null && typeof document !== 'undefined') {
    violationListener = (event) => {
      if (event.effectiveDirective !== 'connect-src') return;
      const artifactId = artifactIdFromUrl(event.blockedURI);
      if (artifactId) notify(artifactId);
    };
    document.addEventListener('securitypolicyviolation', violationListener);
  }
}

export function resetTrackReadGuardForTests(): void {
  failed.clear();
  listeners.clear();
  const current = globalThis.fetch;
  if (originalFetch && typeof current === 'function' && GUARD in current) {
    globalThis.fetch = originalFetch;
  }
  originalFetch = null;
  if (violationListener !== null && typeof document !== 'undefined') {
    document.removeEventListener('securitypolicyviolation', violationListener);
    violationListener = null;
  }
}
