const TRACK_BYTES = /\/browser-tracks\/([^/?#]+)\/bytes\//;
const GUARD = Symbol.for('helixweave.track-read-guard');

type Listener = (artifactId: string) => void;

type TrackRead = { artifactId: string; generation: number };

const listeners = new Set<Listener>();
const failed = new Set<string>();
const failedGeneration = new Map<string, number>();
const activeGeneration = new Map<string, number>();
const signalReads = new WeakMap<AbortSignal, TrackRead>();
let originalFetch: typeof fetch | null = null;
let violationListener: ((event: SecurityPolicyViolationEvent) => void) | null = null;

function artifactIdFromUrl(value: string): string | null {
  try {
    return new URL(value, window.location.origin).href.match(TRACK_BYTES)?.[1] ?? null;
  } catch {
    return null;
  }
}

function isAbortError(error: unknown): boolean {
  return (
    typeof error === 'object' &&
    error !== null &&
    'name' in error &&
    (error as { name?: unknown }).name === 'AbortError'
  );
}

function notify(artifactId: string, generation: number | undefined): void {
  const current = activeGeneration.get(artifactId);
  if (generation !== undefined && generation !== current) return;
  failedGeneration.set(artifactId, generation ?? current ?? 0);
  failed.add(artifactId);
  for (const listener of listeners) listener(artifactId);
}

export function beginTrackRead(artifactId: string, signal: AbortSignal): number {
  const generation = (activeGeneration.get(artifactId) ?? 0) + 1;
  activeGeneration.set(artifactId, generation);
  signalReads.set(signal, { artifactId, generation });
  if (failed.has(artifactId) && failedGeneration.get(artifactId) === 0) {
    failedGeneration.set(artifactId, generation);
  }
  return generation;
}

export function noteTrackReadSucceeded(artifactId: string, generation: number): boolean {
  if (activeGeneration.get(artifactId) !== generation) return false;
  if (failedGeneration.get(artifactId) === generation) return false;
  failed.delete(artifactId);
  failedGeneration.delete(artifactId);
  return true;
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
  failedGeneration.clear();
}

export function ensureTrackReadGuard(): void {
  const current = globalThis.fetch;
  if (typeof current !== 'function' || GUARD in current) return;
  const bound = current.bind(globalThis);
  originalFetch = bound;
  const guarded: typeof fetch = async (input, init) => {
    const boundRead = init?.signal ? signalReads.get(init.signal) : undefined;
    const artifactId = boundRead?.artifactId ?? artifactIdFromUrl(requestUrl(input));
    const generation = boundRead?.generation;
    const recordFailure = () => {
      if (!artifactId || init?.signal?.aborted) return;
      if (
        generation !== undefined &&
        activeGeneration.get(artifactId) !== generation
      ) {
        return;
      }
      notify(artifactId, generation);
    };
    try {
      const response = await bound(input, init);
      if (!response.ok) recordFailure();
      return response;
    } catch (error) {
      if (!isAbortError(error)) recordFailure();
      throw error;
    }
  };
  Object.defineProperty(guarded, GUARD, { value: true });
  globalThis.fetch = guarded;
  if (violationListener === null && typeof document !== 'undefined') {
    violationListener = (event) => {
      if (event.effectiveDirective !== 'connect-src') return;
      const artifactId = artifactIdFromUrl(event.blockedURI);
      if (artifactId) notify(artifactId, activeGeneration.get(artifactId));
    };
    document.addEventListener('securitypolicyviolation', violationListener);
  }
}

export function resetTrackReadGuardForTests(): void {
  failed.clear();
  failedGeneration.clear();
  activeGeneration.clear();
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
