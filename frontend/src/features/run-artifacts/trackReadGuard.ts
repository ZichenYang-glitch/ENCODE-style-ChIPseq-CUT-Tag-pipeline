const TRACK_BYTES = /\/runs\/([^/?#]+)\/browser-tracks\/([^/?#]+)\/bytes\//;
const TRACK_BYTES_BARE = /\/browser-tracks\/([^/?#]+)\/bytes\//;
const GUARD = Symbol.for('helixweave.track-read-guard');

type TrackIdentity = { runId: string | null; artifactId: string };

type Listener = (artifactId: string, detail?: { runId: string | null }) => void;

type TrackRead = { artifactId: string; runId: string | null; generation: number };

type FailureRecord = { artifactId: string; runId: string | null; generation: number };

const listeners = new Set<Listener>();
const failures = new Map<string, FailureRecord>();
const activeGeneration = new Map<string, number>();
const signalReads = new WeakMap<AbortSignal, TrackRead>();
let originalFetch: typeof fetch | null = null;
let violationListener: ((event: SecurityPolicyViolationEvent) => void) | null = null;

export function trackReadKey(runId: string | null, artifactId: string): string {
  return runId ? `${runId}/${artifactId}` : artifactId;
}

function decodePart(value: string): string {
  try {
    return decodeURIComponent(value);
  } catch {
    return value;
  }
}

function identityFromUrl(value: string): TrackIdentity | null {
  try {
    const href = new URL(value, window.location.origin).href;
    const scoped = href.match(TRACK_BYTES);
    if (scoped) {
      return { runId: decodePart(scoped[1]), artifactId: decodePart(scoped[2]) };
    }
    const bare = href.match(TRACK_BYTES_BARE);
    if (!bare) return null;
    return { runId: null, artifactId: decodePart(bare[1]) };
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

function notify(artifactId: string, generation: number | undefined, runId: string | null): void {
  const key = trackReadKey(runId, artifactId);
  const current = activeGeneration.get(key);
  if (generation !== undefined && generation !== current) return;
  failures.set(key, {
    artifactId,
    runId,
    generation: generation ?? current ?? 0,
  });
  for (const listener of listeners) listener(artifactId, { runId });
}

export function beginTrackRead(
  artifactId: string,
  signal: AbortSignal,
  runId?: string | null,
): number {
  const scopedRun = runId ?? null;
  const key = trackReadKey(scopedRun, artifactId);
  const generation = (activeGeneration.get(key) ?? 0) + 1;
  activeGeneration.set(key, generation);
  signalReads.set(signal, { artifactId, runId: scopedRun, generation });
  const existing = failures.get(key);
  if (existing && existing.generation === 0) {
    failures.set(key, { ...existing, generation });
  }
  return generation;
}

export function noteTrackReadSucceeded(
  artifactId: string,
  generation: number,
  runId?: string | null,
): boolean {
  const key = trackReadKey(runId ?? null, artifactId);
  if (activeGeneration.get(key) !== generation) return false;
  if (failures.get(key)?.generation === generation) return false;
  failures.delete(key);
  return true;
}

function requestUrl(input: RequestInfo | URL): string {
  if (typeof input === 'string') return input;
  if (input instanceof URL) return input.toString();
  return input.url;
}

export function hasTrackReadFailure(artifactId: string): boolean {
  for (const record of failures.values()) {
    if (record.artifactId === artifactId) return true;
  }
  return false;
}

export function subscribeTrackReadFailure(listener: Listener): () => void {
  listeners.add(listener);
  for (const record of failures.values()) {
    listener(record.artifactId, { runId: record.runId });
  }
  return () => {
    listeners.delete(listener);
  };
}

export function clearTrackReadFailures(): void {
  failures.clear();
}

export function ensureTrackReadGuard(): void {
  const current = globalThis.fetch;
  if (typeof current !== 'function' || GUARD in current) return;
  const bound = current.bind(globalThis);
  originalFetch = bound;
  const guarded: typeof fetch = async (input, init) => {
    const boundRead = init?.signal ? signalReads.get(init.signal) : undefined;
    const fromUrl = identityFromUrl(requestUrl(input));
    const artifactId = boundRead?.artifactId ?? fromUrl?.artifactId ?? null;
    const runId = boundRead ? boundRead.runId : (fromUrl?.runId ?? null);
    const generation = boundRead?.generation;
    const recordFailure = () => {
      if (!artifactId || init?.signal?.aborted) return;
      if (
        generation !== undefined &&
        activeGeneration.get(trackReadKey(runId, artifactId)) !== generation
      ) {
        return;
      }
      notify(artifactId, generation, runId);
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
      const identity = identityFromUrl(event.blockedURI);
      if (!identity) return;
      const key = trackReadKey(identity.runId, identity.artifactId);
      notify(identity.artifactId, activeGeneration.get(key), identity.runId);
    };
    document.addEventListener('securitypolicyviolation', violationListener);
  }
}

export function resetTrackReadGuardForTests(): void {
  failures.clear();
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
