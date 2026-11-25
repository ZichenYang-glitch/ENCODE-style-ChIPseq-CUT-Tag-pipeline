const TRACK_BYTES = /(?:\/runs\/([^/]+))?\/browser-tracks\/([^/]+)\/bytes\/([^/]+\/[^/]+)$/;
const GUARD = Symbol.for('helixweave.track-read-guard');

type TrackIdentity = { runId: string | null; artifactId: string; version: string | null };

type Listener = (artifactId: string, detail?: { runId: string | null }) => void;

type TrackRead = TrackIdentity & { generation: number };

const listeners = new Set<Listener>();
const failures = new Map<string, TrackRead>();
const activeReads = new Map<string, TrackRead>();
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
    const match = new URL(value, window.location.origin).pathname.match(TRACK_BYTES);
    if (!match) return null;
    return {
      runId: match[1] ? decodePart(match[1]) : null,
      artifactId: decodePart(match[2]),
      version: match[3].replace(/\.tbi$/, ''),
    };
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

function notify(identity: TrackIdentity, generation?: number): void {
  const { artifactId, runId, version } = identity;
  const key = trackReadKey(runId, artifactId);
  const current = activeReads.get(key);
  if (generation !== undefined && generation !== current?.generation) return;
  if (version && current?.version && version !== current.version) return;
  failures.set(key, { ...identity, generation: generation ?? current?.generation ?? 0 });
  for (const listener of listeners) listener(artifactId, { runId });
}

export function beginTrackRead(
  artifactId: string,
  signal: AbortSignal,
  runId?: string | null,
  contentUrl?: string,
): number {
  const scopedRun = runId ?? null;
  const key = trackReadKey(scopedRun, artifactId);
  const generation = (activeReads.get(key)?.generation ?? 0) + 1;
  const version = contentUrl ? identityFromUrl(contentUrl)?.version ?? null : null;
  const read = { artifactId, runId: scopedRun, generation, version };
  activeReads.set(key, read);
  signalReads.set(signal, read);
  const existing = failures.get(key);
  if (existing?.generation === 0 && (!version || existing.version === version)) {
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
  if (activeReads.get(key)?.generation !== generation) return false;
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
    const signal = init?.signal ?? (input instanceof Request ? input.signal : null);
    const boundRead = signal ? signalReads.get(signal) : undefined;
    const identity = boundRead ?? identityFromUrl(requestUrl(input));
    const generation = boundRead?.generation;
    const recordFailure = () => {
      if (identity && !signal?.aborted) notify(identity, generation);
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
      notify(identity);
    };
    document.addEventListener('securitypolicyviolation', violationListener);
  }
}

export function resetTrackReadGuardForTests(): void {
  failures.clear();
  activeReads.clear();
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
