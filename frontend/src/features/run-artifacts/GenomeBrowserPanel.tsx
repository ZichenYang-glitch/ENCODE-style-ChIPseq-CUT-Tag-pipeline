import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useInfiniteQuery, useMutation, useQueries, useQuery, useQueryClient } from '@tanstack/react-query';
import { GenomeHub, type WashuTrackConfig } from 'wuepgg';
import washUStyles from 'wuepgg/style.css?raw';
import {
  getRunBrowserSession,
  prepareRunBrowserTrack,
  summarizeRunBrowserTrackRegion,
} from '../../api/generated/artifacts/artifacts';
import type {
  BrowserSessionResponse,
  BrowserTrackResponse,
  RunSummaryResponse,
} from '../../api/generated/models';
import { listRuns } from '../../api/generated/runs/runs';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import LocalTrackPicker, { type LocalWashuTrack } from '../../viewer/local/LocalTrackPicker';
import LocusSearch from '../../viewer/search/LocusSearch';
import {
  declaredAssembly,
  localTrackHasForeignUrl,
  placeLocalTrack,
  placePlatformRun,
  referenceBinding,
  seenWindow,
  userSpecifiedReference,
  USER_ASSEMBLIES,
  windowAuthority,
  type WindowAuthority,
} from './browserWindow';
import {
  beginTrackRead,
  clearTrackReadFailures,
  ensureTrackReadGuard,
  noteTrackReadSucceeded,
  subscribeTrackReadFailure,
  trackReadKey,
} from './trackReadGuard';

const REGION = /^([A-Za-z0-9_.:-]+):([0-9]+)-([0-9]+)$/;
const DRAW_NOTE =
  'WashU draws at most 100000 rows in one window. A larger count is not fully drawn.';
const EMPTY_RUNS: string[] = [];

type ChosenTrack = {
  artifactId: string;
  revision: string;
};

type TrackSelection = {
  runId: string;
  generation: string | null;
  chosen: ChosenTrack[];
  stale: string[];
};

type RunSelection = {
  generation: string | null;
  chosen: ChosenTrack[];
  stale: string[];
};

type PlatformRow = {
  runId: string;
  generation: string | null;
  track: BrowserTrackResponse;
};

function isGenerationConflict(error: unknown): boolean {
  return (
    error instanceof ApiError &&
    error.status === 409 &&
    error.code === 'BROWSER_TRACK_READ_CONFLICT' &&
    error.issues.some((issue) => issue.path === 'generation')
  );
}

function trackCanBeChosen(track: BrowserTrackResponse): boolean {
  return (
    track.preparation !== 'reference_mismatch' &&
    track.preparation !== 'source_too_large'
  );
}

function presetChosen(
  session: BrowserSessionResponse,
  artifactIds: readonly string[] | null,
): ChosenTrack[] {
  if (artifactIds === null || artifactIds.length === 0) return [];
  const tracks = new Map(session.tracks.map((track) => [track.artifact_id, track]));
  const chosen: ChosenTrack[] = [];
  for (const artifactId of artifactIds) {
    const track = tracks.get(artifactId);
    if (!track || !trackCanBeChosen(track)) continue;
    if (chosen.some((item) => item.artifactId === artifactId)) continue;
    chosen.push({ artifactId, revision: track.revision });
  }
  return chosen;
}

function reconcileChosen(
  current: { chosen: ChosenTrack[]; stale: string[] },
  session: BrowserSessionResponse,
): { chosen: ChosenTrack[]; stale: string[] } | null {
  const tracks = new Map(session.tracks.map((track) => [track.artifact_id, track]));
  const chosen: ChosenTrack[] = [];
  const stale: string[] = [];
  for (const item of current.chosen) {
    const track = tracks.get(item.artifactId);
    if (!track || !trackCanBeChosen(track)) continue;
    if (track.revision !== item.revision) {
      stale.push(item.artifactId);
      continue;
    }
    chosen.push(item);
  }
  const sameChosen =
    chosen.length === current.chosen.length &&
    chosen.every((item, index) => item === current.chosen[index]);
  const sameStale =
    stale.length === current.stale.length &&
    stale.every((artifactId, index) => artifactId === current.stale[index]);
  if (sameChosen && sameStale) return null;
  return { chosen, stale };
}

function reconcileSelection(
  current: TrackSelection,
  session: BrowserSessionResponse,
): TrackSelection {
  const generation = session.artifact_generation ?? null;
  if (current.generation !== generation) {
    return { runId: current.runId, generation, chosen: [], stale: [] };
  }
  const next = reconcileChosen(current, session);
  if (next === null) return current;
  return { ...current, ...next };
}

function sameOriginPath(path: string | null | undefined): string | null {
  if (!path || !path.startsWith('/api/')) return null;
  const configured = import.meta.env.VITE_API_BASE_URL;
  if (typeof configured === 'string' && configured.trim() !== '') {
    try {
      const url = new URL(configured, window.location.origin);
      if (url.origin !== window.location.origin) return null;
    } catch {
      return null;
    }
  }
  return new URL(path, window.location.origin).toString();
}

function messageFromError(error: unknown): string {
  if (error instanceof ApiError && error.message.trim() !== '') return error.message;
  return 'The browser view could not be loaded.';
}

function sameLocal(current: LocalWashuTrack[], next: LocalWashuTrack[]): boolean {
  if (current.length !== next.length) return false;
  return current.every((track, index) => {
    const other = next[index];
    return (
      other !== undefined &&
      track.id === other.id &&
      track.label === other.label &&
      track.assembly === other.assembly &&
      track.model === other.model
    );
  });
}

function sessionReadyTracks(
  session: BrowserSessionResponse,
  selected: readonly string[],
): BrowserTrackResponse[] {
  return session.tracks.filter(
    (track) =>
      selected.includes(track.artifact_id) &&
      track.preparation === 'ready' &&
      sameOriginPath(track.content_path) !== null,
  );
}

function activeChosen(
  selection: RunSelection | null,
  session: BrowserSessionResponse,
): ChosenTrack[] {
  if (selection === null || selection.generation !== (session.artifact_generation ?? null)) {
    return [];
  }
  return selection.chosen.filter((item) =>
    session.tracks.some(
      (track) =>
        track.artifact_id === item.artifactId &&
        track.revision === item.revision &&
        trackCanBeChosen(track),
    ),
  );
}

export default function GenomeBrowserPanel({
  runId = null,
  initialRegion = null,
  artifactIds = null,
  showToolbar = true,
  onReady,
  onFailed,
}: {
  runId?: string | null;
  initialRegion?: string | null;
  artifactIds?: readonly string[] | null;
  showToolbar?: boolean;
  onReady?: (notice: {
    runId: string;
    referenceStatus: 'bound' | 'unbound' | 'unavailable';
    artifactIds: string[];
  }) => void;
  onFailed?: (message: string) => void;
}) {
  ensureTrackReadGuard();
  const queryClient = useQueryClient();
  const openedRunId = runId || null;
  const sessionQuery = useQuery({
    queryKey: ['browser-session', openedRunId],
    enabled: openedRunId !== null,
    queryFn: ({ signal }) => getRunBrowserSession(openedRunId as string, { signal }),
    // A generation conflict is the user's prompt to refresh. Retrying would
    // load the new catalog by itself and hide that choice.
    retry: (failureCount, error) =>
      !isGenerationConflict(error) && failureCount < 3,
    // Tests observe the exhausted retries immediately. Production keeps the backoff.
    retryDelay: import.meta.env.MODE === 'test' ? 0 : undefined,
  });
  const session =
    sessionQuery.data?.run_id === openedRunId ? sessionQuery.data : undefined;
  const [selection, setSelection] = useState<TrackSelection | null>(null);
  const [region, setRegion] = useState('');
  const [loadFailures, setLoadFailures] = useState<Record<string, string>>({});
  const [userAssembly, setUserAssembly] = useState('hg38');
  const [assemblyRegions, setAssemblyRegions] = useState<Record<string, string>>({});
  const [addedRunIds, setAddedRunIds] = useState<string[]>([]);
  const [addedSelections, setAddedSelections] = useState<Record<string, RunSelection>>({});
  const [pendingRunId, setPendingRunId] = useState('');
  const [showRunPicker, setShowRunPicker] = useState(false);
  const browseRuns = openedRunId !== null || showRunPicker;
  const [runNotice, setRunNotice] = useState<string | null>(null);
  const [localTracks, setLocalTracks] = useState<LocalWashuTrack[]>([]);
  const scope = openedRunId ?? '';
  const [scopeKey, setScopeKey] = useState(scope);
  if (scopeKey !== scope) {
    setScopeKey(scope);
    setRegion('');
    setAssemblyRegions({});
    setLocalTracks([]);
    setAddedRunIds([]);
    setAddedSelections({});
    setPendingRunId('');
    setRunNotice(null);
    setLoadFailures({});
  }
  const scopeReady = scopeKey === scope;
  const activeAddedIds = scopeReady ? addedRunIds : EMPTY_RUNS;
  const activeLocal = scopeReady ? localTracks : EMPTY_LOCAL;
  const selectionRef = useRef<TrackSelection | null>(null);
  const pendingToggles = useRef<BrowserTrackResponse[]>([]);
  const presetRef = useRef<readonly string[] | null>(artifactIds ?? null);
  const conflictRef = useRef<string | null>(null);
  const regionKey = useRef<string | null>(null);
  selectionRef.current = selection;
  presetRef.current = artifactIds ?? null;
  if (conflictRef.current !== null && conflictRef.current !== openedRunId) {
    conflictRef.current = null;
  }
  if (isGenerationConflict(sessionQuery.error)) conflictRef.current = openedRunId;

  const onReadyRef = useRef(onReady);
  const onFailedRef = useRef(onFailed);
  onReadyRef.current = onReady;
  onFailedRef.current = onFailed;
  const onLocalTracks = useCallback((next: LocalWashuTrack[]) => {
    setLocalTracks((current) => (sameLocal(current, next) ? current : next));
  }, []);
  const regionHops = useRef<number[]>([]);
  const onNavigate = useCallback((next: string) => {
    for (const id of regionHops.current) window.clearTimeout(id);
    regionHops.current = [];
    setRegion(next.trim());
  }, []);

  useEffect(
    () =>
      subscribeTrackReadFailure((artifactId, detail) => {
        const failedRun = detail?.runId ?? null;
        if (!failedRun) return;
        const key = trackReadKey(failedRun, artifactId);
        setLoadFailures((current) =>
          current[key] ? current : { ...current, [key]: 'Track failed to load.' },
        );
      }),
    [],
  );
  useEffect(() => () => clearTrackReadFailures(), []);

  const browserKey = openedRunId ?? 'standalone';
  useEffect(() => {
    const style = document.createElement('style');
    style.dataset.helixweaveWashu = browserKey;
    style.textContent = washUStyles;
    document.head.appendChild(style);
    return () => {
      style.remove();
      localStorage.removeItem(`persist:helixweave-browser-${browserKey}`);
    };
  }, [browserKey]);

  useEffect(() => {
    if (sessionQuery.isError && !isGenerationConflict(sessionQuery.error)) {
      onFailedRef.current?.(messageFromError(sessionQuery.error));
    }
  }, [sessionQuery.error, sessionQuery.isError]);

  useEffect(() => {
    if (!openedRunId) return;
    if (sessionQuery.isError) {
      const current = selectionRef.current;
      if (
        current &&
        current.runId === openedRunId &&
        (current.chosen.length > 0 || current.stale.length > 0)
      ) {
        const next = { runId: openedRunId, generation: current.generation, chosen: [], stale: [] };
        selectionRef.current = next;
        setSelection(next);
      }
      return;
    }
    if (!session) return;
    const generation = session.artifact_generation ?? null;
    const current = selectionRef.current;
    let next: TrackSelection;
    if (current === null || current.runId !== openedRunId) {
      const skipPreset = conflictRef.current === openedRunId;
      if (skipPreset) conflictRef.current = null;
      next = {
        runId: openedRunId,
        generation,
        chosen: skipPreset ? [] : presetChosen(session, presetRef.current),
        stale: [],
      };
    } else {
      next = reconcileSelection(current, session);
    }
    if (next !== current) {
      selectionRef.current = next;
      setSelection(next);
    }
    const referenceStatus =
      session.reference.status === 'bound' ||
      session.reference.status === 'unbound' ||
      session.reference.status === 'unavailable'
        ? session.reference.status
        : 'unavailable';
    onReadyRef.current?.({
      runId: session.run_id,
      referenceStatus,
      artifactIds: next.chosen.map((item) => item.artifactId),
    });
  }, [openedRunId, session, sessionQuery.error, sessionQuery.isError]);

  useEffect(() => {
    if (!session || !openedRunId) return;
    const key = `${openedRunId}:${session.artifact_generation ?? ''}:${initialRegion ?? ''}`;
    if (regionKey.current === key) return;
    regionKey.current = key;
    const suggested = session.tracks.find(
      (track) => track.provenance?.suggested_region,
    )?.provenance?.suggested_region;
    const target =
      initialRegion && REGION.test(initialRegion)
        ? initialRegion
        : typeof suggested === 'string'
          ? suggested
          : null;
    if (target === null) return;
    // GenomeHub applies the region it receives on the first mount, then only
    // refetches tracks when that string changes. A short hop makes the mounted
    // viewer query the intended window. A later gene or coordinate edit clears
    // both timers, so they cannot pull the window back.
    const parts = target.match(/^(.*:)(\d+)-(\d+)$/);
    const hop = parts
      ? `${parts[1]}${Math.max(0, Number(parts[2]) - 400)}-${Math.max(1, Number(parts[3]) - 400)}`
      : target;
    setRegion(target);
    const move = window.setTimeout(() => setRegion(hop), 600);
    const settle = window.setTimeout(() => setRegion(target), 1100);
    regionHops.current = [move, settle];
    return () => {
      window.clearTimeout(move);
      window.clearTimeout(settle);
      regionHops.current = regionHops.current.filter((id) => id !== move && id !== settle);
    };
  }, [initialRegion, openedRunId, session]);

  const runsQuery = useInfiniteQuery({
    queryKey: ['browser-authorized-runs'],
    enabled: browseRuns,
    initialPageParam: undefined as string | undefined,
    queryFn: async ({ pageParam, signal }) => {
      const page = await listRuns(
        { limit: 100, after: typeof pageParam === 'string' ? pageParam : undefined },
        { signal },
      );
      if (!page || page.ok !== true || !Array.isArray(page.runs)) {
        throw new Error('Authorized runs could not be loaded.');
      }
      return page;
    },
    getNextPageParam: (lastPage, pages, lastParam) => {
      void pages;
      const cursor = lastPage.next_cursor;
      if (typeof cursor !== 'string' || cursor === '' || cursor === lastParam) return undefined;
      return cursor;
    },
    retry: import.meta.env.MODE === 'test' ? false : 1,
  });

  const addedQueryDefs = useMemo(
    () =>
      activeAddedIds.map((id) => ({
        queryKey: ['browser-session', id] as const,
        queryFn: ({ signal }: { signal: AbortSignal }) => getRunBrowserSession(id, { signal }),
        retry: (failureCount: number, error: unknown) =>
          !isGenerationConflict(error) && failureCount < 3,
        retryDelay: import.meta.env.MODE === 'test' ? 0 : undefined,
      })),
    [activeAddedIds],
  );
  const addedSessionQueries = useQueries({ queries: addedQueryDefs });
  const addedSessionKey = activeAddedIds
    .map((id, index) => {
      const query = addedSessionQueries[index];
      const data = query?.data?.run_id === id ? query.data : undefined;
      if (!data) return `${id}:${query?.isError ? 'error' : 'pending'}`;
      return [
        id,
        data.artifact_generation ?? '',
        data.reference.status,
        data.reference.identity_sha256 ?? '',
        data.reference.viewer_scaffold ?? '',
        data.tracks
          .map(
            (track) =>
              `${track.artifact_id}:${track.revision}:${track.preparation}:${track.content_path ?? ''}:${track.index_path ?? ''}`,
          )
          .join(','),
      ].join(':');
    })
    .join('|');
  const addedSessionsRef = useRef(addedSessionQueries);
  addedSessionsRef.current = addedSessionQueries;

  useEffect(() => {
    const queries = addedSessionsRef.current;
    setAddedSelections((current) => {
      let changed = false;
      const next = { ...current };
      activeAddedIds.forEach((id, index) => {
        const data = queries[index]?.data;
        if (!data || data.run_id !== id) return;
        const generation = data.artifact_generation ?? null;
        const existing = current[id];
        if (!existing || existing.generation !== generation) {
          next[id] = { generation, chosen: [], stale: [] };
          changed = true;
          return;
        }
        const reconciled = reconcileChosen(existing, data);
        if (reconciled === null) return;
        next[id] = { ...existing, ...reconciled };
        changed = true;
      });
      for (const id of Object.keys(next)) {
        if (!activeAddedIds.includes(id)) {
          delete next[id];
          changed = true;
        }
      }
      return changed ? next : current;
    });
  }, [activeAddedIds, addedSessionKey]);

  const generationConflict = isGenerationConflict(sessionQuery.error);
  const activeChoice = useMemo(() => {
    if (
      generationConflict ||
      sessionQuery.isError ||
      session === undefined ||
      session.ok !== true ||
      selection === null ||
      selection.runId !== openedRunId ||
      selection.generation !== (session.artifact_generation ?? null)
    ) {
      return [];
    }
    return selection.chosen.filter((item) =>
      session.tracks.some(
        (track) =>
          track.artifact_id === item.artifactId &&
          track.revision === item.revision &&
          trackCanBeChosen(track),
      ),
    );
  }, [generationConflict, openedRunId, selection, session, sessionQuery.isError]);
  const selected = useMemo(
    () => activeChoice.map((item) => item.artifactId),
    [activeChoice],
  );
  const stale =
    session !== undefined &&
    selection !== null &&
    selection.runId === openedRunId &&
    selection.generation === (session.artifact_generation ?? null)
      ? selection.stale
      : [];

  function toggleTrack(track: BrowserTrackResponse) {
    pendingToggles.current.push(track);
    if (pendingToggles.current.length !== 1) return;
    // Clicks in one turn share a mount. Adding a bigWig to a live hub crashes it,
    // and the index fallback is only attempted from the hub that first receives the file.
    queueMicrotask(() => {
      const batch = pendingToggles.current.splice(0);
      setSelection((current) => {
        const generation = session?.artifact_generation ?? null;
        let next: TrackSelection =
          current !== null && current.runId === openedRunId
            ? current
            : { runId: openedRunId ?? '', generation, chosen: [], stale: [] };
        for (const item of batch) {
          const exists = next.chosen.some((chosen) => chosen.artifactId === item.artifact_id);
          next = {
            ...next,
            chosen: exists
              ? next.chosen.filter((chosen) => chosen.artifactId !== item.artifact_id)
              : [...next.chosen, { artifactId: item.artifact_id, revision: item.revision }],
            stale: next.stale.filter((artifactId) => artifactId !== item.artifact_id),
          };
        }
        selectionRef.current = next;
        return next;
      });
    });
  }

  const readyTracks = useMemo(
    () => (session ? sessionReadyTracks(session, selected) : []),
    [selected, session],
  );
  const authority = useMemo(() => {
    const primary =
      session && openedRunId
        ? { runId: session.run_id, reference: referenceBinding(session.reference) }
        : null;
    return windowAuthority(primary, userAssembly);
  }, [openedRunId, session, userAssembly]);
  const primaryPlacement =
    session && openedRunId
      ? placePlatformRun(authority, session.run_id, referenceBinding(session.reference))
      : null;
  const primaryRows = useMemo<PlatformRow[]>(() => {
    if (!session || !primaryPlacement?.ok) return [];
    return readyTracks.map((track) => ({
      runId: session.run_id,
      generation: session.artifact_generation ?? null,
      track,
    }));
  }, [primaryPlacement?.ok, readyTracks, session]);
  const addedReady = useMemo<PlatformRow[]>(() => {
    const rows: PlatformRow[] = [];
    activeAddedIds.forEach((id, index) => {
      const query = addedSessionQueries[index];
      const data = query?.data;
      if (!data || data.run_id !== id || !data.ok || isGenerationConflict(query?.error)) return;
      const placement = placePlatformRun(authority, id, referenceBinding(data.reference));
      if (!placement.ok) return;
      const chosen = activeChosen(addedSelections[id] ?? null, data);
      for (const track of sessionReadyTracks(
        data,
        chosen.map((item) => item.artifactId),
      )) {
        rows.push({ runId: id, generation: data.artifact_generation ?? null, track });
      }
    });
    return rows;
    // addedSessionKey stands in for the query payloads.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeAddedIds, addedSelections, addedSessionKey, authority]);
  const platformRows = useMemo(
    () => [...primaryRows, ...addedReady],
    [addedReady, primaryRows],
  );

  useEffect(() => {
    const controllers = new Map<string, AbortController>();
    let cancelled = false;
    async function preflight(row: PlatformRow) {
      const path = sameOriginPath(row.track.content_path);
      if (!path) return;
      const key = trackReadKey(row.runId, row.track.artifact_id);
      const controller = new AbortController();
      controllers.set(key, controller);
      const generation = beginTrackRead(row.track.artifact_id, controller.signal, row.runId);
      try {
        const response = await fetch(path, {
          method: 'GET',
          headers: { Range: 'bytes=0-0' },
          credentials: 'same-origin',
          signal: controller.signal,
        });
        if (!response.ok && response.status !== 206) {
          throw await trackReadError(response);
        }
        const indexPath = sameOriginPath(row.track.index_path);
        if (indexPath) {
          const indexResponse = await fetch(indexPath, {
            method: 'GET',
            headers: { Range: 'bytes=0-1' },
            credentials: 'same-origin',
            signal: controller.signal,
          });
          if (!indexResponse.ok && indexResponse.status !== 206) {
            throw await trackReadError(indexResponse);
          }
        }
        if (!cancelled && noteTrackReadSucceeded(row.track.artifact_id, generation, row.runId)) {
          setLoadFailures((current) => {
            if (!current[key]) return current;
            const next = { ...current };
            delete next[key];
            return next;
          });
        }
      } catch (error) {
        if (controller.signal.aborted || cancelled) return;
        if (error instanceof DOMException && error.name === 'AbortError') return;
        setLoadFailures((current) => ({
          ...current,
          [key]:
            error instanceof Error && error.message.startsWith('Track failed to load (HTTP ')
              ? error.message
              : 'Track failed to load.',
        }));
      }
    }
    for (const row of platformRows) void preflight(row);
    return () => {
      cancelled = true;
      for (const controller of controllers.values()) controller.abort();
    };
  }, [platformRows, region]);

  const visibleRows = useMemo(
    () => platformRows.filter(row => !loadFailures[trackReadKey(row.runId, row.track.artifact_id)]),
    [platformRows, loadFailures],
  );
  const failedRows = platformRows.filter(
    (row) => loadFailures[trackReadKey(row.runId, row.track.artifact_id)],
  );

  const prepare = useMutation({
    mutationFn: async (track: BrowserTrackResponse) => {
      if (!openedRunId || !session?.artifact_generation) {
        throw new ApiError(409, 'BROWSER_TRACK_READ_CONFLICT', 'Reload the browser view.');
      }
      return prepareRunBrowserTrack(openedRunId, track.artifact_id, {
        generation: session.artifact_generation,
        revision: track.revision,
      });
    },
    onSuccess: async () => {
      await queryClient.invalidateQueries({ queryKey: ['browser-session', openedRunId] });
    },
  });

  const showAssembly = !(
    session &&
    session.reference.status === 'bound' &&
    session.reference.viewer_scaffold
  );
  const authorizedRuns = useMemo(() => {
    const seen = new Set<string>();
    const runs: RunSummaryResponse[] = [];
    for (const page of runsQuery.data?.pages ?? []) {
      for (const run of page.runs) {
        if (seen.has(run.run_id)) continue;
        seen.add(run.run_id);
        runs.push(run);
      }
    }
    return runs;
  }, [runsQuery.data]);

  function addRun() {
    if (!pendingRunId) return;
    if (pendingRunId === openedRunId || activeAddedIds.includes(pendingRunId)) {
      setRunNotice(`${pendingRunId} is already in this window.`);
      return;
    }
    setRunNotice(null);
    setAddedRunIds((current) =>
      current.includes(pendingRunId) ? current : [...current, pendingRunId],
    );
    setPendingRunId('');
  }

  function removeRun(id: string) {
    setAddedRunIds((current) => current.filter((item) => item !== id));
    setAddedSelections((current) => {
      if (!(id in current)) return current;
      const next = { ...current };
      delete next[id];
      return next;
    });
  }

  function toggleAddedTrack(id: string, sessionForRun: BrowserSessionResponse, track: BrowserTrackResponse) {
    setAddedSelections((current) => {
      const generation = sessionForRun.artifact_generation ?? null;
      const existing = current[id] ?? { generation, chosen: [], stale: [] };
      const exists = existing.chosen.some((item) => item.artifactId === track.artifact_id);
      return {
        ...current,
        [id]: {
          generation: existing.generation ?? generation,
          chosen: exists
            ? existing.chosen.filter((item) => item.artifactId !== track.artifact_id)
            : [...existing.chosen, { artifactId: track.artifact_id, revision: track.revision }],
          stale: existing.stale.filter((artifactId) => artifactId !== track.artifact_id),
        },
      };
    });
  }

  const primaryView: 'loading' | 'conflict' | 'error' | 'ready' | 'none' = !openedRunId
    ? 'none'
    : generationConflict
      ? 'conflict'
      : sessionQuery.isError || (session !== undefined && session.ok !== true)
        ? 'error'
        : session === undefined
          ? 'loading'
          : 'ready';

  const localPlaced = useMemo(() => activeLocal.map(track => ({
    track,
    placement: placeLocalTrack(authority, track),
    foreign: localTrackHasForeignUrl(track.model, window.location.origin),
  })), [activeLocal, authority]);
  const hubLocal = useMemo(() => localPlaced.flatMap(item =>
    item.placement.ok && !item.foreign ? [{ id: item.track.id, model: item.track.model }] : [],
  ), [localPlaced]);

  return (
    <section className="min-w-0 space-y-3" aria-label="Genome browser">
      {session && primaryView === 'ready' ? <ReferenceNotice session={session} /> : null}
      {session && primaryView === 'ready' ? <InteractionNotice session={session} /> : null}
      {showAssembly ? (
        <label className="block text-sm">
          Assembly
          <select
            className="mt-1 block border border-[var(--color-border)] px-2 py-1"
            aria-label="Assembly"
            value={userAssembly}
            onChange={(event) => {
              const next = event.target.value;
              setAssemblyRegions(current => ({ ...current, [userAssembly]: region }));
              setUserAssembly(next);
              onNavigate(assemblyRegions[next] ?? '');
            }}
          >
            {USER_ASSEMBLIES.map((item) => (
              <option key={item.scaffold} value={item.scaffold}>
                {item.scaffold} · {item.declared}
              </option>
            ))}
          </select>
        </label>
      ) : null}
      {showAssembly ? (
        <p role="status">当前声明 {declaredAssembly(userAssembly)}</p>
      ) : null}
      {authority.source === 'user' ? (
        <p role="status">{userSpecifiedReference(authority.scaffold)}</p>
      ) : null}
      {openedRunId && primaryView === 'loading' ? (
        <p role="status">Loading genome browser…</p>
      ) : null}
      {openedRunId && primaryView !== 'loading' ? (
        <article
          className="space-y-2"
          aria-label={`Run ${session?.run_id ?? openedRunId}`}
          data-testid="browser-run"
          data-run-id={session?.run_id ?? openedRunId}
          data-run-source="platform"
          data-generation={primaryView === 'ready' ? (session?.artifact_generation ?? '') : ''}
          data-reference-identity={
            primaryView === 'ready' ? (session?.reference.identity_sha256 ?? '') : ''
          }
          data-scaffold={primaryView === 'ready' ? (session?.reference.viewer_scaffold ?? '') : ''}
        >
          {primaryView === 'conflict' ? (
            <div role="alert">
              <p>
                These files changed. Refresh the list and choose the tracks again. New
                files are not read until you choose them.
              </p>
            </div>
          ) : null}
          {primaryView === 'error' ? (
            <div role="alert">
              <p>{messageFromError(sessionQuery.error)}</p>
            </div>
          ) : null}
          {primaryView === 'ready' && session ? (
            <>
              <p>
                Source {session.workflow_id} · run {session.run_id} · generation{' '}
                {session.artifact_generation ?? 'unknown'}
                {session.reference.revision_id
                  ? ` · reference revision ${session.reference.revision_id}`
                  : ''}
              </p>
              <TrackChooser
                runId={session.run_id}
                generation={session.artifact_generation ?? null}
                tracks={session.tracks}
                unavailable={session.unavailable}
                selected={primaryPlacement?.ok ? selected : []}
                stale={stale}
                selectable={primaryPlacement?.ok === true}
                onToggle={toggleTrack}
              />
              <div className="flex flex-wrap gap-2">
                {session.tracks
                  .filter(
                    (track) =>
                      selected.includes(track.artifact_id) && track.preparation === 'needed',
                  )
                  .map((track) => (
                    <Button
                      key={track.artifact_id}
                      type="button"
                      variant="secondary"
                      onClick={() => prepare.mutate(track)}
                      disabled={prepare.isPending}
                    >
                      Prepare {track.label}
                    </Button>
                  ))}
              </div>
              {prepare.isError ? <p role="alert">{messageFromError(prepare.error)}</p> : null}
            </>
          ) : null}
          {primaryView === 'conflict' || primaryView === 'ready' ? (
            <Button type="button" variant="secondary" onClick={() => void sessionQuery.refetch()}>
              Refresh list
            </Button>
          ) : null}
          {primaryView === 'error' ? (
            <Button type="button" variant="secondary" onClick={() => void sessionQuery.refetch()}>
              Retry browser view
            </Button>
          ) : null}
        </article>
      ) : null}
      <div data-testid="local-track-slot">
        <LocalTrackPicker assembly={authority.scaffold} onTracksChange={onLocalTracks} />
      </div>
      {browseRuns ? (
        <RunPicker
          runs={authorizedRuns}
          pendingRunId={pendingRunId}
          onPending={setPendingRunId}
          onAdd={addRun}
          onMore={() => void runsQuery.fetchNextPage()}
          hasMore={runsQuery.hasNextPage}
          loading={runsQuery.isPending}
          loadingMore={runsQuery.isFetchingNextPage}
          error={runsQuery.isError}
          notice={runNotice}
        />
      ) : (
        <Button type="button" variant="secondary" onClick={() => setShowRunPicker(true)}>
          叠加已授权运行
        </Button>
      )}
      {activeAddedIds.map((id, index) => (
        <AddedRunBlock
          key={id}
          runId={id}
          query={addedSessionQueries[index]}
          selection={addedSelections[id] ?? null}
          authority={authority}
          onToggle={(track, data) => toggleAddedTrack(id, data, track)}
          onRemove={() => removeRun(id)}
        />
      ))}
      {localPlaced.length > 0 ? (
        <ul className="space-y-2">
          {localPlaced.map((item) => {
            const rejected = !item.placement.ok || item.foreign;
            const reason = item.foreign
              ? `${item.track.label} names another website. It stays in this browser and is not added to the window.`
              : item.placement.ok
                ? item.placement.warning
                : item.placement.reason;
            return (
              <li
                key={item.track.id}
                data-testid={rejected ? 'reference-conflict' : 'local-browser-track'}
                data-source="local"
                data-track-id={item.track.id}
                data-assembly={item.track.assembly}
                role="status"
              >
                <span className="font-medium">{item.track.label}</span>
                <span className="block">
                  Source local file · {item.track.assembly}
                  {reason ? ` · ${reason}` : ''}
                </span>
              </li>
            );
          })}
        </ul>
      ) : null}
      <LocusSearch assembly={authority.scaffold} region={region} onNavigate={onNavigate} />
      <RegionNotes region={region} rows={platformRows} />
      <FailureList rows={failedRows} failures={loadFailures} />
      <Viewer
        scaffold={authority.scaffold}
        region={region}
        rows={visibleRows}
        localTracks={hubLocal}
        storeId={`helixweave-browser-${browserKey}`}
        showToolbar={showToolbar}
      />
    </section>
  );
}

const EMPTY_LOCAL: LocalWashuTrack[] = [];

function RunPicker({
  runs,
  pendingRunId,
  onPending,
  onAdd,
  onMore,
  hasMore,
  loading,
  loadingMore,
  error,
  notice,
}: {
  runs: RunSummaryResponse[];
  pendingRunId: string;
  onPending: (runId: string) => void;
  onAdd: () => void;
  onMore: () => void;
  hasMore: boolean;
  loading: boolean;
  loadingMore: boolean;
  error: boolean;
  notice: string | null;
}) {
  return (
    <div className="space-y-2">
      <label className="block text-sm">
        Authorized runs
        <select
          className="mt-1 block w-full max-w-xl border border-[var(--color-border)] px-2 py-1"
          aria-label="Authorized runs"
          value={pendingRunId}
          onChange={(event) => onPending(event.target.value)}
        >
          <option value="">Select a run</option>
          {runs.map((run) => (
            <option key={run.run_id} value={run.run_id}>
              {run.workflow_id} · {run.run_id} · {run.status} · {run.updated_at}
            </option>
          ))}
        </select>
      </label>
      <div className="flex flex-wrap gap-2">
        <Button type="button" variant="secondary" onClick={onAdd} disabled={pendingRunId === ''}>
          Add run
        </Button>
        {hasMore ? (
          <Button type="button" variant="secondary" onClick={onMore} disabled={loadingMore}>
            More authorized runs
          </Button>
        ) : null}
      </div>
      {loading ? <p role="status">Loading authorized runs…</p> : null}
      {error ? (
        <p role="status">Authorized runs could not be loaded. Tracks already open stay available.</p>
      ) : null}
      {!loading && !error && runs.length === 0 ? (
        <p role="status">No authorized runs are available.</p>
      ) : null}
      {notice ? <p role="status">{notice}</p> : null}
    </div>
  );
}

function AddedRunBlock({
  runId,
  query,
  selection,
  authority,
  onToggle,
  onRemove,
}: {
  runId: string;
  query:
    | {
        data?: BrowserSessionResponse;
        error: unknown;
        isError: boolean;
        isPending: boolean;
        refetch: () => Promise<unknown>;
      }
    | undefined;
  selection: RunSelection | null;
  authority: WindowAuthority;
  onToggle: (track: BrowserTrackResponse, session: BrowserSessionResponse) => void;
  onRemove: () => void;
}) {
  const queryClient = useQueryClient();
  const data = query?.data?.run_id === runId ? query.data : undefined;
  const conflict = isGenerationConflict(query?.error);
  const prepare = useMutation({
    mutationFn: async (track: BrowserTrackResponse) => {
      if (!data?.artifact_generation) {
        throw new ApiError(409, 'BROWSER_TRACK_READ_CONFLICT', 'Reload the browser view.');
      }
      return prepareRunBrowserTrack(runId, track.artifact_id, {
        generation: data.artifact_generation,
        revision: track.revision,
      });
    },
    onSuccess: async () => {
      await queryClient.invalidateQueries({ queryKey: ['browser-session', runId] });
    },
  });
  const placement =
    data && data.ok && !conflict
      ? placePlatformRun(authority, runId, referenceBinding(data.reference))
      : null;
  const chosen = data && placement?.ok ? activeChosen(selection, data) : [];
  const selected = chosen.map((item) => item.artifactId);
  const stale =
    data && selection && selection.generation === (data.artifact_generation ?? null)
      ? selection.stale
      : [];

  return (
    <article
      className="space-y-2"
      aria-label={`Run ${runId}`}
      data-testid="browser-run"
      data-run-id={runId}
      data-run-source="platform"
      data-generation={data?.artifact_generation ?? ''}
      data-reference-identity={data?.reference.identity_sha256 ?? ''}
      data-scaffold={data?.reference.viewer_scaffold ?? ''}
    >
      <div className="flex flex-wrap items-center justify-between gap-2">
        <p>
          {data
            ? `Source ${data.workflow_id} · run ${data.run_id} · generation ${data.artifact_generation ?? 'unknown'}${data.reference.revision_id ? ` · reference revision ${data.reference.revision_id}` : ''}`
            : `Run ${runId}`}
        </p>
        <Button type="button" variant="secondary" onClick={onRemove}>
          Remove run {runId}
        </Button>
      </div>
      {query?.isPending ? <p role="status">Loading run {runId}…</p> : null}
      {conflict ? (
        <p role="status">
          These files changed. Refresh the list and choose the tracks again. New files are not
          read until you choose them.
        </p>
      ) : null}
      {query?.isError && !conflict ? (
        <p role="status">
          {messageFromError(query.error)} Tracks from other runs stay available.
        </p>
      ) : null}
      {data && !data.ok ? (
        <p role="status">This run's browser view could not be loaded. Other tracks stay available.</p>
      ) : null}
      {placement && !placement.ok ? (
        <p role="status" data-testid="reference-conflict">
          {placement.reason}
        </p>
      ) : null}
      {placement?.ok && placement.warning ? (
        <p role="status" data-testid="reference-note">
          {placement.warning}
        </p>
      ) : null}
      {data && data.ok && !conflict ? (
        <TrackChooser
          runId={runId}
          generation={data.artifact_generation ?? null}
          tracks={data.tracks}
          unavailable={data.unavailable}
          selected={placement?.ok ? selected : []}
          stale={stale}
          selectable={placement?.ok === true}
          onToggle={(track) => onToggle(track, data)}
        />
      ) : null}
      {data && data.ok ? (
        <Button type="button" variant="secondary" onClick={() => void query?.refetch()}>
          Refresh run {runId}
        </Button>
      ) : query?.isError ? (
        <Button type="button" variant="secondary" onClick={() => void query.refetch()}>
          Retry run {runId}
        </Button>
      ) : null}
      <div className="flex flex-wrap gap-2">
        {data && placement?.ok
          ? data.tracks
              .filter(
                (track) => selected.includes(track.artifact_id) && track.preparation === 'needed',
              )
              .map((track) => (
                <Button
                  key={track.artifact_id}
                  type="button"
                  variant="secondary"
                  onClick={() => prepare.mutate(track)}
                  disabled={prepare.isPending}
                >
                  Prepare {track.label}
                </Button>
              ))
          : null}
      </div>
      {prepare.isError ? <p role="status">{messageFromError(prepare.error)}</p> : null}
    </article>
  );
}

function ReferenceNotice({ session }: { session: BrowserSessionResponse }) {
  const reference = session.reference;
  if (reference.status === 'unbound') {
    return (
      <p role="status">
        This run has no reference profile binding, so tracks are not placed on a
        guessed assembly.
      </p>
    );
  }
  if (reference.status !== 'bound' || !reference.viewer_scaffold) {
    return (
      <p role="status">
        The reference identity is present, but this WashU build has no chromosome
        scaffold for assembly {reference.assembly_label ?? 'unknown'}.
      </p>
    );
  }
  return (
    <p role="status">
      Reference {reference.profile_id} revision {reference.revision_id}, identity{' '}
      {reference.identity_sha256}. Navigation uses the built-in {reference.viewer_scaffold}{' '}
      chromosome list and does not prove the reference files match that public genome.
    </p>
  );
}

function InteractionNotice({ session }: { session: BrowserSessionResponse }) {
  if (session.interaction_status === 'available') return null;
  return (
    <p role="status">
      This run has no interaction artifact. Signal or peak tracks are not a substitute
      for Hi-TrAC PET or matrix results.
    </p>
  );
}

function TrackChooser({
  runId,
  generation,
  tracks,
  unavailable,
  selected,
  stale,
  selectable,
  onToggle,
}: {
  runId: string;
  generation: string | null;
  tracks: BrowserTrackResponse[];
  unavailable: BrowserSessionResponse['unavailable'];
  selected: string[];
  stale: string[];
  selectable: boolean;
  onToggle: (track: BrowserTrackResponse) => void;
}) {
  if (tracks.length === 0 && unavailable.length === 0) {
    return <p role="status">This run has no artifacts the genome browser can display.</p>;
  }
  return (
    <ul className="space-y-2">
      {tracks.map((track) => (
        <li
          key={track.artifact_id}
          data-testid="browser-track"
          data-run-id={runId}
          data-artifact-id={track.artifact_id}
          data-revision={track.revision}
          data-generation={generation ?? ''}
          data-source="platform"
        >
          <label className="flex items-start gap-2 text-sm">
            <input
              type="checkbox"
              checked={selectable && selected.includes(track.artifact_id)}
              onChange={() => {
                if (selectable && trackCanBeChosen(track)) onToggle(track);
              }}
              disabled={!selectable || !trackCanBeChosen(track)}
            />
            <span>
              <span className="font-medium">{track.label}</span>
              <span className="block text-[var(--color-text-muted)]">{track.detail}</span>
              <span className="block">
                Source {runId} · revision {track.revision}
              </span>
              {track.preparation === 'needed' ? (
                <span className="block" role="status">
                  An index is not prepared yet. Choose this file, then prepare it.
                </span>
              ) : null}
              {track.preparation === 'source_too_large' ? (
                <span className="block" role="status">
                  This file is too large to index here and is not shown as empty data.
                </span>
              ) : null}
              {track.preparation === 'reference_mismatch' ? (
                <span className="block" role="status">
                  This artifact does not match the run reference identity.
                </span>
              ) : null}
              {stale.includes(track.artifact_id) ? (
                <span className="block" role="status">
                  This file changed. Choose it again before the new version is read.
                </span>
              ) : null}
            </span>
          </label>
        </li>
      ))}
      {unavailable.map((item) => (
        <li key={item.artifact_id}>
          <label className="flex items-start gap-2 text-sm">
            <input type="checkbox" checked={false} disabled onChange={() => undefined} />
            <span>
              <span className="font-medium">{item.name}</span>
              <span className="block" role="status">
                Unsupported format: {item.name}. {item.detail}
              </span>
            </span>
          </label>
        </li>
      ))}
    </ul>
  );
}

function RegionNotes({ region, rows }: { region: string; rows: PlatformRow[] }) {
  const seen = seenWindow(region);
  const countable = rows.filter(
    (row) =>
      row.generation &&
      row.track.kind !== 'signal' &&
      row.track.kind !== 'contact_matrix',
  );
  const summary = useQuery({
    queryKey: [
      'browser-region',
      region,
      countable.map(
        (row) =>
          `${row.runId}:${row.generation}:${row.track.artifact_id}:${row.track.revision}`,
      ),
    ],
    enabled: Boolean(seen && countable.length > 0),
    queryFn: async ({ signal }) => {
      if (!seen) return [];
      return Promise.allSettled(
        countable.map((row) =>
          summarizeRunBrowserTrackRegion(
            row.runId,
            row.track.artifact_id,
            row.generation ?? '',
            row.track.revision,
            {
              chrom: seen.chrom,
              start: seen.start,
              end: seen.end,
            },
            { signal },
          ),
        ),
      );
    },
  });
  if (!summary.data) return null;
  return (
    <div>
      {summary.data.map((result, index) => {
        const row = countable[index];
        if (!row) return null;
        if (result.status === 'rejected') {
          return (
            <p key={`${row.runId}:${row.track.artifact_id}`} role="status">
              {row.track.label}: This window count failed. The track stays available.
            </p>
          );
        }
        const item = result.value;
        return (
          <p key={`${row.runId}:${row.track.artifact_id}`} role="status">
            {row.track.label}: {item.detail}
            {item.records != null ? ` (${item.records} records)` : ''}
            {item.records != null && item.records > 100000 ? ` ${DRAW_NOTE}` : ''}
          </p>
        );
      })}
    </div>
  );
}

function FailureList({
  rows,
  failures,
}: {
  rows: PlatformRow[];
  failures: Record<string, string>;
}) {
  if (rows.length === 0) return null;
  return (
    <div role="alert">
      {rows.map((row) => (
        <p key={trackReadKey(row.runId, row.track.artifact_id)}>
          {row.track.label} ({row.runId}): {failures[trackReadKey(row.runId, row.track.artifact_id)]}
        </p>
      ))}
    </div>
  );
}

async function trackReadError(response: Response): Promise<Error> {
  const status = `Track failed to load (HTTP ${response.status}).`;
  if (response.status === 401 || response.status === 403) {
    return new Error(`${status} Check your session and access to this run.`);
  }
  if (response.status === 404) return new Error(`${status} This artifact is unavailable.`);
  if (response.status === 409) {
    const body = await response.json().catch(() => null);
    const code = body?.issues?.[0]?.context?.reason_code;
    if (['generation_changed', 'revision_changed', 'derivative_changed'].includes(code)) {
      return new Error(`${status} The selected version changed (${code}). Reload the run and select the current track.`);
    }
  }
  return new Error(status);
}

function Viewer({
  scaffold,
  region,
  rows,
  localTracks,
  storeId,
  showToolbar,
}: {
  scaffold: string;
  region: string;
  rows: PlatformRow[];
  localTracks: Array<{ id: string; model: WashuTrackConfig }>;
  storeId: string;
  showToolbar: boolean;
}) {
  const dataTracks = useMemo(() => {
    const models: WashuTrackConfig[] = [{ type: 'ruler', name: 'Ruler' }];
    const names = new Map<string, number>();
    const take = (model: WashuTrackConfig) => {
      const name = typeof model.name === 'string' ? model.name : '';
      const seen = names.get(name) ?? 0;
      names.set(name, seen + 1);
      models.push(
        seen === 0 || name === '' ? model : { ...model, name: `${name} (${seen + 1})` },
      );
    };
    for (const row of rows) {
      take({
        type: row.track.wash_u_type,
        name: row.track.label,
        url: sameOriginPath(row.track.content_path) ?? undefined,
      });
    }
    for (const track of localTracks) take(track.model);
    // 62.0.2 uses this flag for a complete initial track selection. Without
    // it the vendor's track-list effect can dispatch updates back to itself.
    return models.map(model => ({ ...model, waitToUpdate: true }));
  }, [localTracks, rows]);
  const trackKey = [
    scaffold,
    ...rows.map(row => `${row.runId}:${row.generation}:${row.track.artifact_id}:${row.track.revision}`),
    ...localTracks.map((track) => `local:${track.id}`),
  ].join('|');
  const visible = Boolean(scaffold && REGION.test(region) && dataTracks.length > 1);
  // Adding a bigWig to an already mounted hub crashes that hub, so the key
  // includes the chosen files. GenomeHub also draws the region it receives on
  // mount and drops a later viewRegion change while its track-load guard is
  // closed. A coordinate typed in that moment moved the field and left the
  // ruler behind, so the key includes the region and the new hub mounts on the
  // window it should draw. A zoom stays inside the current hub and does not
  // change this key. The store id stays on this browser: the vendor cache
  // never drops a store, and a new id on each coordinate would keep the old
  // tracks for the life of the page.
  const hubKey = visible ? JSON.stringify([trackKey, region]) : '';
  const storeConfig = useMemo(
    () => ({ storeId, enablePersistence: false }),
    [storeId],
  );
  if (!visible) return null;
  return (
    <div className="min-h-[28rem] overflow-hidden" data-testid="washu-genome-viewer">
      <GenomeHub
        key={hubKey}
        storeConfig={storeConfig}
        genomeName={scaffold}
        viewRegion={region}
        tracks={dataTracks}
        showGenomeNavigator={false}
        showNavBar={false}
        showToolBar={showToolbar}
        showDisclosure={false}
      />
    </div>
  );
}
