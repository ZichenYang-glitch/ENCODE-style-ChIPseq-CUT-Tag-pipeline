import { useEffect, useMemo, useRef, useState } from 'react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { GenomeHub } from 'wuepgg';
import washUStyles from 'wuepgg/style.css?raw';
import {
  getRunBrowserSession,
  prepareRunBrowserTrack,
  summarizeRunBrowserTrackRegion,
} from '../../api/generated/artifacts/artifacts';
import type {
  BrowserSessionResponse,
  BrowserTrackResponse,
} from '../../api/generated/models';

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

function reconcileSelection(
  current: TrackSelection,
  session: BrowserSessionResponse,
): TrackSelection {
  const generation = session.artifact_generation ?? null;
  if (current.generation !== generation) {
    return { runId: current.runId, generation, chosen: [], stale: [] };
  }
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
  if (sameChosen && sameStale) return current;
  return { ...current, chosen, stale };
}
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import {
  beginTrackRead,
  clearTrackReadFailures,
  ensureTrackReadGuard,
  noteTrackReadSucceeded,
  subscribeTrackReadFailure,
} from './trackReadGuard';

const REGION = /^([A-Za-z0-9_.:-]+):([0-9]+)-([0-9]+)$/;
const DRAW_NOTE =
  'WashU draws at most 100000 rows in one window. A larger count is not fully drawn.';

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

export default function GenomeBrowserPanel({
  runId,
  initialRegion = null,
  artifactIds = null,
  showToolbar = true,
  onReady,
  onFailed,
}: {
  runId: string;
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
  const sessionQuery = useQuery({
    queryKey: ['browser-session', runId],
    queryFn: ({ signal }) => getRunBrowserSession(runId, { signal }),
    // A generation conflict is the user's prompt to refresh. Retrying would
    // load the new catalog by itself and hide that choice.
    retry: (failureCount, error) =>
      !isGenerationConflict(error) && failureCount < 3,
    // Tests observe the exhausted retries immediately. Production keeps the backoff.
    retryDelay: import.meta.env.MODE === 'test' ? 0 : undefined,
  });
  const session =
    sessionQuery.data?.run_id === runId ? sessionQuery.data : undefined;
  const [selection, setSelection] = useState<TrackSelection | null>(null);
  const [region, setRegion] = useState('');
  const [loadFailures, setLoadFailures] = useState<Record<string, string>>({});
  const selectionRef = useRef<TrackSelection | null>(null);
  const pendingToggles = useRef<BrowserTrackResponse[]>([]);
  const presetRef = useRef<readonly string[] | null>(artifactIds ?? null);
  const conflictRef = useRef<string | null>(null);
  const regionKey = useRef<string | null>(null);
  selectionRef.current = selection;
  presetRef.current = artifactIds ?? null;
  if (conflictRef.current !== null && conflictRef.current !== runId) {
    conflictRef.current = null;
  }
  if (isGenerationConflict(sessionQuery.error)) conflictRef.current = runId;

  const onReadyRef = useRef(onReady);
  const onFailedRef = useRef(onFailed);
  onReadyRef.current = onReady;
  onFailedRef.current = onFailed;

  useEffect(
    () =>
      subscribeTrackReadFailure((artifactId) => {
        setLoadFailures((current) =>
          current[artifactId]
            ? current
            : { ...current, [artifactId]: 'Track failed to load.' },
        );
      }),
    [],
  );
  useEffect(() => () => clearTrackReadFailures(), []);

  useEffect(() => {
    const style = document.createElement('style');
    style.dataset.helixweaveWashu = runId;
    style.textContent = washUStyles;
    document.head.appendChild(style);
    return () => {
      style.remove();
      localStorage.removeItem(`persist:helixweave-browser-${runId}`);
    };
  }, [runId]);

  useEffect(() => {
    if (sessionQuery.isError && !isGenerationConflict(sessionQuery.error)) {
      onFailedRef.current?.(messageFromError(sessionQuery.error));
    }
  }, [sessionQuery.error, sessionQuery.isError]);

  useEffect(() => {
    if (!session || isGenerationConflict(sessionQuery.error)) return;
    const generation = session.artifact_generation ?? null;
    const current = selectionRef.current;
    let next: TrackSelection;
    if (current === null || current.runId !== runId) {
      const skipPreset = conflictRef.current === runId;
      if (skipPreset) conflictRef.current = null;
      next = {
        runId,
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
  }, [runId, session, sessionQuery.error]);

  useEffect(() => {
    if (!session) return;
    const key = `${runId}:${session.artifact_generation ?? ''}:${initialRegion ?? ''}`;
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
    // viewer query the intended window.
    const parts = target.match(/^(.*:)(\d+)-(\d+)$/);
    const hop = parts
      ? `${parts[1]}${Math.max(0, Number(parts[2]) - 400)}-${Math.max(1, Number(parts[3]) - 400)}`
      : target;
    setRegion(target);
    const move = window.setTimeout(() => setRegion(hop), 600);
    const settle = window.setTimeout(() => setRegion(target), 1100);
    return () => {
      window.clearTimeout(move);
      window.clearTimeout(settle);
    };
  }, [initialRegion, runId, session]);

  const generationConflict = isGenerationConflict(sessionQuery.error);
  const activeChoice = useMemo(() => {
    if (
      generationConflict ||
      session === undefined ||
      selection === null ||
      selection.runId !== runId ||
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
  }, [generationConflict, runId, selection, session]);
  const selected = useMemo(
    () => activeChoice.map((item) => item.artifactId),
    [activeChoice],
  );
  const stale =
    session !== undefined &&
    selection !== null &&
    selection.runId === runId &&
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
          current !== null && current.runId === runId
            ? current
            : { runId, generation, chosen: [], stale: [] };
        for (const item of batch) {
          const exists = next.chosen.some(
            (chosen) => chosen.artifactId === item.artifact_id,
          );
          next = {
            ...next,
            chosen: exists
              ? next.chosen.filter((chosen) => chosen.artifactId !== item.artifact_id)
              : [
                  ...next.chosen,
                  { artifactId: item.artifact_id, revision: item.revision },
                ],
            stale: next.stale.filter((artifactId) => artifactId !== item.artifact_id),
          };
        }
        selectionRef.current = next;
        return next;
      });
    });
  }

  const readyTracks = useMemo(
    () =>
      (session?.tracks ?? []).filter(
        (track) =>
          selected.includes(track.artifact_id) &&
          track.preparation === 'ready' &&
          sameOriginPath(track.content_path) !== null,
      ),
    [selected, session],
  );
  const visibleTracks = useMemo(
    () => readyTracks.filter((track) => !loadFailures[track.artifact_id]),
    [loadFailures, readyTracks],
  );

  useEffect(() => {
    const controllers = new Map<string, AbortController>();
    let cancelled = false;
    async function preflight(track: BrowserTrackResponse) {
      const path = sameOriginPath(track.content_path);
      if (!path) return;
      const controller = new AbortController();
      controllers.set(track.artifact_id, controller);
      const generation = beginTrackRead(track.artifact_id, controller.signal);
      try {
        const response = await fetch(path, {
          method: 'GET',
          headers: { Range: 'bytes=0-0' },
          credentials: 'same-origin',
          signal: controller.signal,
        });
        if (!response.ok && response.status !== 206) {
          throw new Error(`HTTP ${response.status}`);
        }
        const indexPath = sameOriginPath(track.index_path);
        if (indexPath) {
          const indexResponse = await fetch(indexPath, {
            method: 'GET',
            headers: { Range: 'bytes=0-1' },
            credentials: 'same-origin',
            signal: controller.signal,
          });
          if (!indexResponse.ok && indexResponse.status !== 206) {
            throw new Error(`HTTP ${indexResponse.status}`);
          }
        }
        if (!cancelled && noteTrackReadSucceeded(track.artifact_id, generation)) {
          setLoadFailures((current) => {
            if (!current[track.artifact_id]) return current;
            const next = { ...current };
            delete next[track.artifact_id];
            return next;
          });
        }
      } catch (error) {
        if (controller.signal.aborted || cancelled) return;
        if (error instanceof DOMException && error.name === 'AbortError') return;
        setLoadFailures((current) => ({
          ...current,
          [track.artifact_id]:
            error instanceof Error && error.message.startsWith('HTTP ')
              ? `Track failed to load (${error.message}).`
              : 'Track failed to load.',
        }));
      }
    }
    for (const track of readyTracks) void preflight(track);
    return () => {
      cancelled = true;
      for (const controller of controllers.values()) controller.abort();
    };
  }, [readyTracks]);

  const prepare = useMutation({
    mutationFn: async (track: BrowserTrackResponse) => {
      if (!session?.artifact_generation) {
        throw new ApiError(409, 'BROWSER_TRACK_READ_CONFLICT', 'Reload the browser view.');
      }
      return prepareRunBrowserTrack(runId, track.artifact_id, {
        generation: session.artifact_generation,
        revision: track.revision,
      });
    },
    onSuccess: async () => {
      await queryClient.invalidateQueries({ queryKey: ['browser-session', runId] });
    },
  });

  if (generationConflict) {
    return (
      <div role="alert">
        <p>
          These files changed. Refresh the list and choose the tracks again. New
          files are not read until you choose them.
        </p>
        <Button variant="secondary" onClick={() => void sessionQuery.refetch()}>
          Refresh list
        </Button>
      </div>
    );
  }
  if (sessionQuery.isError || (session !== undefined && !session.ok)) {
    return (
      <div role="alert">
        <p>{messageFromError(sessionQuery.error)}</p>
        <Button variant="secondary" onClick={() => void sessionQuery.refetch()}>
          Retry browser view
        </Button>
      </div>
    );
  }
  if (sessionQuery.isPending || session === undefined) {
    return <p role="status">Loading genome browser…</p>;
  }

  return (
    <section className="min-w-0 space-y-3" aria-label="Genome browser">
      <ReferenceNotice session={session} />
      <InteractionNotice session={session} />
      <TrackChooser
        tracks={session.tracks}
        unavailable={session.unavailable}
        selected={selected}
        stale={stale}
        onToggle={toggleTrack}
      />
      <Button variant="secondary" onClick={() => void sessionQuery.refetch()}>
        Refresh list
      </Button>
      <label className="block text-sm">
        Genomic region
        <input
          className="mt-1 w-full max-w-md border border-[var(--color-border)] px-2 py-1"
          value={region}
          onChange={(event) => setRegion(event.target.value.trim())}
          aria-label="Genomic region"
        />
      </label>
      <div className="flex flex-wrap gap-2">
        {session.tracks
          .filter(
            (track) =>
              selected.includes(track.artifact_id) && track.preparation === 'needed',
          )
          .map((track) => (
            <Button
              key={track.artifact_id}
              variant="secondary"
              onClick={() => prepare.mutate(track)}
              disabled={prepare.isPending}
            >
              Prepare {track.label}
            </Button>
          ))}
      </div>
      {prepare.isError ? <p role="alert">{messageFromError(prepare.error)}</p> : null}
      <RegionNotes
        runId={runId}
        generation={session.artifact_generation}
        region={region}
        tracks={readyTracks}
      />
      <FailureList tracks={readyTracks} failures={loadFailures} />
      <Viewer
        session={session}
        region={region}
        tracks={visibleTracks}
        trackKey={readyTracks
          .map((track) => `${track.artifact_id}:${track.revision}`)
          .join('|')}
        showToolbar={showToolbar}
      />
    </section>
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
  tracks,
  unavailable,
  selected,
  stale,
  onToggle,
}: {
  tracks: BrowserTrackResponse[];
  unavailable: BrowserSessionResponse['unavailable'];
  selected: string[];
  stale: string[];
  onToggle: (track: BrowserTrackResponse) => void;
}) {
  if (tracks.length === 0 && unavailable.length === 0) {
    return <p role="status">This run has no artifacts the genome browser can display.</p>;
  }
  return (
    <ul className="space-y-2">
      {tracks.map((track) => (
        <li key={track.artifact_id}>
          <label className="flex items-start gap-2 text-sm">
            <input
              type="checkbox"
              checked={selected.includes(track.artifact_id)}
              onChange={() => onToggle(track)}
              disabled={!trackCanBeChosen(track)}
            />
            <span>
              <span className="font-medium">{track.label}</span>
              <span className="block text-[var(--color-text-muted)]">{track.detail}</span>
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

function RegionNotes({
  runId,
  generation,
  region,
  tracks,
}: {
  runId: string;
  generation: string | null | undefined;
  region: string;
  tracks: BrowserTrackResponse[];
}) {
  const parsed = REGION.exec(region);
  const countable = tracks.filter(
    (track) => track.kind !== 'signal' && track.kind !== 'contact_matrix',
  );
  const summary = useQuery({
    queryKey: [
      'browser-region',
      runId,
      region,
      countable.map((track) => track.artifact_id),
    ],
    enabled: Boolean(generation && parsed && countable.length > 0),
    queryFn: async ({ signal }) => {
      if (!parsed || !generation) return [];
      return Promise.all(
        countable.map((track) =>
          summarizeRunBrowserTrackRegion(
            runId,
            track.artifact_id,
            generation,
            track.revision,
            {
              chrom: parsed[1],
              start: Number(parsed[2]),
              end: Number(parsed[3]),
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
      {summary.data.map((item, index) => (
        <p key={countable[index]?.artifact_id ?? index} role="status">
          {countable[index]?.label}: {item.detail}
          {item.records != null ? ` (${item.records} records)` : ''}
          {item.records != null && item.records > 100000 ? ` ${DRAW_NOTE}` : ''}
        </p>
      ))}
    </div>
  );
}

function FailureList({
  tracks,
  failures,
}: {
  tracks: BrowserTrackResponse[];
  failures: Record<string, string>;
}) {
  const failed = tracks.filter((track) => failures[track.artifact_id]);
  if (failed.length === 0) return null;
  return (
    <div role="alert">
      {failed.map((track) => (
        <p key={track.artifact_id}>
          {track.label}: {failures[track.artifact_id]}
        </p>
      ))}
    </div>
  );
}

function Viewer({
  session,
  region,
  tracks,
  trackKey,
  showToolbar,
}: {
  session: BrowserSessionResponse;
  region: string;
  tracks: BrowserTrackResponse[];
  trackKey: string;
  showToolbar: boolean;
}) {
  const scaffold = session.reference.viewer_scaffold;
  const model = useMemo(
    () => [
      { type: 'ruler', name: 'Ruler' },
      ...tracks.map((track) => ({
        type: track.wash_u_type,
        name: track.label,
        url: sameOriginPath(track.content_path) ?? undefined,
      })),
    ],
    [tracks],
  );
  const storeConfig = useMemo(
    () => ({
      storeId: `helixweave-browser-${session.run_id}`,
      enablePersistence: false,
    }),
    [session.run_id],
  );
  if (!scaffold || !REGION.test(region) || tracks.length === 0) return null;
  const crossOrigin = tracks.some((track) => sameOriginPath(track.content_path) === null);
  if (crossOrigin) {
    return (
      <p role="alert">
        Genome tracks must be read from this site. A separate API origin cannot send
        the session cookie.
      </p>
    );
  }
  return (
    <div className="min-h-[28rem] overflow-hidden" data-testid="washu-genome-viewer">
      {/* Adding a bigWig to an already mounted hub crashes that hub. Remount
          when the chosen files change and let the new hub take the region. */}
      <GenomeHub
        key={trackKey}
        storeConfig={storeConfig}
        genomeName={scaffold}
        viewRegion={region}
        tracks={model}
        showGenomeNavigator={false}
        showNavBar={false}
        showToolBar={showToolbar}
        showDisclosure={false}
      />
    </div>
  );
}
