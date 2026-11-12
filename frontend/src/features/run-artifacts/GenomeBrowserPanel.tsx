import { useEffect, useMemo, useState } from 'react';
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
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';

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

export default function GenomeBrowserPanel({ runId }: { runId: string }) {
  const queryClient = useQueryClient();
  const sessionQuery = useQuery({
    queryKey: ['browser-session', runId],
    queryFn: ({ signal }) => getRunBrowserSession(runId, { signal }),
  });
  const session = sessionQuery.data;
  const [selected, setSelected] = useState<string[]>([]);
  const [region, setRegion] = useState('');
  const [loadFailures, setLoadFailures] = useState<Record<string, string>>({});

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
    if (!session) return;
    const initial = session.tracks
      .filter((track) => track.preparation === 'ready' || track.preparation === 'needed')
      .map((track) => track.artifact_id);
    setSelected(initial);
    const suggested = session.tracks.find(
      (track) => track.provenance?.suggested_region,
    )?.provenance?.suggested_region;
    if (typeof suggested === 'string') setRegion(suggested);
  }, [session]);

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

  useEffect(() => {
    const controller = new AbortController();
    let cancelled = false;
    async function preflight(track: BrowserTrackResponse) {
      const path = sameOriginPath(track.content_path);
      if (!path) return;
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
        if (!cancelled) {
          setLoadFailures((current) => {
            const next = { ...current };
            delete next[track.artifact_id];
            return next;
          });
        }
      } catch (error) {
        if (controller.signal.aborted || cancelled) return;
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
      controller.abort();
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

  if (sessionQuery.isPending) {
    return <p role="status">Loading genome browser…</p>;
  }
  if (sessionQuery.isError || !session?.ok) {
    return (
      <div role="alert">
        <p>{messageFromError(sessionQuery.error)}</p>
        <Button variant="secondary" onClick={() => void sessionQuery.refetch()}>
          Retry browser view
        </Button>
      </div>
    );
  }

  return (
    <section className="min-w-0 space-y-3" aria-label="Genome browser">
      <ReferenceNotice session={session} />
      <InteractionNotice session={session} />
      <UnavailableList session={session} />
      <TrackChooser
        tracks={session.tracks}
        selected={selected}
        onToggle={(artifactId) =>
          setSelected((current) =>
            current.includes(artifactId)
              ? current.filter((item) => item !== artifactId)
              : [...current, artifactId],
          )
        }
      />
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
        tracks={readyTracks.filter((track) => !loadFailures[track.artifact_id])}
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

function UnavailableList({ session }: { session: BrowserSessionResponse }) {
  if (session.unavailable.length === 0) return null;
  return (
    <ul>
      {session.unavailable.map((item) => (
        <li key={item.artifact_id} role="status">
          Unsupported format: {item.name}. {item.detail}
        </li>
      ))}
    </ul>
  );
}

function TrackChooser({
  tracks,
  selected,
  onToggle,
}: {
  tracks: BrowserTrackResponse[];
  selected: string[];
  onToggle: (artifactId: string) => void;
}) {
  if (tracks.length === 0) {
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
              onChange={() => onToggle(track.artifact_id)}
              disabled={
                track.preparation === 'reference_mismatch' ||
                track.preparation === 'source_too_large'
              }
            />
            <span>
              <span className="font-medium">{track.label}</span>
              <span className="block text-[var(--color-text-muted)]">{track.detail}</span>
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
}: {
  session: BrowserSessionResponse;
  region: string;
  tracks: BrowserTrackResponse[];
}) {
  const scaffold = session.reference.viewer_scaffold;
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
      <GenomeHub
        storeConfig={{
          storeId: `helixweave-browser-${session.run_id}`,
          enablePersistence: false,
        }}
        genomeName={scaffold}
        viewRegion={region}
        tracks={[
          { type: 'ruler', name: 'Ruler' },
          ...tracks.map((track) => ({
            type: track.wash_u_type,
            name: track.label,
            url: sameOriginPath(track.content_path) ?? undefined,
          })),
        ]}
        showGenomeNavigator={false}
        showNavBar={false}
        showToolBar={true}
        showDisclosure={false}
      />
    </div>
  );
}
