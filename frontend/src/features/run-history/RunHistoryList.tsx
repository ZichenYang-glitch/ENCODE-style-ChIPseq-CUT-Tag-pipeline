import { ExternalLink, FileArchive, Gauge } from 'lucide-react';
import { Link } from 'react-router-dom';
import type { RunSummaryResponse } from '../../api/generated/models';
import { Button } from '../../components/Button';
import { RunStatusBadge } from '../run-progress/RunStatusBadge';

interface RunHistoryListProps {
  runs: RunSummaryResponse[];
  hasNextPage: boolean;
  isFetchingNextPage: boolean;
  onLoadMore: () => void;
}
const DATE_TIME_FORMAT = new Intl.DateTimeFormat(undefined, {
  dateStyle: 'medium',
  timeStyle: 'short',
});

function formatTimestamp(value: string): string {
  return DATE_TIME_FORMAT.format(new Date(value));
}

function formatDuration(run: RunSummaryResponse): string {
  if (run.started_at === null) return 'Not started';
  const end = run.ended_at ?? run.updated_at;
  const seconds = Math.max(
    0,
    Math.round((Date.parse(end) - Date.parse(run.started_at)) / 1000),
  );
  if (seconds < 60) return `${seconds}s`;
  const minutes = Math.floor(seconds / 60);
  const remainder = seconds % 60;
  if (minutes < 60) return `${minutes}m ${remainder}s`;
  const hours = Math.floor(minutes / 60);
  return `${hours}h ${minutes % 60}m`;
}

function ResultLinks({ run }: { run: RunSummaryResponse }) {
  if (run.status !== 'succeeded') return null;
  const root = `/runs/${encodeURIComponent(run.run_id)}`;
  return (
    <span className="flex flex-wrap items-center gap-1.5">
      <Link
        className="inline-flex items-center gap-1 rounded text-xs font-medium text-[var(--color-accent-hover)] underline-offset-2 hover:underline focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
        to={`${root}?view=qc`}
        aria-label={`Open QC for run ${run.run_id}`}
      >
        <Gauge aria-hidden="true" size={13} />
        QC
      </Link>
      <Link
        className="inline-flex items-center gap-1 rounded text-xs font-medium text-[var(--color-accent-hover)] underline-offset-2 hover:underline focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
        to={`${root}?view=artifacts`}
        aria-label={`Open artifacts for run ${run.run_id}`}
      >
        <FileArchive aria-hidden="true" size={13} />
        Artifacts
      </Link>
    </span>
  );
}

function RunIdentity({ run }: { run: RunSummaryResponse }) {
  return (
    <div className="min-w-0">
      <Link
        className="inline-flex max-w-full items-center gap-1 rounded font-medium text-[var(--color-accent-hover)] underline-offset-2 hover:underline focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
        to={`/runs/${encodeURIComponent(run.run_id)}`}
        title={run.run_id}
      >
        <span className="min-w-0 break-all">{run.run_id}</span>
        <ExternalLink className="shrink-0" aria-hidden="true" size={13} />
      </Link>
      <span
        className="mt-1 block break-all text-[11px] text-[var(--color-text-muted)]"
        title={run.workflow_id}
      >
        {run.workflow_id}
      </span>
    </div>
  );
}

export function RunHistoryList({
  runs,
  hasNextPage,
  isFetchingNextPage,
  onLoadMore,
}: RunHistoryListProps) {
  return (
    <div className="min-w-0" data-testid="run-history-list">
      <div className="hidden min-w-0 md:block">
        <table
          className="w-full table-fixed border-collapse text-left text-xs"
          aria-label="Run history"
        >
          <thead>
            <tr className="border-y border-[var(--color-border)] text-[var(--color-text-muted)]">
              <th className="w-[29%] px-2 py-2 font-semibold">Run / workflow</th>
              <th className="w-[13%] px-2 py-2 font-semibold">Status</th>
              <th className="w-[20%] px-2 py-2 font-semibold">Created</th>
              <th className="w-[14%] px-2 py-2 font-semibold">Duration</th>
              <th className="w-[12%] px-2 py-2 font-semibold">Stage</th>
              <th className="w-[12%] px-2 py-2 font-semibold">Results</th>
            </tr>
          </thead>
          <tbody>
            {runs.map((run) => (
              <tr
                key={run.run_id}
                className="border-b border-[var(--color-border)] align-top hover:bg-[var(--color-bg)]"
              >
                <td className="min-w-0 px-2 py-2"><RunIdentity run={run} /></td>
                <td className="px-2 py-2"><RunStatusBadge status={run.status} /></td>
                <td className="px-2 py-2 text-[var(--color-text-muted)]">
                  <time dateTime={run.created_at}>{formatTimestamp(run.created_at)}</time>
                </td>
                <td className="px-2 py-2 text-[var(--color-text-muted)]">
                  {formatDuration(run)}
                  {run.ended_at && (
                    <time className="mt-1 block" dateTime={run.ended_at}>
                      Ended {formatTimestamp(run.ended_at)}
                    </time>
                  )}
                </td>
                <td className="break-words px-2 py-2 text-[var(--color-text-muted)]">
                  {run.current_stage ?? '—'}
                </td>
                <td className="px-2 py-2"><ResultLinks run={run} /></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <ul
        className="divide-y divide-[var(--color-border)] border-y border-[var(--color-border)] md:hidden"
        aria-label="Run history"
      >
        {runs.map((run) => (
          <li key={run.run_id} className="min-w-0 py-3">
            <div className="flex min-w-0 flex-wrap items-start justify-between gap-2">
              <RunIdentity run={run} />
              <RunStatusBadge status={run.status} />
            </div>
            <dl className="mt-2 grid min-w-0 grid-cols-[5.5rem_minmax(0,1fr)] gap-x-2 gap-y-1 text-xs">
              <dt className="text-[var(--color-text-muted)]">Created</dt>
              <dd><time dateTime={run.created_at}>{formatTimestamp(run.created_at)}</time></dd>
              <dt className="text-[var(--color-text-muted)]">Duration</dt>
              <dd>{formatDuration(run)}</dd>
              <dt className="text-[var(--color-text-muted)]">Stage</dt>
              <dd className="break-words">{run.current_stage ?? '—'}</dd>
              {run.status === 'succeeded' && (
                <>
                  <dt className="text-[var(--color-text-muted)]">Results</dt>
                  <dd><ResultLinks run={run} /></dd>
                </>
              )}
            </dl>
          </li>
        ))}
      </ul>

      {hasNextPage && (
        <div className="pt-3">
          <Button
            variant="secondary"
            onClick={onLoadMore}
            disabled={isFetchingNextPage}
            aria-label="Load more runs"
          >
            {isFetchingNextPage ? 'Loading more…' : 'Load more'}
          </Button>
        </div>
      )}
    </div>
  );
}
