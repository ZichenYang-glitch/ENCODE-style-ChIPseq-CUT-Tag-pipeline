import { useState } from 'react';
import { useQuery } from '@tanstack/react-query';
import { catalogListArchivedTracks } from '../../api/archivedTracks';
import { Button } from '../../components/Button';
import { ArchivedTrackFacts } from './ArchivedTrackFacts';
import {
  ARCHIVE_PAGE_LIMIT,
  archivedTrackDataProps,
  archivedTrackOpenName,
  archivedViewerHref,
  archiveErrorFromUnknown,
} from './archiveView';

export function ArchivedTrackList({
  projectId,
  experimentId,
}: {
  projectId: string;
  experimentId: string;
}) {
  const scope = `${projectId}:${experimentId}`;
  const [scopeKey, setScopeKey] = useState(scope);
  const [offset, setOffset] = useState(0);
  const [previousOffsets, setPreviousOffsets] = useState<number[]>([]);
  if (scopeKey !== scope) {
    setScopeKey(scope);
    setOffset(0);
    setPreviousOffsets([]);
  }

  const query = useQuery({
    queryKey: ['archived-tracks', projectId, experimentId, offset],
    retry: false,
    queryFn: ({ signal }) =>
      catalogListArchivedTracks(
        projectId,
        experimentId,
        { limit: ARCHIVE_PAGE_LIMIT, offset },
        { signal },
      ),
  });

  const error = query.isError ? archiveErrorFromUnknown(query.error) : null;
  const page = query.data ?? null;
  const nextOffset = page?.next_offset ?? null;

  return (
    <section
      aria-label="已封存轨道"
      lang="zh-CN"
      className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
      data-testid="archived-track-list"
      data-project-id={projectId}
      data-experiment-id={experimentId}
    >
      <div>
        <h3 className="text-lg font-semibold">已封存轨道</h3>
        <p className="text-sm text-[var(--color-text-muted)]">
          选择一个已封存的结果版本。页面不会改打开较新的同名轨道。
        </p>
      </div>
      {query.isPending ? <p role="status">正在读取已封存版本…</p> : null}
      {error !== null ? (
        <p role="alert" data-testid="archived-track-error" data-error-code={error.code}>
          {error.message}
        </p>
      ) : null}
      {page !== null && page.items.length === 0 ? (
        <p role="status">这条实验没有已封存轨道。</p>
      ) : null}
      {page !== null && page.items.length > 0 ? (
        <ul aria-label="已封存版本" className="space-y-3">
          {page.items.map((track) => (
            <li
              key={track.archive_version_id}
              className="space-y-2 rounded-[4px] border border-[var(--color-border)] p-3"
              data-testid="archived-track-version"
              {...archivedTrackDataProps(track)}
            >
              <ArchivedTrackFacts track={track} />
              <a
                className="inline-block text-sm text-[var(--color-link)] underline"
                data-testid="archived-track-open"
                href={archivedViewerHref(track)}
                aria-label={archivedTrackOpenName(track)}
              >
                在 WashU 查看
              </a>
            </li>
          ))}
        </ul>
      ) : null}
      {previousOffsets.length > 0 || nextOffset !== null ? (
        <div className="flex flex-wrap gap-2">
          <Button
            type="button"
            disabled={previousOffsets.length === 0 || query.isFetching}
            onClick={() => {
              const previous = previousOffsets[previousOffsets.length - 1];
              if (previous === undefined) return;
              setPreviousOffsets((current) => current.slice(0, -1));
              setOffset(previous);
            }}
          >
            上一页
          </Button>
          <Button
            type="button"
            disabled={nextOffset === null || query.isFetching}
            onClick={() => {
              if (nextOffset === null) return;
              setPreviousOffsets((current) => [...current, offset]);
              setOffset(nextOffset);
            }}
          >
            下一页
          </Button>
        </div>
      ) : null}
    </section>
  );
}
