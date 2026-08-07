import { ExternalLink } from 'lucide-react';
import { Link } from 'react-router-dom';
import type { ArtifactPublicationResponse } from '../../api/generated/models';
import { Button } from '../../components/Button';
import { derivedArtifactGenerationStatus } from './artifactPublicationState';

interface ArtifactPublicationListProps {
  publications: ArtifactPublicationResponse[];
  hasNextPage: boolean;
  isFetchingNextPage: boolean;
  onLoadMore: () => void;
}

const DATE_TIME_FORMAT = new Intl.DateTimeFormat(undefined, {
  dateStyle: 'medium',
  timeStyle: 'short',
});

function detailPath(publication: ArtifactPublicationResponse): string {
  return `/artifacts/${encodeURIComponent(publication.run_id)}/${encodeURIComponent(publication.artifact_id)}?generation=${encodeURIComponent(publication.artifact_generation)}`;
}

function formatTimestamp(value: string): string {
  return DATE_TIME_FORMAT.format(new Date(value));
}

function PublicationLink({ publication }: { publication: ArtifactPublicationResponse }) {
  return (
    <Link
      className="inline-flex max-w-full items-center gap-1 rounded font-medium text-[var(--color-accent-hover)] underline-offset-2 hover:underline focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
      to={detailPath(publication)}
      title={publication.artifact_id}
    >
      <span className="min-w-0 break-all">{publication.artifact_id}</span>
      <ExternalLink className="shrink-0" aria-hidden="true" size={13} />
    </Link>
  );
}

function GenerationState({ publication }: { publication: ArtifactPublicationResponse }) {
  const status = derivedArtifactGenerationStatus(publication);
  return (
    <span
      className={`inline-flex rounded-full border px-2 py-0.5 text-[11px] font-semibold ${
        status === 'current'
          ? 'border-emerald-300 bg-emerald-50 text-emerald-800'
          : 'border-slate-300 bg-slate-100 text-slate-700'
      }`}
    >
      {status === 'current' ? 'Current' : 'Superseded'}
    </span>
  );
}

function AssociatedRunSampleSummary({
  publication,
}: {
  publication: ArtifactPublicationResponse;
}) {
  const samples = publication.run_sample_binding.associated_run_samples;
  if (samples.length === 0) {
    return (
      <span className="text-[var(--color-text-muted)]">
        No associated run samples were recorded for this run binding.
      </span>
    );
  }
  const visibleSamples = samples.slice(0, 2);
  const remaining = samples.length - visibleSamples.length;
  return (
    <ul className="space-y-1">
      {visibleSamples.map((sample) => (
        <li key={sample.sample_revision_id} className="min-w-0">
          <span className="block break-all" title={sample.sample_id}>
            {sample.sample_id}
          </span>
          <span
            className="block break-all font-mono text-[11px] text-[var(--color-text-muted)]"
            title={sample.sample_revision_id}
          >
            {sample.sample_revision_id}
          </span>
        </li>
      ))}
      {remaining > 0 && (
        <li className="text-[11px] text-[var(--color-text-muted)]">
          +{remaining} more associated run samples
        </li>
      )}
    </ul>
  );
}

export function ArtifactPublicationList({
  publications,
  hasNextPage,
  isFetchingNextPage,
  onLoadMore,
}: ArtifactPublicationListProps) {
  return (
    <div className="min-w-0" data-testid="artifact-publication-list">
      <p className="mb-2 text-xs text-[var(--color-text-muted)]">
        These are input samples associated with the run; they are not per-artifact sample attribution.
      </p>
      <div className="hidden min-w-0 lg:block">
        <table
          className="w-full table-fixed border-collapse text-left text-xs"
          aria-label="Artifact publications"
        >
          <thead>
            <tr className="border-y border-[var(--color-border)] text-[var(--color-text-muted)]">
              <th className="w-[11%] px-2 py-2 font-semibold">Project</th>
              <th className="w-[11%] px-2 py-2 font-semibold">Run</th>
              <th className="w-[11%] px-2 py-2 font-semibold">Workflow</th>
              <th className="w-[13%] px-2 py-2 font-semibold">Artifact type</th>
              <th className="w-[12%] px-2 py-2 font-semibold">Published</th>
              <th className="w-[18%] px-2 py-2 font-semibold">Generation / revision</th>
              <th className="w-[16%] px-2 py-2 font-semibold">Associated run samples</th>
              <th className="w-[8%] px-2 py-2 font-semibold">State</th>
            </tr>
          </thead>
          <tbody>
            {publications.map((publication) => (
              <tr
                key={`${publication.run_id}:${publication.artifact_id}:${publication.artifact_generation}`}
                className="border-b border-[var(--color-border)] align-top hover:bg-[var(--color-bg)]"
              >
                <td className="break-all px-2 py-2" title={publication.project_id}>
                  {publication.project_id}
                </td>
                <td className="break-all px-2 py-2" title={publication.run_id}>
                  {publication.run_id}
                </td>
                <td className="break-all px-2 py-2" title={publication.workflow_id}>
                  {publication.workflow_id}
                </td>
                <td className="min-w-0 px-2 py-2">
                  <span className="block break-all text-[var(--color-text-muted)]">
                    {publication.output_type}
                  </span>
                  <PublicationLink publication={publication} />
                </td>
                <td className="px-2 py-2 text-[var(--color-text-muted)]">
                  <time dateTime={publication.published_at}>
                    {formatTimestamp(publication.published_at)}
                  </time>
                </td>
                <td className="px-2 py-2 font-mono text-[11px] text-[var(--color-text-muted)]">
                  <span className="block break-all" title={publication.artifact_generation}>
                    {publication.artifact_generation}
                  </span>
                  <span className="mt-1 block break-all" title={publication.artifact_revision}>
                    {publication.artifact_revision}
                  </span>
                  {derivedArtifactGenerationStatus(publication) === 'superseded' && (
                    <span className="mt-1 block break-all" title={publication.current_artifact_generation}>
                      Current generation: {publication.current_artifact_generation}
                    </span>
                  )}
                </td>
                <td className="px-2 py-2">
                  <AssociatedRunSampleSummary publication={publication} />
                </td>
                <td className="px-2 py-2"><GenerationState publication={publication} /></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <ul
        className="divide-y divide-[var(--color-border)] border-y border-[var(--color-border)] lg:hidden"
        aria-label="Artifact publications"
      >
        {publications.map((publication) => (
          <li
            key={`${publication.run_id}:${publication.artifact_id}:${publication.artifact_generation}`}
            className="min-w-0 py-3"
          >
            <div className="flex min-w-0 flex-wrap items-start justify-between gap-2">
              <div className="min-w-0">
                <p className="break-all text-sm font-medium">{publication.output_type}</p>
                <PublicationLink publication={publication} />
              </div>
              <GenerationState publication={publication} />
            </div>
            <dl className="mt-2 grid min-w-0 grid-cols-[6rem_minmax(0,1fr)] gap-x-2 gap-y-1 text-xs">
              <dt className="text-[var(--color-text-muted)]">Project</dt>
              <dd className="break-all">{publication.project_id}</dd>
              <dt className="text-[var(--color-text-muted)]">Run</dt>
              <dd className="break-all">{publication.run_id}</dd>
              <dt className="text-[var(--color-text-muted)]">Workflow</dt>
              <dd className="break-all">{publication.workflow_id}</dd>
              <dt className="text-[var(--color-text-muted)]">Published</dt>
              <dd>
                <time dateTime={publication.published_at}>
                  {formatTimestamp(publication.published_at)}
                </time>
              </dd>
              <dt className="text-[var(--color-text-muted)]">Generation</dt>
              <dd className="break-all font-mono text-[11px]">{publication.artifact_generation}</dd>
              <dt className="text-[var(--color-text-muted)]">Revision</dt>
              <dd className="break-all font-mono text-[11px]">{publication.artifact_revision}</dd>
              {derivedArtifactGenerationStatus(publication) === 'superseded' && (
                <>
                  <dt className="text-[var(--color-text-muted)]">Current generation</dt>
                  <dd className="break-all font-mono text-[11px]">
                    {publication.current_artifact_generation}
                  </dd>
                </>
              )}
              <dt className="text-[var(--color-text-muted)]">Associated run samples</dt>
              <dd><AssociatedRunSampleSummary publication={publication} /></dd>
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
            aria-label="Load more artifact publications"
          >
            {isFetchingNextPage ? 'Loading more…' : 'Load more'}
          </Button>
        </div>
      )}
    </div>
  );
}
