import type { ArtifactReferenceResponse } from '../../api/generated/models';
import { Button } from '../../components/Button';
import { formatBytes, formatProducedTime } from './artifactState';

interface ArtifactListProps {
  artifacts: ArtifactReferenceResponse[];
  selectedArtifactId: string | null;
  onSelect: (artifactId: string) => void;
  hasNextPage: boolean;
  isFetchingNextPage: boolean;
  onLoadMore: () => void;
}

function optionalText(value: string | null | undefined): string {
  return value && value.trim() ? value : '—';
}

function sampleExperiment(artifact: ArtifactReferenceResponse): string {
  return [artifact.metadata.sample_id, artifact.metadata.experiment_id]
    .filter((value): value is string => Boolean(value))
    .join(' · ') || '—';
}

function ArtifactName({ artifact }: { artifact: ArtifactReferenceResponse }) {
  return (
    <div className="min-w-0">
      <span className="block break-words font-medium text-[var(--color-text)]">
        {artifact.output_type}
      </span>
      <span
        className="block break-words text-[var(--color-text-muted)] [overflow-wrap:anywhere]"
        title={artifact.name}
      >
        {artifact.name}
      </span>
    </div>
  );
}

export function ArtifactList({
  artifacts,
  selectedArtifactId,
  onSelect,
  hasNextPage,
  isFetchingNextPage,
  onLoadMore,
}: ArtifactListProps) {
  return (
    <div id="artifact-list" className="min-w-0" data-testid="artifact-list">
      <div className="hidden min-w-0 md:block">
        <table className="w-full table-fixed border-collapse text-left text-xs">
          <caption className="sr-only">Indexed run artifacts</caption>
          <thead>
            <tr className="border-y border-[var(--color-border)] text-[var(--color-text-muted)]">
              <th className="w-[18%] px-2 py-2 font-semibold">Artifact</th>
              <th className="w-[27%] px-2 py-2 font-semibold">Relative path</th>
              <th className="w-[9%] px-2 py-2 font-semibold">Size</th>
              <th className="w-[10%] px-2 py-2 font-semibold">Scope</th>
              <th className="w-[18%] px-2 py-2 font-semibold">Sample / experiment</th>
              <th className="w-[18%] px-2 py-2 font-semibold">Assay / produced</th>
            </tr>
          </thead>
          <tbody>
            {artifacts.map((artifact) => {
              const selected = artifact.artifact_id === selectedArtifactId;
              return (
                <tr
                  key={artifact.artifact_id}
                  className={`border-b border-[var(--color-border)] align-top ${
                    selected ? 'bg-teal-50' : 'hover:bg-[var(--color-bg)]'
                  }`}
                  data-selected={selected || undefined}
                >
                  <td className="min-w-0 px-2 py-2">
                    <button
                      type="button"
                      className="min-w-0 text-left focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
                      onClick={() => onSelect(artifact.artifact_id)}
                      aria-label={`Open artifact ${artifact.name}`}
                      aria-pressed={selected}
                    >
                      <ArtifactName artifact={artifact} />
                    </button>
                  </td>
                  <td className="min-w-0 px-2 py-2">
                    <code
                      className="block break-all text-[11px] text-[var(--color-text-muted)]"
                      title={artifact.relative_path}
                    >
                      {artifact.relative_path}
                    </code>
                  </td>
                  <td className="px-2 py-2 text-[var(--color-text-muted)]">
                    {formatBytes(artifact.size_bytes)}
                  </td>
                  <td className="break-words px-2 py-2 text-[var(--color-text-muted)]">
                    {optionalText(artifact.metadata.scope)}
                  </td>
                  <td className="break-words px-2 py-2 text-[var(--color-text-muted)]">
                    {sampleExperiment(artifact)}
                  </td>
                  <td className="break-words px-2 py-2 text-[var(--color-text-muted)]">
                    <span className="block">{optionalText(artifact.metadata.assay)}</span>
                    <time className="mt-1 block" dateTime={artifact.produced_at}>
                      {formatProducedTime(artifact.produced_at)}
                    </time>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <ul className="divide-y divide-[var(--color-border)] border-y border-[var(--color-border)] md:hidden">
        {artifacts.map((artifact) => {
          const selected = artifact.artifact_id === selectedArtifactId;
          return (
            <li
              key={artifact.artifact_id}
              className={`min-w-0 py-3 ${selected ? 'border-l-2 border-l-[var(--color-accent)] pl-2' : ''}`}
              data-selected={selected || undefined}
            >
              <button
                type="button"
                className="w-full min-w-0 text-left focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
                onClick={() => onSelect(artifact.artifact_id)}
                aria-label={`Open artifact ${artifact.name}`}
                aria-pressed={selected}
              >
                <ArtifactName artifact={artifact} />
              </button>
              <dl className="mt-2 grid min-w-0 grid-cols-[6rem_minmax(0,1fr)] gap-x-2 gap-y-1 text-xs">
                <dt className="text-[var(--color-text-muted)]">Path</dt>
                <dd className="min-w-0">
                  <code className="break-all text-[11px]" title={artifact.relative_path}>
                    {artifact.relative_path}
                  </code>
                </dd>
                <dt className="text-[var(--color-text-muted)]">Size</dt>
                <dd>{formatBytes(artifact.size_bytes)}</dd>
                <dt className="text-[var(--color-text-muted)]">Scope</dt>
                <dd className="break-words">{optionalText(artifact.metadata.scope)}</dd>
                <dt className="text-[var(--color-text-muted)]">Sample / exp.</dt>
                <dd className="break-words">{sampleExperiment(artifact)}</dd>
                <dt className="text-[var(--color-text-muted)]">Assay</dt>
                <dd className="break-words">{optionalText(artifact.metadata.assay)}</dd>
                <dt className="text-[var(--color-text-muted)]">Produced</dt>
                <dd>
                  <time dateTime={artifact.produced_at}>
                    {formatProducedTime(artifact.produced_at)}
                  </time>
                </dd>
              </dl>
            </li>
          );
        })}
      </ul>

      {hasNextPage && (
        <div className="pt-3">
          <Button
            variant="secondary"
            onClick={onLoadMore}
            disabled={isFetchingNextPage}
            aria-label="Load more artifacts"
          >
            {isFetchingNextPage ? 'Loading more…' : 'Load more'}
          </Button>
        </div>
      )}
    </div>
  );
}
