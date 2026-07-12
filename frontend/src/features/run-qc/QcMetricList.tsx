import { ExternalLink } from 'lucide-react';
import type { QcMetricResponse } from '../../api/generated/models';
import { Button } from '../../components/Button';
import { formatQcProducedTime } from './qcState';

interface QcMetricListProps {
  metrics: QcMetricResponse[];
  hasNextPage: boolean;
  isFetchingNextPage: boolean;
  onLoadMore: () => void;
  onOpenSourceArtifact: (artifactId: string) => void;
}

function optionalText(value: string | null | undefined): string {
  return value && value.trim() ? value : '—';
}

function sampleExperiment(metric: QcMetricResponse): string {
  return [metric.sample_id, metric.experiment_id]
    .filter((value): value is string => Boolean(value))
    .join(' · ') || '—';
}

function QcFlag({ value }: { value: QcMetricResponse['qc_flag'] }) {
  if (value === null) {
    return <span className="text-[var(--color-text-muted)]">Not reported</span>;
  }
  const tone = {
    pass: 'border-[var(--color-accent)] text-[var(--color-accent-hover)]',
    warning:
      'border-[var(--color-warning)] bg-[var(--color-warning-bg)] text-[var(--color-warning)]',
    fail: 'border-[var(--color-error)] bg-[var(--color-error-bg)] text-[var(--color-error)]',
  }[value];
  return (
    <span className={`inline-flex rounded border px-1.5 py-0.5 font-medium ${tone}`}>
      {value}
    </span>
  );
}

function MetricIdentity({ metric }: { metric: QcMetricResponse }) {
  return (
    <div className="min-w-0">
      <span className="block break-words font-medium text-[var(--color-text)]">
        {metric.display_name}
      </span>
      <code
        className="block break-all text-[11px] text-[var(--color-text-muted)] [overflow-wrap:anywhere]"
        title={metric.metric_key}
      >
        {metric.metric_key}
      </code>
    </div>
  );
}

function SourceArtifactAction({
  metric,
  onOpen,
}: {
  metric: QcMetricResponse;
  onOpen: (artifactId: string) => void;
}) {
  return (
    <div className="min-w-0">
      <code
        className="block break-all text-[11px] text-[var(--color-text-muted)]"
        title={metric.source_artifact_id}
      >
        {metric.source_artifact_id}
      </code>
      <Button
        className="mt-1 px-2 py-1 text-xs"
        variant="secondary"
        onClick={() => onOpen(metric.source_artifact_id)}
        aria-label={`Open source artifact for ${metric.display_name}`}
        title="Open source artifact"
      >
        <ExternalLink className="mr-1 h-3.5 w-3.5" aria-hidden="true" />
        View source
      </Button>
    </div>
  );
}

export function QcMetricList({
  metrics,
  hasNextPage,
  isFetchingNextPage,
  onLoadMore,
  onOpenSourceArtifact,
}: QcMetricListProps) {
  return (
    <div className="min-w-0" data-testid="qc-metric-list">
      <div className="hidden min-w-0 md:block">
        <table
          className="w-full table-fixed border-collapse text-left text-xs"
          aria-label="Indexed run QC metrics"
        >
          <thead>
            <tr className="border-y border-[var(--color-border)] text-[var(--color-text-muted)]">
              <th className="w-[22%] px-2 py-2 font-semibold">Metric</th>
              <th className="w-[13%] px-2 py-2 font-semibold">Value / unit</th>
              <th className="w-[12%] px-2 py-2 font-semibold">Scope</th>
              <th className="w-[18%] px-2 py-2 font-semibold">Sample / experiment</th>
              <th className="w-[13%] px-2 py-2 font-semibold">Assay / flag</th>
              <th className="w-[22%] px-2 py-2 font-semibold">Source / produced</th>
            </tr>
          </thead>
          <tbody>
            {metrics.map((metric) => (
              <tr
                key={metric.metric_id}
                className="border-b border-[var(--color-border)] align-top hover:bg-[var(--color-bg)]"
              >
                <td className="min-w-0 px-2 py-2">
                  <MetricIdentity metric={metric} />
                </td>
                <td className="min-w-0 px-2 py-2">
                  <code className="block break-all font-semibold" title={metric.value}>
                    {metric.value}
                  </code>
                  <span className="mt-1 block text-[var(--color-text-muted)]">
                    {metric.unit}
                  </span>
                </td>
                <td className="break-words px-2 py-2 text-[var(--color-text-muted)]">
                  {metric.scope}
                </td>
                <td className="break-words px-2 py-2 text-[var(--color-text-muted)]">
                  {sampleExperiment(metric)}
                </td>
                <td className="min-w-0 px-2 py-2">
                  <span className="block break-words text-[var(--color-text-muted)]">
                    {optionalText(metric.assay)}
                  </span>
                  <span className="mt-1 block"><QcFlag value={metric.qc_flag} /></span>
                </td>
                <td className="min-w-0 px-2 py-2">
                  <SourceArtifactAction
                    metric={metric}
                    onOpen={onOpenSourceArtifact}
                  />
                  <time
                    className="mt-2 block text-[var(--color-text-muted)]"
                    dateTime={metric.produced_at}
                  >
                    {formatQcProducedTime(metric.produced_at)}
                  </time>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <ul
        className="divide-y divide-[var(--color-border)] border-y border-[var(--color-border)] md:hidden"
        aria-label="Indexed run QC metrics"
      >
        {metrics.map((metric) => (
          <li key={metric.metric_id} className="min-w-0 py-3">
            <MetricIdentity metric={metric} />
            <dl className="mt-2 grid min-w-0 grid-cols-[6.5rem_minmax(0,1fr)] gap-x-2 gap-y-1 text-xs">
              <dt className="text-[var(--color-text-muted)]">Value</dt>
              <dd className="min-w-0">
                <code className="break-all font-semibold" title={metric.value}>
                  {metric.value}
                </code>{' '}
                <span className="text-[var(--color-text-muted)]">{metric.unit}</span>
              </dd>
              <dt className="text-[var(--color-text-muted)]">Scope</dt>
              <dd className="break-words">{metric.scope}</dd>
              <dt className="text-[var(--color-text-muted)]">Sample / exp.</dt>
              <dd className="break-words">{sampleExperiment(metric)}</dd>
              <dt className="text-[var(--color-text-muted)]">Assay</dt>
              <dd className="break-words">{optionalText(metric.assay)}</dd>
              <dt className="text-[var(--color-text-muted)]">QC flag</dt>
              <dd><QcFlag value={metric.qc_flag} /></dd>
              <dt className="text-[var(--color-text-muted)]">Source</dt>
              <dd className="min-w-0">
                <SourceArtifactAction
                  metric={metric}
                  onOpen={onOpenSourceArtifact}
                />
              </dd>
              <dt className="text-[var(--color-text-muted)]">Produced</dt>
              <dd>
                <time dateTime={metric.produced_at}>
                  {formatQcProducedTime(metric.produced_at)}
                </time>
              </dd>
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
            aria-label="Load more QC metrics"
          >
            {isFetchingNextPage ? 'Loading more…' : 'Load more'}
          </Button>
        </div>
      )}
    </div>
  );
}
