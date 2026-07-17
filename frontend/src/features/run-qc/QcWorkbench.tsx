import { useInfiniteQuery, useQueryClient } from '@tanstack/react-query';
import { RefreshCw, RotateCcw } from 'lucide-react';
import { useEffect, useRef, type ReactNode } from 'react';
import type {
  QcMetricResponse,
  RunQcMetricsResponse,
} from '../../api/generated/models';
import { listRunQcMetrics } from '../../api/generated/qc-metrics/qc-metrics';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import { isValidArtifactId } from '../run-artifacts/artifactState';
import { QcMetricList } from './QcMetricList';
import {
  flattenQcMetricPages,
  isValidQcGeneration,
  isValidQcMetricId,
  qcPaginationAnomaly,
  safeNextQcCursor,
  type QcIndexingOutcome,
} from './qcState';

const QC_PAGE_SIZE = 50;
const DECIMAL_TEXT_PATTERN = /^-?(?:0|[1-9]\d{0,25})(?:\.\d{1,12})?$/;
const METRIC_KEY_PATTERN = /^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$/;
const ALLOWED_UNITS = new Set(['count', 'fraction', 'ratio', 'score']);
const ALLOWED_SCOPES = new Set(['run', 'sample', 'experiment']);
const ALLOWED_FLAGS = new Set(['pass', 'warning', 'fail']);

interface QcWorkbenchProps {
  runId: string;
  runStatus: string;
  outcome: QcIndexingOutcome;
  onOpenSourceArtifact: (artifactId: string) => void;
  onRefreshStatus: () => void;
}

export function qcMetricListQueryKey(
  runId: string,
  generation: string | null,
) {
  return ['run-qc-metrics', runId, generation] as const;
}

function isOptionalSafeText(value: unknown): value is string | null {
  return (
    value === null ||
    (typeof value === 'string' &&
      value.length > 0 &&
      value.length <= 255 &&
      !/[\u0000-\u001f\u007f]/.test(value))
  );
}

function isValidMetric(value: unknown): value is QcMetricResponse {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) {
    return false;
  }
  const metric = value as Record<string, unknown>;
  return (
    isValidQcMetricId(
      typeof metric.metric_id === 'string' ? metric.metric_id : null,
    ) &&
    typeof metric.metric_key === 'string' &&
    METRIC_KEY_PATTERN.test(metric.metric_key) &&
    typeof metric.display_name === 'string' &&
    metric.display_name.trim().length > 0 &&
    metric.display_name.length <= 255 &&
    !/[\u0000-\u001f\u007f]/.test(metric.display_name) &&
    typeof metric.value === 'string' &&
    DECIMAL_TEXT_PATTERN.test(metric.value) &&
    typeof metric.unit === 'string' &&
    ALLOWED_UNITS.has(metric.unit) &&
    typeof metric.scope === 'string' &&
    ALLOWED_SCOPES.has(metric.scope) &&
    isOptionalSafeText(metric.sample_id) &&
    isOptionalSafeText(metric.experiment_id) &&
    isOptionalSafeText(metric.assay) &&
    (metric.qc_flag === null ||
      (typeof metric.qc_flag === 'string' && ALLOWED_FLAGS.has(metric.qc_flag))) &&
    typeof metric.source_artifact_id === 'string' &&
    isValidArtifactId(metric.source_artifact_id) &&
    typeof metric.produced_at === 'string' &&
    !Number.isNaN(new Date(metric.produced_at).getTime())
  );
}

function validateListResponse(
  response: RunQcMetricsResponse,
  runId: string,
  generation: string,
): RunQcMetricsResponse {
  const responseGeneration = response.qc_generation;
  if (
    response.ok !== true ||
    response.run_id !== runId ||
    !isValidQcGeneration(responseGeneration) ||
    responseGeneration !== generation ||
    !Array.isArray(response.qc_metrics) ||
    response.qc_metrics.some((metric) => !isValidMetric(metric))
  ) {
    throw new Error('invalid QC metric list response');
  }
  return response;
}

function isCursorNotFound(error: unknown): boolean {
  return (
    error instanceof ApiError && error.code === 'RUN_QC_METRIC_CURSOR_NOT_FOUND'
  );
}

function isGenerationChanged(error: unknown): boolean {
  return (
    error instanceof ApiError &&
    error.code === 'RUN_QC_METRIC_GENERATION_CHANGED'
  );
}

export function QcWorkbench({
  runId,
  runStatus,
  outcome,
  onOpenSourceArtifact,
  onRefreshStatus,
}: QcWorkbenchProps) {
  const queryClient = useQueryClient();
  const indexed = runStatus === 'succeeded' && outcome.kind === 'indexed';
  const generation = outcome.kind === 'indexed' ? outcome.generation : null;
  const queryKey = qcMetricListQueryKey(runId, generation);
  const listQuery = useInfiniteQuery({
    queryKey,
    enabled: indexed,
    initialPageParam: undefined as string | undefined,
    queryFn: async ({ pageParam }) => {
      if (!generation) throw new Error('QC generation is unavailable');
      const parameters = {
        after: pageParam,
        generation,
        limit: QC_PAGE_SIZE,
      };
      return validateListResponse(
        await listRunQcMetrics(runId, parameters),
        runId,
        generation,
      );
    },
    getNextPageParam: (lastPage, allPages, _lastPageParam, allPageParams) =>
      safeNextQcCursor(lastPage, allPages, allPageParams),
    retry: 1,
    retryDelay: 100,
  });
  const generationChangeHandled = useRef<string | null>(null);
  const generationChanged = isGenerationChanged(listQuery.error);

  useEffect(() => {
    if (!generation) {
      queryClient.removeQueries({ queryKey: ['run-qc-metrics', runId] });
      return;
    }
    queryClient.removeQueries({
      queryKey: ['run-qc-metrics', runId],
      predicate: (query) => query.queryKey[2] !== generation,
    });
  }, [generation, queryClient, runId]);

  useEffect(() => {
    if (!generationChanged || !generation) return;
    if (generationChangeHandled.current === generation) return;
    generationChangeHandled.current = generation;
    queryClient.removeQueries({ queryKey, exact: true });
    onRefreshStatus();
  }, [generation, generationChanged, onRefreshStatus, queryClient, queryKey]);

  if (runStatus !== 'succeeded') {
    return (
      <QcStatus title="QC metrics are not available yet">
        QC metrics become available after the run completes successfully.
      </QcStatus>
    );
  }

  if (outcome.kind === 'pending') {
    return (
      <QcStatus title="Indexing QC metrics" busy>
        The workflow succeeded. Waiting for the worker to finish indexing persisted QC metrics.
      </QcStatus>
    );
  }

  if (outcome.kind === 'unconfirmed') {
    return (
      <QcStatus title="QC metric status could not be confirmed">
        The available run history is incomplete or malformed. Refresh the canonical run status before continuing.
        <StatusRefreshButton onClick={onRefreshStatus} />
      </QcStatus>
    );
  }

  if (outcome.kind === 'failed') {
    return (
      <QcStatus title="QC metric indexing failed" tone="error">
        The run remains succeeded, but its QC metric summary could not be indexed.
        {outcome.reasonCode && (
          <span className="mt-1 block font-mono text-xs">
            Reference: {outcome.reasonCode}
          </span>
        )}
        <StatusRefreshButton onClick={onRefreshStatus} />
      </QcStatus>
    );
  }

  const pages = generationChanged ? [] : (listQuery.data?.pages ?? []);
  const pageParams = generationChanged ? [] : (listQuery.data?.pageParams ?? []);
  let metrics: QcMetricResponse[] = [];
  let pageIntegrityFailed = false;
  try {
    metrics = flattenQcMetricPages(pages);
  } catch {
    pageIntegrityFailed = true;
  }
  const hasCachedMetrics = metrics.length > 0;
  const cursorAnomaly = qcPaginationAnomaly(pages, pageParams);
  const cursorWasReplaced = isCursorNotFound(listQuery.error);
  const cursorRecoveryNeeded =
    cursorAnomaly || cursorWasReplaced || pageIntegrityFailed;
  const queryFailed = listQuery.error !== null;

  function reloadFromFirstPage() {
    void queryClient.resetQueries({ queryKey, exact: true });
  }

  if (listQuery.isPending) return <QcListSkeleton />;

  if (queryFailed && !hasCachedMetrics) {
    return (
      <QcStatus title="QC metrics could not be loaded" tone="error">
        The persisted QC metric index is temporarily unavailable. Existing run state was not changed.
        <Button
          className="mt-3"
          variant="secondary"
          onClick={() => void listQuery.refetch()}
        >
          Retry QC metrics
        </Button>
      </QcStatus>
    );
  }

  if (outcome.count === 0 && metrics.length === 0 && !queryFailed) {
    return (
      <QcStatus title="No indexed QC metrics">
        QC metric indexing completed successfully and produced an empty result set.
      </QcStatus>
    );
  }

  if (metrics.length === 0) {
    return (
      <QcStatus title="QC metric index is not ready">
        The indexing event and persisted QC metric list do not yet agree. Retry before treating this run as empty.
        <Button
          className="mt-3"
          variant="secondary"
          onClick={() => void listQuery.refetch()}
        >
          Retry QC metrics
        </Button>
      </QcStatus>
    );
  }

  return (
    <div className="min-w-0 space-y-3" data-testid="qc-workbench">
      <div className="flex min-w-0 flex-wrap items-baseline justify-between gap-2">
        <div>
          <h3 className="text-sm font-semibold">Indexed QC metrics</h3>
          <p className="text-xs text-[var(--color-text-muted)]">
            {metrics.length} loaded
            {outcome.count >= metrics.length ? ` of ${outcome.count}` : ''}
          </p>
        </div>
      </div>

      {cursorRecoveryNeeded && (
        <div
          className="flex flex-wrap items-center justify-between gap-2 border-y border-[var(--color-warning)] bg-[var(--color-warning-bg)] px-2 py-2 text-sm text-[var(--color-warning)]"
          role="status"
        >
          <span>
            {cursorWasReplaced
              ? 'QC results changed while pages were being loaded. Reload from the first page.'
              : 'QC metric pagination stopped because the next cursor was not safe. Reload from the first page.'}
          </span>
          <Button
            variant="secondary"
            onClick={reloadFromFirstPage}
            aria-label="Reload QC metrics from first page"
          >
            <RotateCcw className="mr-1.5 h-4 w-4" aria-hidden="true" />
            Reload first page
          </Button>
        </div>
      )}

      {queryFailed && hasCachedMetrics && !cursorRecoveryNeeded && (
        <div
          className="flex flex-wrap items-center justify-between gap-2 border-y border-[var(--color-warning)] bg-[var(--color-warning-bg)] px-2 py-2 text-sm text-[var(--color-warning)]"
          role="status"
        >
          <span>
            Additional QC data could not be loaded. Cached QC rows are preserved and may be stale.
          </span>
          <Button variant="secondary" onClick={() => void listQuery.refetch()}>
            Retry QC metrics
          </Button>
        </div>
      )}

      {listQuery.isFetching &&
        !listQuery.isFetchingNextPage &&
        hasCachedMetrics &&
        !queryFailed && (
          <p
            className="border-y border-[var(--color-border)] px-2 py-2 text-sm text-[var(--color-text-muted)]"
            role="status"
          >
            Refreshing persisted QC metrics. Cached rows remain visible until the refresh is confirmed.
          </p>
        )}

      <QcMetricList
        metrics={metrics}
        hasNextPage={Boolean(listQuery.hasNextPage) && !cursorRecoveryNeeded}
        isFetchingNextPage={listQuery.isFetchingNextPage}
        onLoadMore={() => void listQuery.fetchNextPage()}
        onOpenSourceArtifact={onOpenSourceArtifact}
      />
    </div>
  );
}

function QcStatus({
  title,
  children,
  tone = 'neutral',
  busy = false,
}: {
  title: string;
  children: ReactNode;
  tone?: 'neutral' | 'error';
  busy?: boolean;
}) {
  return (
    <div
      className={`min-h-32 border-y border-[var(--color-border)] px-2 py-5 ${
        tone === 'error'
          ? 'text-[var(--color-error)]'
          : 'text-[var(--color-text-muted)]'
      }`}
      role="status"
      aria-live="polite"
      aria-busy={busy || undefined}
    >
      <h3 className="text-sm font-semibold text-[var(--color-text)]">{title}</h3>
      <div className="mt-1 max-w-2xl text-sm">{children}</div>
    </div>
  );
}

function StatusRefreshButton({ onClick }: { onClick: () => void }) {
  return (
    <Button
      className="mt-3"
      variant="secondary"
      onClick={onClick}
      aria-label="Refresh QC status"
    >
      <RefreshCw className="mr-1.5 h-4 w-4" aria-hidden="true" />
      Refresh status
    </Button>
  );
}

function QcListSkeleton() {
  return (
    <div
      className="min-h-64 space-y-2"
      role="status"
      aria-label="Loading QC metrics"
    >
      <div className="h-8 animate-pulse rounded bg-gray-100" />
      {Array.from({ length: 5 }, (_, index) => (
        <div
          key={index}
          className="h-12 animate-pulse border-b border-[var(--color-border)] bg-gray-50"
        />
      ))}
    </div>
  );
}
