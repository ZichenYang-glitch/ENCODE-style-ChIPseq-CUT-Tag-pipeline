import type {
  QcMetricResponse,
  RunQcMetricsResponse,
} from '../../api/generated/models';
import type { RunEventResponse } from '../../api/runTypes';

const QC_METRIC_ID_PATTERN = /^qcmetric-[0-9a-f]{64}$/;
const QC_GENERATION_PATTERN = /^qcgen-[0-9a-f]{64}$/;
const QC_CURSOR_PATTERN = /^qccur_[A-Za-z0-9_-]{1,1018}$/;
const REASON_CODE_PATTERN = /^[A-Z][A-Z0-9_]{0,127}$/;

export type QcIndexingOutcome =
  | { kind: 'pending' }
  | { kind: 'unconfirmed' }
  | { kind: 'indexed'; count: number; generation: string }
  | { kind: 'failed'; reasonCode: string | null };

export function isValidQcMetricId(
  value: string | null | undefined,
): value is string {
  return typeof value === 'string' && QC_METRIC_ID_PATTERN.test(value);
}

export function isValidQcGeneration(
  value: unknown,
): value is string {
  return typeof value === 'string' && QC_GENERATION_PATTERN.test(value);
}

export function isValidQcCursor(
  value: string | null | undefined,
): value is string {
  return typeof value === 'string' && QC_CURSOR_PATTERN.test(value);
}

export function qcIndexingOutcome(
  events: RunEventResponse[],
  truncated = false,
): QcIndexingOutcome {
  if (truncated) return { kind: 'unconfirmed' };
  const latest = events
    .filter((item) =>
      [
        'qc_metrics_indexed',
        'qc_metrics_indexing_failed',
        'qc_metrics_invalidated',
      ].includes(item.event_type),
    )
    .sort((left, right) => right.sequence - left.sequence)[0];

  if (!latest) return { kind: 'pending' };
  if (latest.event_type === 'qc_metrics_invalidated') return { kind: 'pending' };
  if (latest.event_type === 'qc_metrics_indexing_failed') {
    const reasonCode = latest.context.reason_code;
    return {
      kind: 'failed',
      reasonCode:
        typeof reasonCode === 'string' && REASON_CODE_PATTERN.test(reasonCode)
          ? reasonCode
          : null,
    };
  }

  const count = latest.context.metric_count;
  const generation = latest.context.qc_generation;
  if (typeof count !== 'number' || !Number.isSafeInteger(count) || count < 0) {
    return { kind: 'unconfirmed' };
  }
  if (!isValidQcGeneration(generation)) {
    return { kind: 'unconfirmed' };
  }
  return { kind: 'indexed', count, generation };
}

export function flattenQcMetricPages(
  pages: RunQcMetricsResponse[],
): QcMetricResponse[] {
  const generations = new Set(pages.map((page) => page.qc_generation));
  if (
    generations.size > 1 ||
    [...generations].some((generation) => !isValidQcGeneration(generation))
  ) {
    throw new Error('QC pagination generation mismatch');
  }
  const seen = new Set<string>();
  const result: QcMetricResponse[] = [];
  for (const page of pages) {
    for (const metric of page.qc_metrics ?? []) {
      if (seen.has(metric.metric_id)) {
        throw new Error('duplicate QC metric in paginated result');
      }
      seen.add(metric.metric_id);
      result.push(metric);
    }
  }
  return result;
}

function pageCursorIsAnomalous(
  page: RunQcMetricsResponse,
  pageIndex: number,
  pageParams: unknown[],
): boolean {
  const next = page.next_cursor;
  if (next == null) return false;
  if (!isValidQcCursor(next)) return true;
  const currentPageParam = pageParams[pageIndex];
  if (currentPageParam === next) return true;
  return pageParams.slice(0, pageIndex).some((value) => value === next);
}

export function qcPaginationAnomaly(
  pages: RunQcMetricsResponse[],
  pageParams: unknown[],
): boolean {
  return pages.some((page, index) =>
    pageCursorIsAnomalous(page, index, pageParams),
  );
}

export function safeNextQcCursor(
  lastPage: RunQcMetricsResponse,
  allPages: RunQcMetricsResponse[],
  pageParams: unknown[],
): string | undefined {
  const pageIndex = allPages.length - 1;
  if (pageCursorIsAnomalous(lastPage, pageIndex, pageParams)) return undefined;
  return lastPage.next_cursor ?? undefined;
}

export function formatQcProducedTime(value: string, locale?: string): string {
  const timestamp = new Date(value);
  if (Number.isNaN(timestamp.getTime())) return 'Unknown time';
  return new Intl.DateTimeFormat(locale, {
    dateStyle: 'medium',
    timeStyle: 'short',
  }).format(timestamp);
}
