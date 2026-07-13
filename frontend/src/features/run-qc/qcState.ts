import type {
  QcMetricResponse,
  RunQcMetricsResponse,
} from '../../api/generated/models';
import type { RunEventResponse } from '../../api/runTypes';

const QC_METRIC_ID_PATTERN = /^qcmetric-[0-9a-f]{64}$/;
const REASON_CODE_PATTERN = /^[A-Z][A-Z0-9_]{0,127}$/;

export type QcIndexingOutcome =
  | { kind: 'pending' }
  | { kind: 'unconfirmed' }
  | { kind: 'indexed'; count: number }
  | { kind: 'failed'; reasonCode: string | null };

export function isValidQcMetricId(
  value: string | null | undefined,
): value is string {
  return typeof value === 'string' && QC_METRIC_ID_PATTERN.test(value);
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
  if (typeof count !== 'number' || !Number.isSafeInteger(count) || count < 0) {
    return { kind: 'unconfirmed' };
  }
  return { kind: 'indexed', count };
}

export function flattenQcMetricPages(
  pages: RunQcMetricsResponse[],
): QcMetricResponse[] {
  const seen = new Set<string>();
  const result: QcMetricResponse[] = [];
  for (const page of pages) {
    for (const metric of page.qc_metrics ?? []) {
      if (seen.has(metric.metric_id)) continue;
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
  if (!isValidQcMetricId(next)) return true;
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
