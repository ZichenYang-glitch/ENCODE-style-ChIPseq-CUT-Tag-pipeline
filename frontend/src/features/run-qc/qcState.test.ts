import { describe, expect, it } from 'vitest';
import type {
  QcMetricResponse,
  RunQcMetricsResponse,
} from '../../api/generated/models';
import type { RunEventResponse } from '../../api/runTypes';
import {
  flattenQcMetricPages,
  formatQcProducedTime,
  isValidQcMetricId,
  qcIndexingOutcome,
  qcPaginationAnomaly,
  safeNextQcCursor,
} from './qcState';

const METRIC_A = `qcmetric-${'a'.repeat(64)}`;
const METRIC_B = `qcmetric-${'b'.repeat(64)}`;

function event(
  sequence: number,
  eventType: string,
  context: Record<string, unknown> = {},
): RunEventResponse {
  return {
    event_id: `event-${sequence}`,
    run_id: 'run-1',
    sequence,
    event_type: eventType,
    timestamp: '2026-07-13T08:00:00Z',
    status: 'succeeded',
    stage: 'qc_summary_indexing',
    message: 'QC outcome changed.',
    context,
    issue: null,
  };
}

function metric(metricId: string): QcMetricResponse {
  return {
    metric_id: metricId,
    metric_key: 'alignment.mapped_reads',
    display_name: 'Mapped reads',
    value: '9007199254740993.125',
    unit: 'count',
    scope: 'sample',
    sample_id: 'sample-1',
    experiment_id: null,
    assay: 'chipseq',
    qc_flag: 'pass',
    source_artifact_id: 'artifact-qc-summary',
    produced_at: '2026-07-13T08:00:00Z',
  };
}

function page(
  metrics: QcMetricResponse[],
  nextCursor: string | null = null,
): RunQcMetricsResponse {
  return {
    ok: true,
    run_id: 'run-1',
    qc_metrics: metrics,
    next_cursor: nextCursor,
    issues: [],
  };
}

describe('qcIndexingOutcome', () => {
  it('uses the latest valid indexed, failed, or invalidated event', () => {
    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexing_failed', {
          reason_code: 'QC_SUMMARY_ADAPTER_FAILED',
        }),
        event(3, 'qc_metrics_indexed', { metric_count: 2 }),
        event(2, 'qc_metrics_invalidated'),
      ]),
    ).toEqual({ kind: 'indexed', count: 2 });

    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexed', { metric_count: 2 }),
        event(2, 'qc_metrics_invalidated'),
      ]),
    ).toEqual({ kind: 'pending' });

    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexed', { metric_count: 2 }),
        event(2, 'qc_metrics_indexing_failed', {
          reason_code: 'QC_SUMMARY_SOURCE_CHANGED',
        }),
      ]),
    ).toEqual({
      kind: 'failed',
      reasonCode: 'QC_SUMMARY_SOURCE_CHANGED',
    });
  });

  it('does not claim an outcome from malformed or incomplete history', () => {
    expect(qcIndexingOutcome([], false)).toEqual({ kind: 'pending' });
    expect(qcIndexingOutcome([], true)).toEqual({ kind: 'unconfirmed' });
    for (const visibleOutcome of [
      event(1, 'qc_metrics_indexed', { metric_count: 0 }),
      event(1, 'qc_metrics_indexing_failed', {
        reason_code: 'QC_SUMMARY_SOURCE_CHANGED',
      }),
      event(1, 'qc_metrics_invalidated'),
    ]) {
      expect(qcIndexingOutcome([visibleOutcome], true)).toEqual({
        kind: 'unconfirmed',
      });
    }
    expect(
      qcIndexingOutcome([event(1, 'qc_metrics_indexed', { metric_count: -1 })]),
    ).toEqual({ kind: 'unconfirmed' });
    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexing_failed', {
          reason_code: '/private/qc.tsv',
        }),
      ]),
    ).toEqual({ kind: 'failed', reasonCode: null });
  });
});

describe('QC keyset helpers', () => {
  it('validates durable metric IDs exactly', () => {
    expect(isValidQcMetricId(METRIC_A)).toBe(true);
    expect(isValidQcMetricId('qcmetric-ABC')).toBe(false);
    expect(isValidQcMetricId('../metric')).toBe(false);
    expect(isValidQcMetricId(null)).toBe(false);
  });

  it('deduplicates metrics without changing first-seen order', () => {
    expect(
      flattenQcMetricPages([
        page([metric(METRIC_A)]),
        page([metric(METRIC_A), metric(METRIC_B)]),
      ]).map((item) => item.metric_id),
    ).toEqual([METRIC_A, METRIC_B]);
  });

  it('continues only with a new valid cursor', () => {
    const first = page([metric(METRIC_A)], METRIC_A);
    expect(safeNextQcCursor(first, [first], [undefined])).toBe(METRIC_A);

    const repeated = page([metric(METRIC_B)], METRIC_A);
    expect(
      safeNextQcCursor(repeated, [first, repeated], [undefined, METRIC_A]),
    ).toBeUndefined();
    expect(
      qcPaginationAnomaly([first, repeated], [undefined, METRIC_A]),
    ).toBe(true);

    const invalid = page([metric(METRIC_A)], '../cursor');
    expect(safeNextQcCursor(invalid, [invalid], [undefined])).toBeUndefined();
    expect(qcPaginationAnomaly([invalid], [undefined])).toBe(true);
    expect(qcPaginationAnomaly([page([metric(METRIC_A)])], [undefined])).toBe(
      false,
    );
  });
});

describe('QC formatting', () => {
  it('uses Intl for valid produced times and a stable invalid fallback', () => {
    expect(formatQcProducedTime('2026-07-13T08:00:00Z', 'en-US')).toContain(
      'Jul',
    );
    expect(formatQcProducedTime('not-a-time', 'en-US')).toBe('Unknown time');
  });
});
