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
const GENERATION_A = `qcgen-${'a'.repeat(64)}`;
const GENERATION_B = `qcgen-${'b'.repeat(64)}`;
const CURSOR_A = `qccur_${'a'.repeat(64)}`;

type GenerationBoundQcMetricsResponse = RunQcMetricsResponse & {
  qc_generation: string;
};

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
  generation = GENERATION_A,
): GenerationBoundQcMetricsResponse {
  return {
    ok: true,
    run_id: 'run-1',
    qc_generation: generation,
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
        event(3, 'qc_metrics_indexed', {
          metric_count: 2,
          qc_generation: GENERATION_A,
        }),
        event(2, 'qc_metrics_invalidated'),
      ]),
    ).toEqual({ kind: 'indexed', count: 2, generation: GENERATION_A });

    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexed', {
          metric_count: 2,
          qc_generation: GENERATION_A,
        }),
        event(2, 'qc_metrics_invalidated'),
      ]),
    ).toEqual({ kind: 'pending' });

    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexed', {
          metric_count: 2,
          qc_generation: GENERATION_A,
        }),
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
      event(1, 'qc_metrics_indexed', {
        metric_count: 0,
        qc_generation: GENERATION_A,
      }),
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
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexed', {
          metric_count: -1,
          qc_generation: GENERATION_A,
        }),
      ]),
    ).toEqual({ kind: 'unconfirmed' });
    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexing_failed', {
          reason_code: '/private/qc.tsv',
        }),
      ]),
    ).toEqual({ kind: 'failed', reasonCode: null });
  });

  it('requires a public-safe generation on every indexed outcome', () => {
    for (const qcGeneration of [
      undefined,
      '',
      '/private/qc-generation',
      `qcgen-${'A'.repeat(64)}`,
      `qcgen-${'a'.repeat(63)}`,
    ]) {
      expect(
        qcIndexingOutcome([
          event(1, 'qc_metrics_indexed', {
            metric_count: 1,
            ...(qcGeneration === undefined
              ? {}
              : { qc_generation: qcGeneration }),
          }),
        ]),
      ).toEqual({ kind: 'unconfirmed' });
    }
  });

  it('distinguishes equal-count indexed outcomes by generation', () => {
    expect(
      qcIndexingOutcome([
        event(1, 'qc_metrics_indexed', {
          metric_count: 1,
          qc_generation: GENERATION_A,
        }),
        event(2, 'qc_metrics_indexed', {
          metric_count: 1,
          qc_generation: GENERATION_B,
        }),
      ]),
    ).toEqual({ kind: 'indexed', count: 1, generation: GENERATION_B });
  });
});

describe('QC keyset helpers', () => {
  it('validates durable metric IDs exactly', () => {
    expect(isValidQcMetricId(METRIC_A)).toBe(true);
    expect(isValidQcMetricId('qcmetric-ABC')).toBe(false);
    expect(isValidQcMetricId('../metric')).toBe(false);
    expect(isValidQcMetricId(null)).toBe(false);
  });

  it('rejects duplicate metric identities across pages instead of hiding them', () => {
    expect(() =>
      flattenQcMetricPages([
        page([metric(METRIC_A)]),
        page([metric(METRIC_A), metric(METRIC_B)]),
      ]),
    ).toThrow(/duplicate|pagination|generation/i);
  });

  it('rejects pages from different QC generations', () => {
    expect(() =>
      flattenQcMetricPages([
        page([metric(METRIC_A)], null, GENERATION_A),
        page([metric(METRIC_B)], null, GENERATION_B),
      ]),
    ).toThrow(/generation/i);
  });

  it('continues only with a new valid cursor', () => {
    const first = page([metric(METRIC_A)], CURSOR_A);
    expect(safeNextQcCursor(first, [first], [undefined])).toBe(CURSOR_A);

    const repeated = page([metric(METRIC_B)], CURSOR_A);
    expect(
      safeNextQcCursor(repeated, [first, repeated], [undefined, CURSOR_A]),
    ).toBeUndefined();
    expect(
      qcPaginationAnomaly([first, repeated], [undefined, CURSOR_A]),
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
