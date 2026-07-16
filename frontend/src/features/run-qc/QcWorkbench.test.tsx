import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { act, render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { ComponentProps } from 'react';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import type {
  QcMetricResponse,
  RunQcMetricsResponse,
} from '../../api/generated/models';
import { listRunQcMetrics } from '../../api/generated/qc-metrics/qc-metrics';
import { ApiError } from '../../api/fetcher';
import { QcWorkbench } from './QcWorkbench';

vi.mock('../../api/generated/qc-metrics/qc-metrics', () => ({
  listRunQcMetrics: vi.fn(),
}));

const listQcMetricsMock = vi.mocked(listRunQcMetrics);
const METRIC_A = `qcmetric-${'a'.repeat(64)}`;
const METRIC_B = `qcmetric-${'b'.repeat(64)}`;

function metric(metricId: string, name: string): QcMetricResponse {
  return {
    metric_id: metricId,
    metric_key: `alignment.${name.toLowerCase().replace(/ /g, '_')}`,
    display_name: name,
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

function renderWorkbench(
  overrides: Partial<ComponentProps<typeof QcWorkbench>> = {},
) {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false }, mutations: { retry: false } },
  });
  const props: ComponentProps<typeof QcWorkbench> = {
    runId: 'run-1',
    runStatus: 'succeeded',
    outcome: { kind: 'indexed', count: 1 },
    onOpenSourceArtifact: vi.fn(),
    onRefreshStatus: vi.fn(),
    ...overrides,
  };
  return {
    props,
    queryClient,
    ...render(
      <QueryClientProvider client={queryClient}>
        <QcWorkbench {...props} />
      </QueryClientProvider>,
    ),
  };
}

beforeEach(() => {
  listQcMetricsMock.mockReset();
});

afterEach(() => {
  vi.clearAllMocks();
});

describe('QcWorkbench honest states', () => {
  it('does not query an active run or a succeeded run still indexing', () => {
    const active = renderWorkbench({
      runStatus: 'running',
      outcome: { kind: 'pending' },
    });
    expect(screen.getByText('QC metrics are not available yet')).toBeInTheDocument();
    active.unmount();

    renderWorkbench({ outcome: { kind: 'pending' } });
    expect(screen.getByText('Indexing QC metrics')).toBeInTheDocument();
    expect(listQcMetricsMock).not.toHaveBeenCalled();
  });

  it('shows unconfirmed and indexing-failed outcomes with status refresh', async () => {
    const user = userEvent.setup();
    const onRefreshStatus = vi.fn();
    const first = renderWorkbench({
      outcome: { kind: 'unconfirmed' },
      onRefreshStatus,
    });
    await user.click(screen.getByRole('button', { name: 'Refresh QC status' }));
    expect(onRefreshStatus).toHaveBeenCalledTimes(1);
    first.unmount();

    renderWorkbench({
      outcome: {
        kind: 'failed',
        reasonCode: 'QC_SUMMARY_SOURCE_CHANGED',
      },
      onRefreshStatus,
    });
    expect(screen.getByText('QC metric indexing failed')).toBeInTheDocument();
    expect(screen.getByText('Reference: QC_SUMMARY_SOURCE_CHANGED')).toBeInTheDocument();
    expect(screen.queryByText(/private|exception/i)).not.toBeInTheDocument();
  });

  it('shows a stable skeleton while the first page is loading', () => {
    listQcMetricsMock.mockReturnValue(new Promise(() => undefined));
    renderWorkbench();
    expect(screen.getByRole('status', { name: 'Loading QC metrics' })).toBeInTheDocument();
  });

  it('claims empty only after indexed zero and a successful empty page', async () => {
    listQcMetricsMock.mockResolvedValue(page([]));
    renderWorkbench({ outcome: { kind: 'indexed', count: 0 } });
    expect(await screen.findByText('No indexed QC metrics')).toBeInTheDocument();
  });

  it('does not treat a nonzero indexed count with an empty page as empty', async () => {
    listQcMetricsMock.mockResolvedValue(page([]));
    renderWorkbench({ outcome: { kind: 'indexed', count: 1 } });
    expect(await screen.findByText('QC metric index is not ready')).toBeInTheDocument();
    expect(screen.queryByText('No indexed QC metrics')).not.toBeInTheDocument();
  });
});

describe('QcWorkbench pagination and recovery', () => {
  it('loads explicit pages of 50 and stops a repeated cursor visibly', async () => {
    const user = userEvent.setup();
    listQcMetricsMock
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Mapped reads')], METRIC_A))
      .mockResolvedValueOnce(page([metric(METRIC_B, 'Duplicate reads')], METRIC_A));
    renderWorkbench({ outcome: { kind: 'indexed', count: 2 } });

    expect(await screen.findAllByText('Mapped reads')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more QC metrics' }));
    expect(await screen.findAllByText('Duplicate reads')).not.toHaveLength(0);
    expect(listQcMetricsMock).toHaveBeenNthCalledWith(1, 'run-1', {
      after: undefined,
      limit: 50,
    });
    expect(listQcMetricsMock).toHaveBeenNthCalledWith(2, 'run-1', {
      after: METRIC_A,
      limit: 50,
    });
    expect(
      screen.getByText(/pagination stopped because the next cursor was not safe/i),
    ).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Load more QC metrics' })).not.toBeInTheDocument();
  });

  it('reloads from the first page when atomic replacement invalidates a cursor', async () => {
    const user = userEvent.setup();
    const cursorError = new ApiError(
      400,
      'RUN_QC_METRIC_CURSOR_NOT_FOUND',
      'Cursor unavailable.',
    );
    listQcMetricsMock
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Mapped reads')], METRIC_A))
      .mockRejectedValueOnce(cursorError)
      .mockRejectedValueOnce(cursorError)
      .mockResolvedValueOnce(page([metric(METRIC_B, 'Current mapped reads')]));
    renderWorkbench({ outcome: { kind: 'indexed', count: 1 } });

    await screen.findAllByText('Mapped reads');
    await user.click(screen.getByRole('button', { name: 'Load more QC metrics' }));
    expect(
      await screen.findByText(/QC results changed while pages were being loaded/i),
    ).toBeInTheDocument();
    await user.click(
      screen.getByRole('button', { name: 'Reload QC metrics from first page' }),
    );
    expect(await screen.findAllByText('Current mapped reads')).not.toHaveLength(0);
    expect(listQcMetricsMock).toHaveBeenLastCalledWith('run-1', {
      after: undefined,
      limit: 50,
    });
  });

  it('retries an initial redacted transport failure', async () => {
    const user = userEvent.setup();
    listQcMetricsMock
      .mockRejectedValueOnce(new Error('/private/qc.tsv'))
      .mockRejectedValueOnce(new Error('/private/qc.tsv'))
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Mapped reads')]));
    renderWorkbench();

    expect(await screen.findByText('QC metrics could not be loaded')).toBeInTheDocument();
    expect(screen.queryByText('/private/qc.tsv')).not.toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Retry QC metrics' }));
    expect(await screen.findAllByText('Mapped reads')).not.toHaveLength(0);
  });

  it('preserves confirmed rows when a later refresh fails', async () => {
    const queryError = new Error('/private/qc.tsv');
    listQcMetricsMock
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Mapped reads')]))
      .mockRejectedValueOnce(queryError)
      .mockRejectedValueOnce(queryError);
    const { queryClient } = renderWorkbench();

    expect(await screen.findAllByText('Mapped reads')).not.toHaveLength(0);
    await act(async () => {
      await queryClient.invalidateQueries({
        queryKey: ['run-qc-metrics', 'run-1'],
      });
    });
    expect(
      await screen.findByText(/Cached QC rows are preserved and may be stale/i),
    ).toBeInTheDocument();
    expect(screen.getAllByText('Mapped reads').length).toBeGreaterThan(0);
    expect(screen.queryByText('/private/qc.tsv')).not.toBeInTheDocument();
  });

  it('fails closed for a malformed successful envelope', async () => {
    listQcMetricsMock.mockResolvedValue({
      ...page([metric(METRIC_A, 'Mapped reads')]),
      run_id: 'run-other',
    });
    renderWorkbench();
    expect(await screen.findByText('QC metrics could not be loaded')).toBeInTheDocument();
    expect(screen.queryByText('Mapped reads')).not.toBeInTheDocument();
  });

  it('accepts the workflow-neutral score unit and rejects unknown units', async () => {
    const scoreMetric = {
      ...metric(METRIC_A, 'TIN mean score'),
      value: '72.125',
      unit: 'score' as const,
    };
    listQcMetricsMock.mockResolvedValueOnce(page([scoreMetric]));
    const accepted = renderWorkbench();
    expect(await screen.findAllByText('TIN mean score')).not.toHaveLength(0);
    expect(screen.getAllByText('score').length).toBeGreaterThan(0);
    accepted.unmount();

    listQcMetricsMock.mockResolvedValueOnce(
      page([
        {
          ...metric(METRIC_B, 'Unknown unit'),
          unit: 'percent',
        } as unknown as QcMetricResponse,
      ]),
    );
    renderWorkbench();
    expect(await screen.findByText('QC metrics could not be loaded')).toBeInTheDocument();
    expect(screen.queryByText('Unknown unit')).not.toBeInTheDocument();
  });
});
