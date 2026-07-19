import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { act, render, screen, waitFor } from '@testing-library/react';
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
const GENERATION_A = `qcgen-${'a'.repeat(64)}`;
const GENERATION_B = `qcgen-${'b'.repeat(64)}`;
const CURSOR_A = `qccur_${'a'.repeat(64)}`;

type GenerationBoundQcMetricsResponse = RunQcMetricsResponse & {
  qc_generation: string;
};
type GenerationBoundIndexedOutcome = Extract<
  ComponentProps<typeof QcWorkbench>['outcome'],
  { kind: 'indexed' }
> & { generation: string };

function indexedOutcome(
  count = 1,
  generation = GENERATION_A,
): GenerationBoundIndexedOutcome {
  return { kind: 'indexed', count, generation };
}

function metric(
  metricId: string,
  name: string,
  value = '9007199254740993.125',
): QcMetricResponse {
  return {
    metric_id: metricId,
    metric_key: `alignment.${name.toLowerCase().replace(/ /g, '_')}`,
    display_name: name,
    value,
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

function renderWorkbench(
  overrides: Partial<ComponentProps<typeof QcWorkbench>> = {},
) {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false }, mutations: { retry: false } },
  });
  const props: ComponentProps<typeof QcWorkbench> = {
    runId: 'run-1',
    runStatus: 'succeeded',
    outcome: indexedOutcome(),
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
    renderWorkbench({ outcome: indexedOutcome(0) });
    expect(await screen.findByText('No indexed QC metrics')).toBeInTheDocument();
  });

  it('does not treat a nonzero indexed count with an empty page as empty', async () => {
    listQcMetricsMock.mockResolvedValue(page([]));
    renderWorkbench({ outcome: indexedOutcome(1) });
    expect(await screen.findByText('QC metric index is not ready')).toBeInTheDocument();
    expect(screen.queryByText('No indexed QC metrics')).not.toBeInTheDocument();
  });
});

describe('QcWorkbench pagination and recovery', () => {
  it('loads explicit pages of 50 and stops a repeated cursor visibly', async () => {
    const user = userEvent.setup();
    listQcMetricsMock
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Mapped reads')], CURSOR_A))
      .mockResolvedValueOnce(page([metric(METRIC_B, 'Duplicate reads')], CURSOR_A));
    renderWorkbench({ outcome: indexedOutcome(2) });

    expect(await screen.findAllByText('Mapped reads')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more QC metrics' }));
    expect(await screen.findAllByText('Duplicate reads')).not.toHaveLength(0);
    expect(listQcMetricsMock).toHaveBeenNthCalledWith(1, 'run-1', {
      after: undefined,
      generation: GENERATION_A,
      limit: 50,
    });
    expect(listQcMetricsMock).toHaveBeenNthCalledWith(2, 'run-1', {
      after: CURSOR_A,
      generation: GENERATION_A,
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
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Mapped reads')], CURSOR_A))
      .mockRejectedValueOnce(cursorError)
      .mockRejectedValueOnce(cursorError)
      .mockResolvedValueOnce(page([metric(METRIC_B, 'Current mapped reads')]));
    renderWorkbench({ outcome: indexedOutcome(1) });

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
      generation: GENERATION_A,
      limit: 50,
    });
  });

  it('keys the metric query by the accepted QC generation', async () => {
    listQcMetricsMock.mockResolvedValue(page([metric(METRIC_A, 'Mapped reads')]));
    const { queryClient } = renderWorkbench({
      outcome: indexedOutcome(1, GENERATION_A),
    });

    expect(await screen.findAllByText('Mapped reads')).not.toHaveLength(0);
    expect(
      queryClient.getQueryCache().find({
        queryKey: ['run-qc-metrics', 'run-1', GENERATION_A],
        exact: true,
      }),
    ).toBeDefined();
    expect(
      queryClient.getQueryCache().find({
        queryKey: ['run-qc-metrics', 'run-1'],
        exact: true,
      }),
    ).toBeUndefined();
  });

  it('clears equal-count cached pages when the generation changes with the same metric ID', async () => {
    listQcMetricsMock
      .mockResolvedValueOnce(
        page(
          [metric(METRIC_A, 'Mapped fragments', '101.125')],
          null,
          GENERATION_A,
        ),
      )
      .mockResolvedValueOnce(
        page(
          [metric(METRIC_A, 'Mapped fragments', '202.25')],
          null,
          GENERATION_B,
        ),
      );
    const rendered = renderWorkbench({
      outcome: indexedOutcome(1, GENERATION_A),
    });
    expect(await screen.findAllByText('101.125')).not.toHaveLength(0);

    const nextProps = {
      ...rendered.props,
      outcome: indexedOutcome(1, GENERATION_B),
    };
    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <QcWorkbench {...nextProps} />
      </QueryClientProvider>,
    );

    await waitFor(() => expect(screen.queryAllByText('101.125')).toHaveLength(0));
    expect(await screen.findAllByText('202.25')).not.toHaveLength(0);
    expect(
      rendered.queryClient
        .getQueryCache()
        .findAll({ queryKey: ['run-qc-metrics', 'run-1'] })
        .map((query) => query.queryKey),
    ).toEqual([['run-qc-metrics', 'run-1', GENERATION_B]]);
  });

  it('clears cached generation pages when the outcome is invalidated', async () => {
    listQcMetricsMock.mockResolvedValue(
      page([metric(METRIC_A, 'Mapped fragments')]),
    );
    const rendered = renderWorkbench({
      outcome: indexedOutcome(1, GENERATION_A),
    });
    expect(await screen.findAllByText('Mapped fragments')).not.toHaveLength(0);

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <QcWorkbench {...rendered.props} outcome={{ kind: 'pending' }} />
      </QueryClientProvider>,
    );

    expect(await screen.findByText('Indexing QC metrics')).toBeInTheDocument();
    await waitFor(() =>
      expect(
        rendered.queryClient.getQueryCache().findAll({
          queryKey: ['run-qc-metrics', 'run-1'],
        }),
      ).toHaveLength(0),
    );
  });

  it('fails closed when a successful response belongs to a different generation', async () => {
    listQcMetricsMock.mockResolvedValue(
      page([metric(METRIC_A, 'Wrong generation metric')], null, GENERATION_B),
    );
    renderWorkbench({ outcome: indexedOutcome(1, GENERATION_A) });

    expect(await screen.findByText('QC metrics could not be loaded')).toBeInTheDocument();
    expect(screen.queryAllByText('Wrong generation metric')).toHaveLength(0);
  });

  it('clears cached pages and refreshes status after a generation-changed response', async () => {
    const user = userEvent.setup();
    const onRefreshStatus = vi.fn();
    const changed = new ApiError(
      409,
      'RUN_QC_METRIC_GENERATION_CHANGED',
      'QC metric generation changed.',
    );
    listQcMetricsMock
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Stale metric')], CURSOR_A))
      .mockRejectedValueOnce(changed)
      .mockRejectedValueOnce(changed);
    renderWorkbench({
      outcome: indexedOutcome(2, GENERATION_A),
      onRefreshStatus,
    });

    expect(await screen.findAllByText('Stale metric')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more QC metrics' }));

    await waitFor(() => expect(onRefreshStatus).toHaveBeenCalledTimes(1));
    expect(screen.queryAllByText('Stale metric')).toHaveLength(0);
  });

  it('retries canonical status after the first 409 refresh fails and resumes only on the new generation', async () => {
    const user = userEvent.setup();
    const changed = new ApiError(
      409,
      'RUN_QC_METRIC_GENERATION_CHANGED',
      'QC metric generation changed.',
    );
    listQcMetricsMock
      .mockResolvedValueOnce(page([metric(METRIC_A, 'Stale metric')], CURSOR_A))
      .mockRejectedValueOnce(changed)
      .mockResolvedValueOnce(
        page([metric(METRIC_A, 'Current metric', '202')], null, GENERATION_B),
      );
    let rendered: ReturnType<typeof renderWorkbench>;
    const onRefreshStatus = vi
      .fn<[], Promise<void>>()
      .mockRejectedValueOnce(new Error('/private/status'))
      .mockImplementationOnce(async () => {
        rendered.rerender(
          <QueryClientProvider client={rendered.queryClient}>
            <QcWorkbench
              {...rendered.props}
              outcome={indexedOutcome(1, GENERATION_B)}
            />
          </QueryClientProvider>,
        );
      });
    rendered = renderWorkbench({
      outcome: indexedOutcome(2, GENERATION_A),
      onRefreshStatus,
    });

    expect(await screen.findAllByText('Stale metric')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more QC metrics' }));
    expect(
      await screen.findByRole('button', { name: 'Refresh QC status' }),
    ).toBeInTheDocument();
    expect(onRefreshStatus).toHaveBeenCalledTimes(1);
    expect(
      listQcMetricsMock.mock.calls.filter(
        ([, parameters]) => parameters?.generation === GENERATION_A,
      ),
    ).toHaveLength(2);

    await user.click(screen.getByRole('button', { name: 'Refresh QC status' }));
    expect(onRefreshStatus).toHaveBeenCalledTimes(2);
    expect(await screen.findAllByText('Current metric')).not.toHaveLength(0);
    expect(screen.queryAllByText('Stale metric')).toHaveLength(0);
    expect(listQcMetricsMock).toHaveBeenLastCalledWith('run-1', {
      after: undefined,
      generation: GENERATION_B,
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
