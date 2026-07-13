import { describe, it, expect, vi, type Mock } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import userEvent from '@testing-library/user-event';
import {
  isRunPollingPaused,
  RunProgressPanel,
  shouldPollRunSnapshot,
  type RunSnapshot,
} from './RunProgressPanel';
import type { RunApiClient } from '../../api/runClient';
import type { ValidateWorkflowResponse, WorkflowInputs } from '../../api/types';
import type { RunEventResponse, RunResponse } from '../../api/runTypes';

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

const validatedInputs: WorkflowInputs = {
  config: { genome: 'hg38' },
  samples: [{ name: 'sample-1', fastq_r1: 's1_R1.fq.gz' }],
  options: { threads: 4 },
};

const successfulValidation: ValidateWorkflowResponse = {
  ok: true,
  workflow_id: WORKFLOW_ID,
  value: null,
  issues: [],
};

const failedValidation: ValidateWorkflowResponse = {
  ok: false,
  workflow_id: WORKFLOW_ID,
  value: null,
  issues: [
    {
      code: 'SAMPLE_INVALID',
      message: 'Sample sheet is invalid.',
      severity: 'error',
      path: 'samples',
      source: 'stub',
      technical_message: null,
      hint: null,
      context: {},
    },
  ],
};

function createMockRunClient(): RunApiClient {
  const runs = new Map<string, { status: string; run_id: string; cancelled: boolean }>();

  return {
    createRun: vi.fn().mockImplementation(async (workflowId, request) => {
      const runId = `run-${runs.size + 1}`;
      runs.set(runId, { status: 'created', run_id: runId, cancelled: false });
      return {
        ok: true,
        run: {
          run_id: runId,
          workflow_id: workflowId,
          inputs: {
            config: request.config,
            samples: request.samples,
            options: request.options ?? {},
          },
          status: 'created',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:00:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: null,
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      };
    }),
    preflightRun: vi.fn().mockImplementation(async (runId) => {
      const record = runs.get(runId);
      if (!record) {
        return {
          ok: false,
          run: null,
          issues: [{ code: 'RUN_NOT_FOUND', message: 'Not found.' }],
        };
      }
      record.status = 'planned';
      return {
        ok: true,
        run: {
          run_id: runId,
          workflow_id: WORKFLOW_ID,
          inputs: validatedInputs,
          status: 'planned',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:01:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: 'preflight',
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      };
    }),
    startRun: vi.fn().mockImplementation(async (runId) => {
      const record = runs.get(runId);
      if (!record) return { ok: false, run: null, issues: [] };
      record.status = 'queued';
      return {
        ok: true,
        run: {
          run_id: runId,
          workflow_id: WORKFLOW_ID,
          inputs: validatedInputs,
          status: 'queued',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:02:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: 'execution',
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      };
    }),
    getRun: vi.fn().mockImplementation(async (runId) => {
      const record = runs.get(runId);
      if (!record) {
        return { ok: false, run: null, issues: [{ code: 'RUN_NOT_FOUND', message: 'Not found.' }] };
      }
      return {
        ok: true,
        run: {
          run_id: runId,
          workflow_id: WORKFLOW_ID,
          inputs: validatedInputs,
          status: record.status,
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: record.cancelled
            ? '2026-07-04T12:01:00.000Z'
            : '2026-07-04T12:00:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: null,
          cancellation_reason: record.cancelled ? 'User requested cancellation.' : null,
          error: null,
          tags: {},
        },
        issues: [],
      };
    }),
    listRunEvents: vi.fn().mockImplementation(async (runId) => {
      const record = runs.get(runId);
      if (!record) {
        return { ok: false, run_id: runId, events: [], next_cursor: null, issues: [] };
      }
      const events: RunEventResponse[] = [
        {
          event_id: `${runId}-evt-1`,
          run_id: runId,
          sequence: 1,
          event_type: 'status_changed',
          timestamp: '2026-07-04T12:00:00.000Z',
          status: 'created',
          stage: null,
          message: 'Run created.',
          context: { previous_status: null, new_status: 'created' },
          issue: null,
        },
      ];
      if (record.cancelled) {
        events.push({
          event_id: `${runId}-evt-2`,
          run_id: runId,
          sequence: 2,
          event_type: 'status_changed',
          timestamp: '2026-07-04T12:01:00.000Z',
          status: 'cancelled',
          stage: null,
          message: 'Run cancelled.',
          context: {
            previous_status: 'created',
            new_status: 'cancelled',
            cancellation_reason: 'User requested cancellation.',
          },
          issue: null,
        });
      }
      return { ok: true, run_id: runId, events, next_cursor: null, issues: [] };
    }),
    listRunLogs: vi.fn().mockImplementation(async (runId, options = {}) => {
      const streamName = options.streamName ?? 'stdout';
      const chunks = streamName === 'stderr'
        ? []
        : [
            {
              chunk_id: `${runId}-log-1`,
              run_id: runId,
              stream_name: 'stdout',
              sequence: 1,
              timestamp: '2026-07-04T12:00:00.000Z',
              lines: ['[stub] validating inputs'],
            },
          ];
      return {
        ok: true,
        run_id: runId,
        stream_name: streamName,
        chunks,
        next_cursor: null,
        issues: [],
      };
    }),
    cancelRun: vi.fn().mockImplementation(async (runId) => {
      const record = runs.get(runId);
      if (!record) {
        return { ok: false, run: null, issues: [{ code: 'RUN_NOT_FOUND', message: 'Not found.' }] };
      }
      record.status = 'cancelled';
      record.cancelled = true;
      return {
        ok: true,
        run: {
          run_id: runId,
          workflow_id: WORKFLOW_ID,
          inputs: validatedInputs,
          status: 'cancelled',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:01:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: null,
          cancellation_reason: 'User requested cancellation.',
          error: null,
          tags: {},
        },
        issues: [],
      };
    }),
  };
}

function renderPanel(props: Partial<Parameters<typeof RunProgressPanel>[0]> = {}) {
  const runClient = createMockRunClient();
  const defaultProps = {
    workflowId: WORKFLOW_ID,
    validationResult: null,
    validatedInputs: null,
    runClient,
  };
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });

  function Wrapper({ children }: { children: React.ReactNode }) {
    return <QueryClientProvider client={queryClient}>{children}</QueryClientProvider>;
  }

  const result = render(
    <Wrapper>
      <RunProgressPanel {...defaultProps} {...props} />
    </Wrapper>,
  );

  return {
    runClient,
    ...result,
    rerender(nextProps: Partial<Parameters<typeof RunProgressPanel>[0]> = {}) {
      result.rerender(
        <Wrapper>
          <RunProgressPanel {...defaultProps} {...nextProps} />
        </Wrapper>,
      );
    },
  };
}

describe('RunProgressPanel', () => {
  it('marks active polling stale only after the bounded observation window', () => {
    const startedAt = 1_000;

    expect(isRunPollingPaused('running', startedAt, startedAt + 899_999)).toBe(false);
    expect(isRunPollingPaused('running', startedAt, startedAt + 900_000)).toBe(true);
    expect(isRunPollingPaused('succeeded', startedAt, startedAt + 900_000)).toBe(false);
  });

  it('uses one polling eligibility for lifecycle and succeeded indexing states', () => {
    const succeeded: RunSnapshot = {
      run: {
        run_id: 'run-1',
        workflow_id: WORKFLOW_ID,
        inputs: {},
        status: 'succeeded',
        created_at: '2026-07-12T12:00:00Z',
        updated_at: '2026-07-12T12:01:00Z',
        started_at: '2026-07-12T12:00:10Z',
        ended_at: '2026-07-12T12:01:00Z',
        current_stage: 'execution',
        cancellation_reason: null,
        error: null,
        tags: {},
      },
      events: [],
      eventsTruncated: false,
      stdoutLogs: [],
      stderrLogs: [],
      issues: [],
      truncated: false,
    };

    expect(shouldPollRunSnapshot(succeeded, 'activity')).toBe(false);
    expect(shouldPollRunSnapshot(succeeded, 'artifacts')).toBe(true);
    expect(shouldPollRunSnapshot(succeeded, 'qc')).toBe(true);
    expect(
      shouldPollRunSnapshot(
        {
          ...succeeded,
          events: [
            {
              event_id: 'event-indexed',
              run_id: 'run-1',
              sequence: 1,
              event_type: 'artifacts_indexed',
              timestamp: '2026-07-12T12:01:00Z',
              status: 'succeeded',
              stage: 'artifact_extraction',
              message: 'Artifacts indexed.',
              context: { artifact_count: 0 },
              issue: null,
            },
          ],
        },
        'artifacts',
      ),
    ).toBe(false);
    expect(
      shouldPollRunSnapshot(
        { ...succeeded, truncated: true, eventsTruncated: false },
        'artifacts',
      ),
    ).toBe(true);
    expect(
      shouldPollRunSnapshot(
        { ...succeeded, truncated: true, eventsTruncated: true },
        'artifacts',
      ),
    ).toBe(false);
    expect(
      shouldPollRunSnapshot(
        {
          ...succeeded,
          events: [
            {
              event_id: 'event-qc-indexed',
              run_id: 'run-1',
              sequence: 2,
              event_type: 'qc_metrics_indexed',
              timestamp: '2026-07-12T12:01:01Z',
              status: 'succeeded',
              stage: 'qc_summary_indexing',
              message: 'QC metrics indexed.',
              context: { metric_count: 0 },
              issue: null,
            },
          ],
        },
        'qc',
      ),
    ).toBe(false);
    expect(
      shouldPollRunSnapshot(
        {
          ...succeeded,
          events: [
            {
              event_id: 'event-qc-invalidated',
              run_id: 'run-1',
              sequence: 3,
              event_type: 'qc_metrics_invalidated',
              timestamp: '2026-07-12T12:01:02Z',
              status: 'succeeded',
              stage: 'qc_summary_indexing',
              message: 'QC metrics invalidated.',
              context: {},
              issue: null,
            },
          ],
        },
        'qc',
      ),
    ).toBe(true);
    expect(
      shouldPollRunSnapshot(
        {
          ...succeeded,
          events: [
            {
              event_id: 'event-qc-failed',
              run_id: 'run-1',
              sequence: 4,
              event_type: 'qc_metrics_indexing_failed',
              timestamp: '2026-07-12T12:01:03Z',
              status: 'succeeded',
              stage: 'qc_summary_indexing',
              message: 'QC indexing failed.',
              context: { reason_code: 'QC_SUMMARY_SOURCE_CHANGED' },
              issue: null,
            },
          ],
        },
        'qc',
      ),
    ).toBe(false);
    expect(
      shouldPollRunSnapshot(
        { ...succeeded, eventsTruncated: true },
        'qc',
      ),
    ).toBe(false);
  });

  it('shows the no-run state and disables Create run before validation', () => {
    renderPanel();

    expect(
      screen.getByText(/Validate inputs before creating a run record/i),
    ).toBeInTheDocument();
    expect(screen.getByTestId('create-run-button')).toBeDisabled();
  });

  it('keeps Create run disabled after failed validation', () => {
    renderPanel({ validationResult: failedValidation, validatedInputs });

    expect(screen.getByTestId('create-run-button')).toBeDisabled();
  });

  it('enables Create run after successful validation', () => {
    renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    expect(screen.getByTestId('create-run-button')).toBeEnabled();
  });

  it('creates a run record and shows status, events, and stdout logs', async () => {
    const user = userEvent.setup();
    const { runClient } = renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));

    await waitFor(() => {
      expect(runClient.createRun).toHaveBeenCalledTimes(1);
    });
    expect(runClient.createRun).toHaveBeenCalledWith(WORKFLOW_ID, validatedInputs);

    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent('planned');
    expect(runClient.preflightRun).toHaveBeenCalledWith('run-1');
    expect(screen.getByTestId('run-event-feed')).toBeInTheDocument();
    expect(screen.getByText(/status_changed/i)).toBeInTheDocument();
    expect(screen.getByText(/\[stub\] validating inputs/i)).toBeInTheDocument();
  });

  it('reports a malformed successful create response', async () => {
    const runClient = createMockRunClient();
    vi.mocked(runClient.createRun).mockResolvedValue({
      ok: true,
      run: null,
      issues: [],
    });
    const user = userEvent.setup();

    renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
      runClient,
    });
    await user.click(screen.getByTestId('create-run-button'));

    expect(await screen.findByRole('alert')).toHaveTextContent(
      'The run record was missing.',
    );
    expect(runClient.preflightRun).not.toHaveBeenCalled();
  });

  it('rejects validated samples that are not compatible with run creation', async () => {
    const user = userEvent.setup();
    const invalidValidatedInputs: WorkflowInputs = {
      ...validatedInputs,
      samples: [
        { name: 'sample-1', replicate: 1 },
      ] as unknown as WorkflowInputs['samples'],
    };
    const { runClient } = renderPanel({
      validationResult: successfulValidation,
      validatedInputs: invalidValidatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));

    await waitFor(() => {
      expect(runClient.createRun).not.toHaveBeenCalled();
    });
    expect(screen.getByRole('alert')).toHaveTextContent(
      /run record could not be created/i,
    );
  });

  it('cancels a run and updates the status and event feed', async () => {
    const user = userEvent.setup();
    const { runClient } = renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));
    await waitFor(() => {
      expect(runClient.createRun).toHaveBeenCalledTimes(1);
    });

    await user.click(await screen.findByTestId('cancel-run-button'));

    await waitFor(() => {
      expect(runClient.cancelRun).toHaveBeenCalledTimes(1);
    });
    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent('cancelled');
    expect(screen.getAllByText(/status_changed/i)).toHaveLength(2);
  });

  it('refreshes run, events, and logs on Refresh click', async () => {
    const user = userEvent.setup();
    const { runClient } = renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));
    await waitFor(() => {
      expect(runClient.createRun).toHaveBeenCalledTimes(1);
    });

    await user.click(screen.getByTestId('refresh-run-button'));

    await waitFor(() => {
      expect(vi.mocked(runClient.getRun).mock.calls.length).toBeGreaterThanOrEqual(2);
      expect(vi.mocked(runClient.listRunEvents).mock.calls.length).toBeGreaterThanOrEqual(2);
      expect(vi.mocked(runClient.listRunLogs).mock.calls.length).toBeGreaterThanOrEqual(4);
    });
  });

  it('fetches stdout and stderr logs when creating a run', async () => {
    const user = userEvent.setup();
    const { runClient } = renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));

    await waitFor(() => {
      expect(vi.mocked(runClient.listRunLogs).mock.calls.length).toBeGreaterThanOrEqual(2);
    });

    const stdoutCalls = (runClient.listRunLogs as Mock).mock.calls.filter(
      ([, options]) => options?.streamName === 'stdout',
    );
    const stderrCalls = (runClient.listRunLogs as Mock).mock.calls.filter(
      ([, options]) => options?.streamName === 'stderr',
    );

    expect(stdoutCalls.length).toBeGreaterThanOrEqual(1);
    expect(stderrCalls.length).toBeGreaterThanOrEqual(1);
  });

  it('resets log stream selector to stdout when workflow or inputs change', async () => {
    const user = userEvent.setup();
    const newRunClient = createMockRunClient();
    const { rerender } = renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));
    await waitFor(() => {
      expect(screen.getByTestId('stdout-tab')).toHaveAttribute('aria-selected', 'true');
    });

    await user.click(screen.getByTestId('stderr-tab'));
    expect(screen.getByTestId('stderr-tab')).toHaveAttribute('aria-selected', 'true');

    rerender({
      workflowId: 'other-workflow',
      validationResult: successfulValidation,
      validatedInputs: { ...validatedInputs, samples: [{ name: 'sample-2', fastq_r1: 's2_R1.fq.gz' }] },
      runClient: newRunClient,
    });

    // Inputs changed, so the existing run is cleared. Creating a new run should
    // show the selector reset to stdout.
    await user.click(screen.getByTestId('create-run-button'));
    await waitFor(() => {
      expect(screen.getByTestId('stdout-tab')).toHaveAttribute('aria-selected', 'true');
    });
  });

  it('resets log stream selector to stdout when a new run is created with the same workflow and inputs', async () => {
    const user = userEvent.setup();
    renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    await user.click(screen.getByTestId('create-run-button'));
    await waitFor(() => {
      expect(screen.getByTestId('stdout-tab')).toHaveAttribute('aria-selected', 'true');
    });

    await user.click(screen.getByTestId('stderr-tab'));
    expect(screen.getByTestId('stderr-tab')).toHaveAttribute('aria-selected', 'true');

    await user.click(screen.getByTestId('create-run-button'));
    await waitFor(() => {
      expect(screen.getByTestId('stdout-tab')).toHaveAttribute('aria-selected', 'true');
    });
  });

  it('loads an existing run when runId is provided', async () => {
    const runClient = createMockRunClient();
    await runClient.createRun(WORKFLOW_ID, {
      config: validatedInputs.config,
      samples: validatedInputs.samples as Record<string, string>[],
      options: validatedInputs.options,
    });

    renderPanel({
      runId: 'run-1',
      runClient,
    });

    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent(
      'created',
    );
    expect(screen.getByTestId('run-event-feed')).toBeInTheDocument();
    expect(screen.getByText(/\[stub\] validating inputs/i)).toBeInTheDocument();
  });

  it('displays a safe issue message when getRun returns ok=false', async () => {
    const runClient: RunApiClient = {
      createRun: vi.fn(),
      preflightRun: vi.fn(),
      startRun: vi.fn(),
      getRun: vi.fn().mockResolvedValue({
        ok: false,
        run: null,
        issues: [
          {
            code: 'RUN_NOT_FOUND',
            message: 'Run was not found.',
            severity: 'error',
            path: 'run_id',
            source: 'stub',
            technical_message: null,
            hint: null,
            context: { run_id: 'missing-run' },
          },
        ],
      } as RunResponse),
      listRunEvents: vi.fn(),
      listRunLogs: vi.fn(),
      cancelRun: vi.fn(),
    };

    renderPanel({
      runId: 'missing-run',
      runClient,
    });

    expect(await screen.findByTestId('run-progress-error')).toHaveTextContent(
      /RUN_NOT_FOUND: Run was not found/i,
    );
    expect(screen.queryByText(/Validate inputs before creating a run record/i)).not.toBeInTheDocument();
  });

  it('does not treat failed events/logs as successful empty arrays', async () => {
    const runClient: RunApiClient = {
      createRun: vi.fn(),
      preflightRun: vi.fn(),
      startRun: vi.fn(),
      getRun: vi.fn().mockResolvedValue({
        ok: true,
        run: {
          run_id: 'run-1',
          workflow_id: WORKFLOW_ID,
          inputs: validatedInputs as unknown as Record<string, unknown>,
          status: 'created',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:00:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: null,
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      } as unknown as RunResponse),
      listRunEvents: vi.fn().mockResolvedValue({
        ok: false,
        run_id: 'run-1',
        events: [],
        next_cursor: null,
        issues: [
          {
            code: 'RUN_EVENTS_UNAVAILABLE',
            message: 'Events unavailable.',
            severity: 'error',
            path: 'events',
            source: 'stub',
            technical_message: null,
            hint: null,
            context: {},
          },
        ],
      }),
      listRunLogs: vi.fn().mockResolvedValue({
        ok: true,
        run_id: 'run-1',
        stream_name: 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      }),
      cancelRun: vi.fn(),
    };

    renderPanel({
      runId: 'run-1',
      runClient,
    });

    expect(await screen.findByTestId('run-progress-error')).toHaveTextContent(
      /RUN_EVENTS_UNAVAILABLE: Events unavailable/i,
    );
    expect(screen.getByTestId('run-event-feed-empty')).toBeInTheDocument();
  });

  it('does not treat failed logs as successful empty arrays', async () => {
    const runClient: RunApiClient = {
      createRun: vi.fn(),
      preflightRun: vi.fn(),
      startRun: vi.fn(),
      getRun: vi.fn().mockResolvedValue({
        ok: true,
        run: {
          run_id: 'run-1',
          workflow_id: WORKFLOW_ID,
          inputs: validatedInputs as unknown as Record<string, unknown>,
          status: 'created',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:00:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: null,
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      } as unknown as RunResponse),
      listRunEvents: vi.fn().mockResolvedValue({
        ok: true,
        run_id: 'run-1',
        events: [],
        next_cursor: null,
        issues: [],
      }),
      listRunLogs: vi.fn().mockImplementation(async (_, options = {}) => {
        const streamName = options.streamName ?? 'stdout';
        if (streamName === 'stderr') {
          return {
            ok: false,
            run_id: 'run-1',
            stream_name: 'stderr',
            chunks: [],
            next_cursor: null,
            issues: [
              {
                code: 'RUN_LOGS_UNAVAILABLE',
                message: 'Stderr logs unavailable.',
                severity: 'error',
                path: 'logs',
                source: 'stub',
                technical_message: null,
                hint: null,
                context: {},
              },
            ],
          };
        }
        return {
          ok: true,
          run_id: 'run-1',
          stream_name: 'stdout',
          chunks: [],
          next_cursor: null,
          issues: [],
        };
      }),
      cancelRun: vi.fn(),
    };

    renderPanel({
      runId: 'run-1',
      runClient,
    });

    expect(await screen.findByTestId('run-progress-error')).toHaveTextContent(
      /RUN_LOGS_UNAVAILABLE: Stderr logs unavailable/i,
    );
  });

  it('ignores stale responses after runId changes', async () => {
    const runClient: RunApiClient = {
      createRun: vi.fn(),
      preflightRun: vi.fn(),
      startRun: vi.fn(),
      getRun: vi.fn().mockImplementation(async (runId) => {
        const delay = runId === 'slow-run' ? 300 : 10;
        await new Promise((resolve) => setTimeout(resolve, delay));
        return {
          ok: true,
          run: {
            run_id: runId,
            workflow_id: WORKFLOW_ID,
            inputs: validatedInputs as unknown as Record<string, unknown>,
            status: 'created',
            created_at: '2026-07-04T12:00:00.000Z',
            updated_at: '2026-07-04T12:00:00.000Z',
            started_at: null,
            ended_at: null,
            current_stage: null,
            cancellation_reason: null,
            error: null,
            tags: {},
          },
          issues: [],
        } as unknown as RunResponse;
      }),
      listRunEvents: vi.fn().mockImplementation(async (runId) => {
        const delay = runId === 'slow-run' ? 300 : 10;
        await new Promise((resolve) => setTimeout(resolve, delay));
        return {
          ok: true,
          run_id: runId,
          events: [
            {
              event_id: `${runId}-evt-1`,
              run_id: runId,
              sequence: 1,
              event_type: 'status_changed',
              timestamp: '2026-07-04T12:00:00.000Z',
              status: 'created',
              stage: null,
              message: `Run ${runId} created.`,
              context: { previous_status: null, new_status: 'created' },
              issue: null,
            },
          ],
          next_cursor: null,
          issues: [],
        };
      }),
      listRunLogs: vi.fn().mockResolvedValue({
        ok: true,
        run_id: 'run-1',
        stream_name: 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      }),
      cancelRun: vi.fn(),
    };

    const { rerender } = renderPanel({
      runId: 'slow-run',
      runClient,
    });

    // Immediately switch to a faster run before slow-run responses arrive.
    rerender({ runId: 'fast-run', runClient });

    await waitFor(() => {
      expect(screen.getByTestId('run-status-badge')).toHaveTextContent('created');
    });

    await waitFor(() => {
      expect(screen.getByText(/Run fast-run created/i)).toBeInTheDocument();
    });
    expect(screen.queryByText(/Run slow-run created/i)).not.toBeInTheDocument();

    // Wait long enough that the stale slow-run responses would have landed.
    await new Promise((resolve) => setTimeout(resolve, 400));
    expect(screen.queryByText(/Run slow-run created/i)).not.toBeInTheDocument();
  });

  it('calls getRun only once per runId on deep link', async () => {
    const runClient = createMockRunClient();
    await runClient.createRun(WORKFLOW_ID, {
      config: validatedInputs.config,
      samples: validatedInputs.samples as Record<string, string>[],
      options: validatedInputs.options,
    });

    renderPanel({
      runId: 'run-1',
      runClient,
    });

    await waitFor(() => {
      expect(screen.getByTestId('run-status-badge')).toBeInTheDocument();
    });

    expect(runClient.getRun).toHaveBeenCalledTimes(1);
  });

  it('polls an active preflight until it reaches planned', async () => {
    const baseClient = createMockRunClient();
    let getRunCalls = 0;
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => {
        getRunCalls += 1;
        return {
          ok: true,
          run: {
            run_id: 'run-polling',
            workflow_id: WORKFLOW_ID,
            inputs: {
              config: validatedInputs.config,
              samples: validatedInputs.samples,
              options: validatedInputs.options,
            },
            status: getRunCalls === 1 ? 'validating' : 'planned',
            created_at: '2026-07-04T12:00:00.000Z',
            updated_at: '2026-07-04T12:01:00.000Z',
            started_at: null,
            ended_at: null,
            current_stage: 'preflight',
            cancellation_reason: null,
            error: null,
            tags: {},
          },
          issues: [],
        };
      }),
      listRunEvents: vi.fn().mockResolvedValue({
        ok: true,
        run_id: 'run-polling',
        events: [],
        next_cursor: null,
        issues: [],
      }),
      listRunLogs: vi.fn().mockImplementation(async (_, options = {}) => ({
        ok: true,
        run_id: 'run-polling',
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      })),
    };

    renderPanel({ runId: 'run-polling', runClient });

    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent(
      'validating',
    );
    await waitFor(
      () => {
        expect(screen.getByTestId('run-status-badge')).toHaveTextContent(
          'planned',
        );
      },
      { timeout: 2500 },
    );

    const callsAtTerminalStatus = vi.mocked(runClient.getRun).mock.calls.length;
    await new Promise((resolve) => setTimeout(resolve, 1100));
    expect(runClient.getRun).toHaveBeenCalledTimes(callsAtTerminalStatus);
  }, 4000);

  it('pauses preflight polling while cancellation is in flight', async () => {
    const baseClient = createMockRunClient();
    let status = 'validating';
    const getRun = vi.fn(async () => ({
      ok: true,
      run: {
        run_id: 'run-cancelling',
        workflow_id: WORKFLOW_ID,
        inputs: {
          config: validatedInputs.config,
          samples: validatedInputs.samples,
          options: validatedInputs.options,
        },
        status,
        created_at: '2026-07-04T12:00:00.000Z',
        updated_at: '2026-07-04T12:01:00.000Z',
        started_at: null,
        ended_at: null,
        current_stage: 'preflight',
        cancellation_reason:
          status === 'cancelled' ? 'User requested cancellation.' : null,
        error: null,
        tags: {},
      },
      issues: [],
    }));
    const runClient: RunApiClient = {
      ...baseClient,
      getRun,
      listRunEvents: vi.fn().mockResolvedValue({
        ok: true,
        run_id: 'run-cancelling',
        events: [],
        next_cursor: null,
        issues: [],
      }),
      listRunLogs: vi.fn().mockImplementation(async (_, options = {}) => ({
        ok: true,
        run_id: 'run-cancelling',
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      })),
      cancelRun: vi.fn(async () => {
        await new Promise((resolve) => setTimeout(resolve, 1200));
        status = 'cancelled';
        return getRun();
      }),
    };
    const user = userEvent.setup();

    renderPanel({ runId: 'run-cancelling', runClient });
    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent(
      'validating',
    );
    await user.click(screen.getByTestId('cancel-run-button'));

    await waitFor(
      () => {
        expect(screen.getByTestId('run-status-badge')).toHaveTextContent(
          'cancelled',
        );
      },
      { timeout: 2500 },
    );
    expect(runClient.cancelRun).toHaveBeenCalledTimes(1);
  }, 4000);

  it('starts only a planned run through the real start client operation', async () => {
    const runClient = createMockRunClient();
    await runClient.createRun(WORKFLOW_ID, {
      config: validatedInputs.config,
      samples: validatedInputs.samples as Record<string, string>[],
      options: validatedInputs.options,
    });
    await runClient.preflightRun('run-1');
    const user = userEvent.setup();

    renderPanel({ runId: 'run-1', runClient });
    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent('planned');
    await user.click(screen.getByRole('button', { name: 'Start run' }));

    await waitFor(() => expect(runClient.startRun).toHaveBeenCalledWith('run-1'));
    expect(screen.getByTestId('run-status-badge')).toHaveTextContent('queued');
    expect(screen.queryByRole('button', { name: 'Start run' })).not.toBeInTheDocument();
  });

  it('treats an unconfirmed start as unknown and clears it after canonical advancement', async () => {
    const baseClient = createMockRunClient();
    await baseClient.createRun(WORKFLOW_ID, {
      config: validatedInputs.config,
      samples: validatedInputs.samples as Record<string, string>[],
      options: validatedInputs.options,
    });
    await baseClient.preflightRun('run-1');
    const originalGetRun = baseClient.getRun;
    let status = 'planned';
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => {
        const response = await originalGetRun('run-1');
        return {
          ...response,
          run: response.run ? { ...response.run, status } : null,
        } as RunResponse;
      }),
      startRun: vi.fn().mockRejectedValue(new Error('response lost')),
    };
    const user = userEvent.setup();

    renderPanel({ runId: 'run-1', runClient });
    await user.click(await screen.findByRole('button', { name: 'Start run' }));
    expect(
      await screen.findByText(/could not confirm whether the run was submitted/i),
    ).toBeInTheDocument();

    status = 'queued';
    await user.click(screen.getByRole('button', { name: 'Refresh run progress' }));
    await waitFor(() => {
      expect(screen.getByTestId('run-status-badge')).toHaveTextContent('queued');
      expect(
        screen.queryByText(/could not confirm whether the run was submitted/i),
      ).not.toBeInTheDocument();
    });
  });

  it('refetches canonical terminal state after a start conflict envelope', async () => {
    const baseClient = createMockRunClient();
    await baseClient.createRun(WORKFLOW_ID, {
      config: validatedInputs.config,
      samples: validatedInputs.samples as Record<string, string>[],
      options: validatedInputs.options,
    });
    await baseClient.preflightRun('run-1');
    const originalGetRun = baseClient.getRun;
    let status = 'planned';
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => {
        const response = await originalGetRun('run-1');
        return {
          ...response,
          run: response.run
            ? {
                ...response.run,
                status,
                ended_at:
                  status === 'cancelled' ? '2026-07-04T12:02:00.000Z' : null,
                cancellation_reason:
                  status === 'cancelled' ? 'User requested cancellation.' : null,
              }
            : null,
        } as RunResponse;
      }),
      startRun: vi.fn(async () => {
        status = 'cancelled';
        return {
          ok: false,
          run: null,
          issues: [
            {
              code: 'RUN_START_CONFLICT',
              message: 'Run can no longer be started.',
              severity: 'error',
              path: 'run.status',
              source: 'api',
              technical_message: null,
              hint: 'Refresh the run status.',
              context: {},
            },
          ],
        } as RunResponse;
      }),
    };
    const user = userEvent.setup();

    renderPanel({ runId: 'run-1', runClient });
    await user.click(await screen.findByRole('button', { name: 'Start run' }));

    await waitFor(() => {
      expect(screen.getByTestId('run-status-badge')).toHaveTextContent('cancelled');
      expect(screen.queryByRole('button', { name: 'Start run' })).not.toBeInTheDocument();
      expect(screen.queryByText(/run can no longer be started/i)).not.toBeInTheDocument();
    });
    expect(runClient.getRun).toHaveBeenCalledTimes(2);
  });

  it('keeps a running run RUNNING after cancellation is only requested', async () => {
    const baseClient = createMockRunClient();
    const runningRun = {
      run_id: 'run-running',
      workflow_id: WORKFLOW_ID,
      inputs: validatedInputs as unknown as Record<string, unknown>,
      status: 'running',
      created_at: '2026-07-04T12:00:00.000Z',
      updated_at: '2026-07-04T12:01:00.000Z',
      started_at: '2026-07-04T12:01:00.000Z',
      ended_at: null,
      current_stage: 'execution',
      cancellation_reason: null,
      error: null,
      tags: {},
    };
    let requested = false;
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => ({ ok: true, run: runningRun, issues: [] })),
      listRunEvents: vi.fn(async () => ({
        ok: true,
        run_id: runningRun.run_id,
        events: requested
          ? [
              {
                event_id: 'evt-cancel-requested',
                run_id: runningRun.run_id,
                sequence: 1,
                event_type: 'cancellation_requested',
                timestamp: runningRun.updated_at,
                status: 'running',
                stage: 'execution',
                message: 'Cancellation requested.',
                context: {},
                issue: null,
              },
            ]
          : [],
        next_cursor: null,
        issues: [],
      })),
      listRunLogs: vi.fn(async (_, options = {}) => ({
        ok: true,
        run_id: runningRun.run_id,
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      })),
      cancelRun: vi.fn(async () => {
        requested = true;
        return { ok: true, run: runningRun, issues: [] };
      }),
    };
    const user = userEvent.setup();

    renderPanel({ runId: runningRun.run_id, runClient });
    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent('running');
    await user.click(screen.getByRole('button', { name: 'Cancel run' }));

    expect(await screen.findByTestId('cancellation-requested')).toHaveTextContent(
      /remains running until the worker confirms/i,
    );
    expect(screen.getByTestId('run-status-badge')).toHaveTextContent('running');
    expect(screen.getByRole('button', { name: 'Retry cancellation' })).toBeEnabled();
    expect(screen.queryByText('cancelled')).not.toBeInTheDocument();
  });

  it('treats an unconfirmed cancellation as unknown and clears it at canonical terminal', async () => {
    const baseClient = createMockRunClient();
    await baseClient.createRun(WORKFLOW_ID, {
      config: validatedInputs.config,
      samples: validatedInputs.samples as Record<string, string>[],
      options: validatedInputs.options,
    });
    const originalGetRun = baseClient.getRun;
    let status = 'running';
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => {
        const response = await originalGetRun('run-1');
        return {
          ...response,
          run: response.run
            ? {
                ...response.run,
                status,
                started_at: '2026-07-04T12:01:00.000Z',
                ended_at:
                  status === 'cancelled' ? '2026-07-04T12:02:00.000Z' : null,
                cancellation_reason:
                  status === 'cancelled' ? 'User requested cancellation.' : null,
              }
            : null,
        } as RunResponse;
      }),
      cancelRun: vi.fn().mockRejectedValue(new Error('response lost')),
    };
    const user = userEvent.setup();

    renderPanel({ runId: 'run-1', runClient });
    await user.click(await screen.findByRole('button', { name: 'Cancel run' }));
    expect(await screen.findByText(/could not confirm cancellation/i)).toBeInTheDocument();

    status = 'cancelled';
    await user.click(screen.getByRole('button', { name: 'Refresh run progress' }));
    await waitFor(() => {
      expect(screen.getByTestId('run-status-badge')).toHaveTextContent('cancelled');
      expect(screen.queryByText(/could not confirm cancellation/i)).not.toBeInTheDocument();
    });
  });

  it('restores cancellation-requested feedback from persisted events after reload', async () => {
    const baseClient = createMockRunClient();
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => ({
        ok: true,
        run: {
          run_id: 'run-reloaded',
          workflow_id: WORKFLOW_ID,
          inputs: {},
          status: 'running',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:01:00.000Z',
          started_at: '2026-07-04T12:01:00.000Z',
          ended_at: null,
          current_stage: 'execution',
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      })),
      listRunEvents: vi.fn(async () => ({
        ok: true,
        run_id: 'run-reloaded',
        events: [
          {
            event_id: 'evt-request',
            run_id: 'run-reloaded',
            sequence: 2,
            event_type: 'cancellation_requested',
            timestamp: '2026-07-04T12:01:00.000Z',
            status: 'running',
            stage: 'execution',
            message: 'Cancellation requested.',
            context: { internal_path: '/private/worker' },
            issue: null,
          },
        ],
        next_cursor: null,
        issues: [],
      })),
      listRunLogs: vi.fn(async (_, options = {}) => ({
        ok: true,
        run_id: 'run-reloaded',
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      })),
    };

    renderPanel({ runId: 'run-reloaded', runClient });

    expect(await screen.findByTestId('cancellation-requested')).toBeInTheDocument();
    expect(screen.queryByText(/private\/worker/)).not.toBeInTheDocument();
  });

  it('renders an unknown backend status neutrally without unsafe actions', async () => {
    const baseClient = createMockRunClient();
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => ({
        ok: true,
        run: {
          run_id: 'run-future',
          workflow_id: WORKFLOW_ID,
          inputs: {},
          status: 'future_state',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:00:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: null,
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      })),
      listRunEvents: vi.fn(async () => ({
        ok: true,
        run_id: 'run-future',
        events: [],
        next_cursor: null,
        issues: [],
      })),
      listRunLogs: vi.fn(async (_, options = {}) => ({
        ok: true,
        run_id: 'run-future',
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      })),
    };

    renderPanel({ runId: 'run-future', runClient });

    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent('unknown');
    expect(screen.queryByTestId('start-run-button')).not.toBeInTheDocument();
    expect(screen.queryByTestId('cancel-run-button')).not.toBeInTheDocument();
    expect(screen.getByTestId('refresh-run-button')).toBeEnabled();
  });

  it('follows event cursors but stops after a bounded sequence', async () => {
    const baseClient = createMockRunClient();
    const listRunEvents = vi
      .fn()
      .mockResolvedValueOnce({
        ok: true,
        run_id: 'run-pages',
        events: [],
        next_cursor: 'cursor-1',
        issues: [],
      })
      .mockResolvedValueOnce({
        ok: true,
        run_id: 'run-pages',
        events: [],
        next_cursor: null,
        issues: [],
      });
    const runClient: RunApiClient = {
      ...baseClient,
      getRun: vi.fn(async () => ({
        ok: true,
        run: {
          run_id: 'run-pages',
          workflow_id: WORKFLOW_ID,
          inputs: {},
          status: 'planned',
          created_at: '2026-07-04T12:00:00.000Z',
          updated_at: '2026-07-04T12:00:00.000Z',
          started_at: null,
          ended_at: null,
          current_stage: 'preflight',
          cancellation_reason: null,
          error: null,
          tags: {},
        },
        issues: [],
      })),
      listRunEvents,
      listRunLogs: vi.fn(async (_, options = {}) => ({
        ok: true,
        run_id: 'run-pages',
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      })),
    };

    renderPanel({ runId: 'run-pages', runClient });
    await screen.findByTestId('run-status-badge');

    expect(listRunEvents).toHaveBeenNthCalledWith(1, 'run-pages', {
      after: undefined,
      limit: 100,
    });
    expect(listRunEvents).toHaveBeenNthCalledWith(2, 'run-pages', {
      after: 'cursor-1',
      limit: 100,
    });
  });

  it('does not render execution-like wording', () => {
    renderPanel({
      validationResult: successfulValidation,
      validatedInputs,
    });

    const panel = screen.getByTestId('run-progress-panel');
    const text = panel.textContent?.toLowerCase() ?? '';
    const forbidden = [
      'execute',
      'execution',
      'snakemake',
      'subprocess',
      'pipeline',
      'launch',
      'submit',
    ];
    for (const word of forbidden) {
      expect(text).not.toContain(word);
    }
  });
});
