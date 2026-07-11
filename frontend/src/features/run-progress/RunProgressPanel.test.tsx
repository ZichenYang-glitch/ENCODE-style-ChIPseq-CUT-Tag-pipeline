import { describe, it, expect, vi, type Mock } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import userEvent from '@testing-library/user-event';
import { RunProgressPanel } from './RunProgressPanel';
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
      'Run record missing.',
    );
    expect(runClient.preflightRun).not.toHaveBeenCalled();
  });

  it('rejects validated samples that are not compatible with run creation', async () => {
    const user = userEvent.setup();
    const invalidValidatedInputs: WorkflowInputs = {
      ...validatedInputs,
      samples: [{ name: 'sample-1', replicate: 1 }],
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
      /Validated samples are not compatible with run creation/i,
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
      expect(runClient.getRun).toHaveBeenCalledTimes(2);
      expect(runClient.listRunEvents).toHaveBeenCalledTimes(2);
      expect(runClient.listRunLogs).toHaveBeenCalledTimes(4); // 2 streams x 2 refreshes
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
      expect(runClient.listRunLogs).toHaveBeenCalledTimes(2);
    });

    const stdoutCalls = (runClient.listRunLogs as Mock).mock.calls.filter(
      ([, options]) => options?.streamName === 'stdout',
    );
    const stderrCalls = (runClient.listRunLogs as Mock).mock.calls.filter(
      ([, options]) => options?.streamName === 'stderr',
    );

    expect(stdoutCalls).toHaveLength(1);
    expect(stderrCalls).toHaveLength(1);
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
