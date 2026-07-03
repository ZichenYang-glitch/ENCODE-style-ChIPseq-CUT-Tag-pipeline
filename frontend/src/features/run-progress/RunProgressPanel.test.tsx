import { describe, it, expect, vi } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { RunProgressPanel } from './RunProgressPanel';
import type { RunApiClient } from '../../api/runClient';
import type { ValidateWorkflowResponse, WorkflowInputs } from '../../api/types';
import type { RunEventResponse } from '../../api/runTypes';

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
    listRunLogs: vi.fn().mockResolvedValue({
      ok: true,
      run_id: '',
      stream_name: 'stdout',
      chunks: [],
      next_cursor: null,
      issues: [],
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
  return {
    runClient,
    ...render(<RunProgressPanel {...defaultProps} {...props} />),
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

  it('creates a run record and shows status, events, and empty logs', async () => {
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

    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent('created');
    expect(screen.getByTestId('run-event-feed')).toBeInTheDocument();
    expect(screen.getByText(/status_changed/i)).toBeInTheDocument();
    expect(screen.getByTestId('run-log-panel-empty')).toHaveTextContent(
      /No log entries yet/i,
    );
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
      expect(runClient.listRunLogs).toHaveBeenCalledTimes(2);
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
