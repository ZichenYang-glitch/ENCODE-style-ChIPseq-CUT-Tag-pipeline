import { useEffect, useState } from 'react';
import { useQuery } from '@tanstack/react-query';
import { RefreshCw, Ban } from 'lucide-react';
import { useClients } from '../../api/client-context';
import type { RunApiClient } from '../../api/runClient';
import type {
  RunCreateRequest,
  RunEventResponse,
  RunLogChunkResponse,
  RunRecordResponse,
} from '../../api/runTypes';
import type { ValidateWorkflowResponse, WorkflowInputs } from '../../api/types';
import { Button } from '../../components/Button';
import { RunEventFeed } from './RunEventFeed';
import { RunLogPanel } from './RunLogPanel';
import { RunStatusBadge } from './RunStatusBadge';

interface RunProgressPanelProps {
  workflowId?: string | null;
  validationResult?: ValidateWorkflowResponse | null;
  validatedInputs?: WorkflowInputs | null;
  runClient?: RunApiClient;
  runId?: string | null;
  onRunCreated?: (runId: string) => void;
}

const terminalStatuses = ['succeeded', 'failed', 'cancelled'];

function isStringRecord(value: unknown): value is Record<string, string> {
  return (
    typeof value === 'object' &&
    value !== null &&
    !Array.isArray(value) &&
    Object.values(value).every((entry) => typeof entry === 'string')
  );
}

function normalizeSamplesForRun(
  samples: WorkflowInputs['samples'],
): RunCreateRequest['samples'] {
  if (samples === null || typeof samples === 'string') {
    return samples;
  }

  if (samples.every(isStringRecord)) {
    return samples.map((sample) => ({ ...sample }));
  }

  throw new Error('Validated samples are not compatible with run creation.');
}

function toRunCreateRequest(inputs: WorkflowInputs): RunCreateRequest {
  return {
    config: inputs.config,
    samples: normalizeSamplesForRun(inputs.samples),
    options: inputs.options,
  };
}

export function RunProgressPanel({
  workflowId = null,
  validationResult = null,
  validatedInputs = null,
  runClient: runClientProp,
  runId = null,
  onRunCreated,
}: RunProgressPanelProps) {
  const { runClient: clientFromContext } = useClients();
  const runClient = runClientProp ?? clientFromContext;

  const [run, setRun] = useState<RunRecordResponse | null>(null);
  const [events, setEvents] = useState<RunEventResponse[]>([]);
  const [stdoutLogs, setStdoutLogs] = useState<RunLogChunkResponse[]>([]);
  const [stderrLogs, setStderrLogs] = useState<RunLogChunkResponse[]>([]);
  const [activeLogStream, setActiveLogStream] = useState<'stdout' | 'stderr'>('stdout');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const {
    data: runQueryData,
    isLoading: runQueryLoading,
    error: runQueryError,
  } = useQuery({
    queryKey: ['run', runId],
    queryFn: async () => {
      if (!runId) return null;
      const response = await runClient.getRun(runId);
      return response;
    },
    enabled: runId !== null,
    retry: false,
  });

  useEffect(() => {
    if (runId && runQueryData?.ok && runQueryData.run) {
      setRun(runQueryData.run);
      void refreshRun(runQueryData.run.run_id);
    }
  }, [runId, runQueryData]);

  useEffect(() => {
    if (runQueryError) {
      setError(
        runQueryError instanceof Error
          ? runQueryError.message
          : 'Failed to load run record.',
      );
    }
  }, [runQueryError]);

  useEffect(() => {
    setRun(null);
    setEvents([]);
    setStdoutLogs([]);
    setStderrLogs([]);
    setActiveLogStream('stdout');
    setError(null);
  }, [workflowId, validatedInputs, runId]);

  const canCreateRun =
    workflowId !== null &&
    validationResult?.ok === true &&
    validatedInputs !== null;

  async function refreshRun(runIdToRefresh: string) {
    const [runResponse, eventsResponse, stdoutResponse, stderrResponse] =
      await Promise.all([
        runClient.getRun(runIdToRefresh),
        runClient.listRunEvents(runIdToRefresh, { limit: 50 }),
        runClient.listRunLogs(runIdToRefresh, { streamName: 'stdout', limit: 50 }),
        runClient.listRunLogs(runIdToRefresh, { streamName: 'stderr', limit: 50 }),
      ]);

    if (runResponse.run) {
      setRun(runResponse.run);
    }
    setEvents(eventsResponse.events);
    setStdoutLogs(stdoutResponse.chunks);
    setStderrLogs(stderrResponse.chunks);
  }

  async function handleCreateRun() {
    if (!canCreateRun || !workflowId || !validatedInputs) return;

    setLoading(true);
    setError(null);

    try {
      const response = await runClient.createRun(
        workflowId,
        toRunCreateRequest(validatedInputs),
      );

      if (!response.ok) {
        const message =
          response.issues[0]?.message ?? 'Failed to create run record.';
        setError(message);
        setLoading(false);
        return;
      }

      const createdRun = response.run;
      if (createdRun) {
        setRun(createdRun);
        setStdoutLogs([]);
        setStderrLogs([]);
        setActiveLogStream('stdout');
        await refreshRun(createdRun.run_id);
        onRunCreated?.(createdRun.run_id);
      }
    } catch (err) {
      setError(
        err instanceof Error ? err.message : 'Failed to create run record.',
      );
    } finally {
      setLoading(false);
    }
  }

  async function handleCancelRun() {
    if (!run) return;

    setLoading(true);
    setError(null);

    try {
      await runClient.cancelRun(run.run_id);
      await refreshRun(run.run_id);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to cancel run.');
    } finally {
      setLoading(false);
    }
  }

  async function handleRefresh() {
    if (!run) return;

    setLoading(true);
    setError(null);

    try {
      await refreshRun(run.run_id);
    } catch (err) {
      setError(
        err instanceof Error ? err.message : 'Failed to refresh run progress.',
      );
    } finally {
      setLoading(false);
    }
  }

  const canCancelRun = run !== null && !terminalStatuses.includes(run.status);
  const canRefresh = run !== null;

  return (
    <div className="space-y-4" data-testid="run-progress-panel">
      {!run && !runQueryLoading && (
        <p className="text-sm text-[var(--color-text-muted)]">
          Validate inputs before creating a run record.
        </p>
      )}

      {runQueryLoading && (
        <p className="text-sm text-[var(--color-text-muted)]">Loading run…</p>
      )}

      {error && (
        <div
          className="rounded border border-[var(--color-error)] bg-[var(--color-error-bg)] p-2 text-xs text-[var(--color-error)]"
          role="alert"
          data-testid="run-progress-error"
        >
          {error}
        </div>
      )}

      {run && (
        <div className="space-y-3">
          <div className="flex flex-wrap items-center gap-3 text-sm">
            <span className="font-medium text-[var(--color-text)]">Run ID:</span>
            <code className="text-xs text-[var(--color-text-muted)]">
              {run.run_id}
            </code>
            <RunStatusBadge status={run.status} />
          </div>
          <div className="grid grid-cols-2 gap-2 text-xs text-[var(--color-text-muted)]">
            <div>Created: {new Date(run.created_at).toLocaleString()}</div>
            <div>Updated: {new Date(run.updated_at).toLocaleString()}</div>
          </div>

          <div>
            <h4 className="mb-2 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
              Events
            </h4>
            <RunEventFeed events={events} />
          </div>

          <div>
            <h4 className="mb-2 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
              Logs
            </h4>
            <RunLogPanel
              stdoutChunks={stdoutLogs}
              stderrChunks={stderrLogs}
              activeStream={activeLogStream}
              onStreamChange={setActiveLogStream}
            />
          </div>
        </div>
      )}

      <div className="flex flex-wrap gap-2">
        <Button
          variant="primary"
          onClick={handleCreateRun}
          disabled={!canCreateRun || loading}
          data-testid="create-run-button"
        >
          {loading && !run ? 'Creating run record…' : 'Create run'}
        </Button>
        {run && (
          <>
            <Button
              variant="secondary"
              onClick={handleRefresh}
              disabled={!canRefresh || loading}
              data-testid="refresh-run-button"
            >
              <RefreshCw className="mr-1.5 h-4 w-4" aria-hidden="true" />
              Refresh
            </Button>
            <Button
              variant="secondary"
              onClick={handleCancelRun}
              disabled={!canCancelRun || loading}
              data-testid="cancel-run-button"
            >
              <Ban className="mr-1.5 h-4 w-4" aria-hidden="true" />
              Cancel run
            </Button>
          </>
        )}
      </div>
    </div>
  );
}
