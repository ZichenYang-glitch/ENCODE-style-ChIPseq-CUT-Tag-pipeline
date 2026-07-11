import { useEffect, useRef, useState } from 'react';
import { useQuery, useQueryClient } from '@tanstack/react-query';
import { RefreshCw, Ban } from 'lucide-react';
import type { RunApiClient } from '../../api/runClient';
import type {
  RunCreateRequest,
  RunEventResponse,
  RunLogChunkResponse,
  RunRecordResponse,
} from '../../api/runTypes';
import type { Issue, ValidateWorkflowResponse, WorkflowInputs } from '../../api/types';
import { Button } from '../../components/Button';
import { RunEventFeed } from './RunEventFeed';
import { RunLogPanel } from './RunLogPanel';
import { RunStatusBadge } from './RunStatusBadge';

interface RunProgressPanelProps {
  workflowId?: string | null;
  validationResult?: ValidateWorkflowResponse | null;
  validatedInputs?: WorkflowInputs | null;
  runClient: RunApiClient;
  runId?: string | null;
  onRunCreated?: (runId: string) => void;
}

const terminalStatuses = ['succeeded', 'failed', 'cancelled'];
const preflightPollingStatuses = ['validating', 'queued', 'running'];
const PREFLIGHT_POLL_INTERVAL_MS = 1000;

interface RunSnapshot {
  run: RunRecordResponse | null;
  events: RunEventResponse[];
  stdoutLogs: RunLogChunkResponse[];
  stderrLogs: RunLogChunkResponse[];
  issue: string | null;
}

function runProgressQueryKey(runId: string | null) {
  return ['run-progress', runId] as const;
}

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

function firstIssueMessage(issues: Issue[] | undefined): string {
  const issue = issues?.[0];
  if (!issue) return 'Request failed.';
  return `${issue.code}: ${issue.message}`;
}

async function loadRunSnapshot(
  runClient: RunApiClient,
  runId: string,
): Promise<RunSnapshot> {
  const runResponse = await runClient.getRun(runId);
  if (!runResponse.ok) {
    return {
      run: null,
      events: [],
      stdoutLogs: [],
      stderrLogs: [],
      issue: firstIssueMessage(runResponse.issues),
    };
  }

  if (!runResponse.run) {
    return {
      run: null,
      events: [],
      stdoutLogs: [],
      stderrLogs: [],
      issue: 'Run record missing.',
    };
  }

  const [eventsResponse, stdoutResponse, stderrResponse] = await Promise.all([
    runClient.listRunEvents(runId, { limit: 50 }),
    runClient.listRunLogs(runId, { streamName: 'stdout', limit: 50 }),
    runClient.listRunLogs(runId, { streamName: 'stderr', limit: 50 }),
  ]);
  const issue = !eventsResponse.ok
    ? firstIssueMessage(eventsResponse.issues)
    : !stdoutResponse.ok
      ? firstIssueMessage(stdoutResponse.issues)
      : !stderrResponse.ok
        ? firstIssueMessage(stderrResponse.issues)
        : null;

  return {
    run: runResponse.run,
    events: eventsResponse.ok ? eventsResponse.events : [],
    stdoutLogs: stdoutResponse.ok ? stdoutResponse.chunks : [],
    stderrLogs: stderrResponse.ok ? stderrResponse.chunks : [],
    issue,
  };
}

export function RunProgressPanel({
  workflowId = null,
  validationResult = null,
  validatedInputs = null,
  runClient,
  runId = null,
  onRunCreated,
}: RunProgressPanelProps) {
  const [localRunId, setLocalRunId] = useState<string | null>(null);
  const [optimisticRun, setOptimisticRun] =
    useState<RunRecordResponse | null>(null);
  const [activeLogStream, setActiveLogStream] = useState<'stdout' | 'stderr'>('stdout');
  const [actionLoading, setActionLoading] = useState(false);
  const [actionError, setActionError] = useState<string | null>(null);
  const requestGeneration = useRef(0);
  const queryClient = useQueryClient();

  const targetRunId = runId ?? localRunId;
  const runQuery = useQuery({
    queryKey: runProgressQueryKey(targetRunId),
    queryFn: () => loadRunSnapshot(runClient, targetRunId!),
    enabled: targetRunId !== null,
    refetchInterval: (query) => {
      if (actionLoading) return false;
      const status = query.state.data?.run?.status;
      return status && preflightPollingStatuses.includes(status)
        ? PREFLIGHT_POLL_INTERVAL_MS
        : false;
    },
  });

  const snapshot = runQuery.data;
  const run = snapshot?.run ?? optimisticRun;
  const events = snapshot?.events ?? [];
  const stdoutLogs = snapshot?.stdoutLogs ?? [];
  const stderrLogs = snapshot?.stderrLogs ?? [];
  const queryError =
    runQuery.error instanceof Error
      ? runQuery.error.message
      : runQuery.error
        ? 'Failed to load run.'
        : snapshot?.issue ?? null;
  const error = actionError ?? queryError;
  const loading = actionLoading || (targetRunId !== null && runQuery.isLoading);

  function startRequestGeneration(): number {
    requestGeneration.current += 1;
    return requestGeneration.current;
  }

  function isCurrentRequest(generation: number): boolean {
    return requestGeneration.current === generation;
  }

  useEffect(() => {
    startRequestGeneration();
    setLocalRunId(null);
    setOptimisticRun(null);
    setActiveLogStream('stdout');
    setActionError(null);
    setActionLoading(false);
  }, [runId, runClient, workflowId, validatedInputs]);

  useEffect(() => {
    setActiveLogStream('stdout');
  }, [targetRunId]);

  const canCreateRun =
    workflowId !== null &&
    validationResult?.ok === true &&
    validatedInputs !== null;

  async function handleCreateRun() {
    if (!canCreateRun || !workflowId || !validatedInputs) return;

    const generation = startRequestGeneration();
    setActionLoading(true);
    setActionError(null);

    try {
      const response = await runClient.createRun(
        workflowId,
        toRunCreateRequest(validatedInputs),
      );

      if (!isCurrentRequest(generation)) return;

      if (!response.ok) {
        setActionError(firstIssueMessage(response.issues));
        return;
      }

      const createdRun = response.run;
      if (!createdRun) {
        setActionError('Run record missing.');
        return;
      }

      setOptimisticRun(createdRun);
      setActiveLogStream('stdout');
      const preflightResponse = await runClient.preflightRun(createdRun.run_id);
      if (!isCurrentRequest(generation)) return;
      if (!preflightResponse.ok) {
        setActionError(firstIssueMessage(preflightResponse.issues));
        return;
      }
      if (preflightResponse.run) setOptimisticRun(preflightResponse.run);
      setLocalRunId(createdRun.run_id);
      if (onRunCreated) {
        onRunCreated(createdRun.run_id);
      }
    } catch (err) {
      if (!isCurrentRequest(generation)) return;
      setActionError(
        err instanceof Error ? err.message : 'Failed to create run record.',
      );
    } finally {
      if (isCurrentRequest(generation)) setActionLoading(false);
    }
  }

  async function handleCancelRun() {
    if (!run) return;

    const generation = startRequestGeneration();
    setActionLoading(true);
    setActionError(null);

    try {
      const response = await runClient.cancelRun(run.run_id);
      if (!isCurrentRequest(generation)) return;
      if (!response.ok) {
        setActionError(firstIssueMessage(response.issues));
        return;
      }
      if (response.run) {
        const cancelledRun = response.run;
        queryClient.setQueryData<RunSnapshot>(
          runProgressQueryKey(cancelledRun.run_id),
          (current) => ({
            run: cancelledRun,
            events: current?.events ?? [],
            stdoutLogs: current?.stdoutLogs ?? [],
            stderrLogs: current?.stderrLogs ?? [],
            issue: current?.issue ?? null,
          }),
        );
      }
      await runQuery.refetch();
    } catch (err) {
      if (!isCurrentRequest(generation)) return;
      setActionError(
        err instanceof Error ? err.message : 'Failed to cancel run.',
      );
    } finally {
      if (isCurrentRequest(generation)) setActionLoading(false);
    }
  }

  async function handleRefresh() {
    if (!run) return;

    const generation = startRequestGeneration();
    setActionLoading(true);
    setActionError(null);

    try {
      await runQuery.refetch();
    } catch (err) {
      if (!isCurrentRequest(generation)) return;
      setActionError(
        err instanceof Error ? err.message : 'Failed to refresh run progress.',
      );
    } finally {
      if (isCurrentRequest(generation)) setActionLoading(false);
    }
  }

  const canCancelRun = run !== null && !terminalStatuses.includes(run.status);
  const canRefresh = run !== null;

  const showNoRunPlaceholder = !run && !loading && !error && !runId;
  const showMissingRunError = !run && !loading && error && runId !== null;

  return (
    <div className="space-y-4" data-testid="run-progress-panel">
      {showNoRunPlaceholder && (
        <p className="text-sm text-[var(--color-text-muted)]">
          Validate inputs before creating a run record.
        </p>
      )}

      {loading && (
        <p className="text-sm text-[var(--color-text-muted)]" data-testid="run-loading">
          Loading run…
        </p>
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

      {showMissingRunError && (
        <p className="text-sm text-[var(--color-text-muted)]" data-testid="run-missing-hint">
          Run could not be loaded.
        </p>
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
          {loading && !run ? 'Starting preflight…' : 'Start preflight'}
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
