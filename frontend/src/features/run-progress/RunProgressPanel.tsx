import { useEffect, useRef, useState } from 'react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { Ban, Play, RefreshCw } from 'lucide-react';
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
import { RunIssuePanel } from './RunIssuePanel';
import { RunLogPanel } from './RunLogPanel';
import { RunStatusBadge } from './RunStatusBadge';

interface RunProgressPanelProps {
  workflowId?: string | null;
  validationResult?: ValidateWorkflowResponse | null;
  validatedInputs?: WorkflowInputs | null;
  runClient: RunApiClient;
  runId?: string | null;
  onRunCreated?: (runId: string) => void;
  beginPreflight?: boolean;
  preflightRequestId?: string | null;
  onPreflightConsumed?: () => void;
}

const pollingStatuses = new Set(['created', 'validating', 'queued', 'running']);
const cancellableStatuses = new Set([
  'created',
  'validating',
  'planned',
  'queued',
  'running',
]);
const POLL_INTERVAL_MS = 1000;
const MAX_POLL_DURATION_MS = 15 * 60 * 1000;
const PAGE_LIMIT = 100;
const MAX_PAGES = 5;
const claimedAutoPreflightRequests = new Set<string>();

export interface RunSnapshot {
  run: RunRecordResponse | null;
  events: RunEventResponse[];
  stdoutLogs: RunLogChunkResponse[];
  stderrLogs: RunLogChunkResponse[];
  issues: Issue[];
  truncated: boolean;
}

export function runProgressQueryKey(runId: string | null) {
  return ['run-progress', runId] as const;
}

function safeUiIssue(code: string, message: string): Issue {
  return {
    code,
    message,
    severity: 'error',
    path: null,
    source: 'ui',
    technical_message: null,
    hint: null,
    context: {},
  };
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
  if (samples === null || typeof samples === 'string') return samples;
  if (samples.every(isStringRecord)) {
    return samples.map((sample) => ({ ...sample }));
  }
  throw new Error('incompatible samples');
}

function toRunCreateRequest(inputs: WorkflowInputs): RunCreateRequest {
  return {
    config: inputs.config,
    samples: normalizeSamplesForRun(inputs.samples),
    options: inputs.options,
  };
}

function emptySnapshot(run: RunRecordResponse): RunSnapshot {
  return {
    run,
    events: [],
    stdoutLogs: [],
    stderrLogs: [],
    issues: [],
    truncated: false,
  };
}

async function loadEventPages(
  runClient: RunApiClient,
  runId: string,
): Promise<{ values: RunEventResponse[]; issues: Issue[]; truncated: boolean }> {
  const values: RunEventResponse[] = [];
  const issues: Issue[] = [];
  const cursors = new Set<string>();
  let after: string | undefined;

  for (let page = 0; page < MAX_PAGES; page += 1) {
    const response = await runClient.listRunEvents(runId, {
      after,
      limit: PAGE_LIMIT,
    });
    if (!response.ok) {
      issues.push(...response.issues);
      return { values, issues, truncated: false };
    }
    values.push(...response.events);
    if (!response.next_cursor) return { values, issues, truncated: false };
    if (cursors.has(response.next_cursor)) {
      issues.push(
        safeUiIssue('RUN_CURSOR_INVALID', 'Run event pagination could not continue safely.'),
      );
      return { values, issues, truncated: true };
    }
    cursors.add(response.next_cursor);
    after = response.next_cursor;
  }

  return { values, issues, truncated: true };
}

async function loadLogPages(
  runClient: RunApiClient,
  runId: string,
  streamName: 'stdout' | 'stderr',
): Promise<{ values: RunLogChunkResponse[]; issues: Issue[]; truncated: boolean }> {
  const values: RunLogChunkResponse[] = [];
  const issues: Issue[] = [];
  const cursors = new Set<string>();
  let after: string | undefined;

  for (let page = 0; page < MAX_PAGES; page += 1) {
    const response = await runClient.listRunLogs(runId, {
      streamName,
      after,
      limit: PAGE_LIMIT,
    });
    if (!response.ok) {
      issues.push(...response.issues);
      return { values, issues, truncated: false };
    }
    values.push(...response.chunks);
    if (!response.next_cursor) return { values, issues, truncated: false };
    if (cursors.has(response.next_cursor)) {
      issues.push(
        safeUiIssue('RUN_CURSOR_INVALID', 'Run log pagination could not continue safely.'),
      );
      return { values, issues, truncated: true };
    }
    cursors.add(response.next_cursor);
    after = response.next_cursor;
  }

  return { values, issues, truncated: true };
}

export async function loadRunSnapshot(
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
      issues: runResponse.issues,
      truncated: false,
    };
  }
  if (!runResponse.run) {
    return {
      run: null,
      events: [],
      stdoutLogs: [],
      stderrLogs: [],
      issues: [safeUiIssue('RUN_RECORD_MISSING', 'The run record was missing.')],
      truncated: false,
    };
  }

  const [events, stdout, stderr] = await Promise.all([
    loadEventPages(runClient, runId),
    loadLogPages(runClient, runId, 'stdout'),
    loadLogPages(runClient, runId, 'stderr'),
  ]);
  return {
    run: runResponse.run,
    events: events.values,
    stdoutLogs: stdout.values,
    stderrLogs: stderr.values,
    issues: [...events.issues, ...stdout.issues, ...stderr.issues],
    truncated: events.truncated || stdout.truncated || stderr.truncated,
  };
}

function mergeRun(current: RunSnapshot | undefined, run: RunRecordResponse): RunSnapshot {
  return current ? { ...current, run } : emptySnapshot(run);
}

export function RunProgressPanel({
  workflowId = null,
  validationResult = null,
  validatedInputs = null,
  runClient,
  runId = null,
  onRunCreated,
  beginPreflight = false,
  preflightRequestId = null,
  onPreflightConsumed,
}: RunProgressPanelProps) {
  const [localRunId, setLocalRunId] = useState<string | null>(null);
  const [activeLogStream, setActiveLogStream] = useState<'stdout' | 'stderr'>('stdout');
  const [actionIssues, setActionIssues] = useState<Issue[]>([]);
  const [cancellationRequested, setCancellationRequested] = useState(false);
  const preflightConsumedFor = useRef<string | null>(null);
  const pollStartedAt = useRef(Date.now());
  const queryClient = useQueryClient();
  const targetRunId = runId ?? localRunId;

  const runQuery = useQuery({
    queryKey: runProgressQueryKey(targetRunId),
    queryFn: () => loadRunSnapshot(runClient, targetRunId!),
    enabled: targetRunId !== null,
    refetchInterval: (query) => {
      const status = query.state.data?.run?.status;
      if (!status || !pollingStatuses.has(status)) return false;
      if (Date.now() - pollStartedAt.current >= MAX_POLL_DURATION_MS) return false;
      return POLL_INTERVAL_MS;
    },
  });

  const updateRun = (run: RunRecordResponse) => {
    queryClient.setQueryData<RunSnapshot>(
      runProgressQueryKey(run.run_id),
      (current) => mergeRun(current, run),
    );
  };

  const preflightMutation = useMutation({
    mutationFn: (id: string) => runClient.preflightRun(id),
    onSuccess: (response) => {
      if (!response.ok) {
        setActionIssues(response.issues);
        return;
      }
      setActionIssues([]);
      if (response.run) updateRun(response.run);
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
    onError: () => {
      setActionIssues([
        safeUiIssue('PREFLIGHT_UNAVAILABLE', 'Preflight could not be started. Try again.'),
      ]);
    },
  });

  const startMutation = useMutation({
    mutationFn: (id: string) => runClient.startRun(id),
    onSuccess: (response) => {
      if (!response.ok) {
        setActionIssues(response.issues);
        return;
      }
      setActionIssues([]);
      if (response.run) updateRun(response.run);
      pollStartedAt.current = Date.now();
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
    onError: () => {
      setActionIssues([
        safeUiIssue('RUN_START_UNAVAILABLE', 'The run could not be started. Try again.'),
      ]);
    },
  });

  const cancelMutation = useMutation({
    mutationFn: (id: string) => runClient.cancelRun(id),
    onSuccess: (response) => {
      if (!response.ok) {
        setActionIssues(response.issues);
        void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
        return;
      }
      setActionIssues([]);
      if (response.run) {
        updateRun(response.run);
        setCancellationRequested(response.run.status === 'running');
      }
      pollStartedAt.current = Date.now();
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
    onError: () => {
      setActionIssues([
        safeUiIssue(
          'RUN_CANCELLATION_UNAVAILABLE',
          'Cancellation could not be requested. The run state was not changed; try again.',
        ),
      ]);
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
  });

  const createMutation = useMutation({
    mutationFn: async () => {
      if (!workflowId || !validatedInputs) {
        throw new Error('invalid create inputs');
      }
      return runClient.createRun(workflowId, toRunCreateRequest(validatedInputs));
    },
  });

  useEffect(() => {
    setLocalRunId(null);
    setActiveLogStream('stdout');
    setActionIssues([]);
    setCancellationRequested(false);
    preflightConsumedFor.current = null;
    pollStartedAt.current = Date.now();
  }, [runId, runClient, workflowId, validatedInputs]);

  useEffect(() => {
    setActiveLogStream('stdout');
  }, [targetRunId]);

  useEffect(() => {
    if (
      beginPreflight &&
      targetRunId &&
      preflightRequestId &&
      preflightConsumedFor.current !== targetRunId &&
      !claimedAutoPreflightRequests.has(preflightRequestId)
    ) {
      preflightConsumedFor.current = targetRunId;
      claimedAutoPreflightRequests.add(preflightRequestId);
      preflightMutation.mutate(targetRunId);
      onPreflightConsumed?.();
    }
  }, [
    beginPreflight,
    onPreflightConsumed,
    preflightMutation,
    preflightRequestId,
    targetRunId,
  ]);

  const snapshot = runQuery.data;
  const run = snapshot?.run ?? null;
  const cancellationEventPresent =
    snapshot?.events.some((event) => event.event_type === 'cancellation_requested') ?? false;
  const showCancellationRequested =
    run?.status === 'running' && (cancellationRequested || cancellationEventPresent);

  useEffect(() => {
    if (run && run.status !== 'running') setCancellationRequested(false);
  }, [run]);

  const canCreateRun =
    workflowId !== null && validationResult?.ok === true && validatedInputs !== null;
  const executionActionPending = startMutation.isPending || cancelMutation.isPending;

  async function handleCreateRun() {
    setActionIssues([]);
    try {
      const response = await createMutation.mutateAsync();
      if (!response.ok) {
        setActionIssues(response.issues);
        return;
      }
      if (!response.run) {
        setActionIssues([safeUiIssue('RUN_RECORD_MISSING', 'The run record was missing.')]);
        return;
      }
      updateRun(response.run);
      setActiveLogStream('stdout');
      if (onRunCreated) {
        onRunCreated(response.run.run_id);
        return;
      }
      setLocalRunId(response.run.run_id);
      await preflightMutation.mutateAsync(response.run.run_id);
    } catch {
      setActionIssues([
        safeUiIssue('RUN_CREATE_UNAVAILABLE', 'The run record could not be created. Try again.'),
      ]);
    }
  }

  function handleRefresh() {
    setActionIssues([]);
    pollStartedAt.current = Date.now();
    void runQuery.refetch();
  }

  const queryIssues = snapshot?.issues ?? [];
  const unexpectedQueryIssues = runQuery.error
    ? [safeUiIssue('RUN_PROGRESS_UNAVAILABLE', 'Run progress could not be loaded. Try again.')]
    : [];
  const visibleIssues = [...actionIssues, ...unexpectedQueryIssues, ...queryIssues];
  const canCancelRun = run !== null && cancellableStatuses.has(run.status);
  const canStartRun = run?.status === 'planned';
  const canPreflightRun = run?.status === 'created';
  const showNoRunPlaceholder = !run && !runQuery.isLoading && visibleIssues.length === 0 && !runId;
  const showMissingRunError =
    !run && !runQuery.isLoading && visibleIssues.length > 0 && runId !== null;

  return (
    <div className="space-y-4" data-testid="run-progress-panel">
      {showNoRunPlaceholder && (
        <p className="text-sm text-[var(--color-text-muted)]">
          Validate inputs before creating a run record.
        </p>
      )}

      {runQuery.isLoading && (
        <p className="text-sm text-[var(--color-text-muted)]" data-testid="run-loading">
          Loading run…
        </p>
      )}

      <RunIssuePanel issues={visibleIssues} />

      {showMissingRunError && (
        <p className="text-sm text-[var(--color-text-muted)]" data-testid="run-missing-hint">
          Run could not be loaded.
        </p>
      )}

      {run && (
        <div className="min-w-0 space-y-3">
          <div className="flex min-w-0 flex-wrap items-center gap-3 text-sm">
            <span className="font-medium text-[var(--color-text)]">Run ID:</span>
            <code className="min-w-0 break-all text-xs text-[var(--color-text-muted)]">
              {run.run_id}
            </code>
            <RunStatusBadge status={run.status} />
          </div>
          <div className="grid grid-cols-1 gap-2 text-xs text-[var(--color-text-muted)] sm:grid-cols-2">
            <div>Created: {new Date(run.created_at).toLocaleString()}</div>
            <div>Updated: {new Date(run.updated_at).toLocaleString()}</div>
          </div>

          {showCancellationRequested && (
            <div
              className="rounded border border-yellow-300 bg-yellow-50 p-2 text-sm text-yellow-800"
              role="status"
              data-testid="cancellation-requested"
            >
              Cancellation requested. The run remains running until the worker confirms termination.
            </div>
          )}

          {run.error && <RunIssuePanel issues={[run.error]} title="Run failure" />}

          {snapshot?.truncated && (
            <p className="text-xs text-[var(--color-text-muted)]" role="status">
              Showing the first {MAX_PAGES * PAGE_LIMIT} entries per stream. Refresh to load the latest bounded view.
            </p>
          )}

          <div>
            <h4 className="mb-2 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
              Events
            </h4>
            <RunEventFeed events={snapshot?.events ?? []} />
          </div>

          <div>
            <h4 className="mb-2 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
              Logs
            </h4>
            <RunLogPanel
              stdoutChunks={snapshot?.stdoutLogs ?? []}
              stderrChunks={snapshot?.stderrLogs ?? []}
              activeStream={activeLogStream}
              onStreamChange={setActiveLogStream}
            />
          </div>
        </div>
      )}

      <div className="flex flex-wrap gap-2">
        {!runId && (
          <Button
            variant="primary"
            onClick={handleCreateRun}
            disabled={!canCreateRun || createMutation.isPending}
            aria-label="Create run"
            data-testid="create-run-button"
          >
            {createMutation.isPending ? 'Creating run…' : 'Create run'}
          </Button>
        )}
        {run && (
          <>
            {canPreflightRun && (
              <Button
                variant="primary"
                onClick={() => preflightMutation.mutate(run.run_id)}
                disabled={preflightMutation.isPending}
                aria-label="Run preflight"
                data-testid="preflight-run-button"
              >
                {preflightMutation.isPending ? 'Starting preflight…' : 'Run preflight'}
              </Button>
            )}
            {canStartRun && (
              <Button
                variant="primary"
                onClick={() => startMutation.mutate(run.run_id)}
                disabled={executionActionPending}
                aria-label="Start run"
                data-testid="start-run-button"
              >
                <Play className="mr-1.5 h-4 w-4" aria-hidden="true" />
                {startMutation.isPending ? 'Starting run…' : 'Start run'}
              </Button>
            )}
            <Button
              variant="secondary"
              onClick={handleRefresh}
              disabled={runQuery.isFetching || executionActionPending}
              aria-label="Refresh run progress"
              data-testid="refresh-run-button"
            >
              <RefreshCw className="mr-1.5 h-4 w-4" aria-hidden="true" />
              Refresh
            </Button>
            {canCancelRun && (
              <Button
                variant="secondary"
                onClick={() => cancelMutation.mutate(run.run_id)}
                disabled={executionActionPending}
                aria-label={showCancellationRequested ? 'Retry cancellation' : 'Cancel run'}
                data-testid="cancel-run-button"
              >
                <Ban className="mr-1.5 h-4 w-4" aria-hidden="true" />
                {cancelMutation.isPending
                  ? 'Requesting cancellation…'
                  : showCancellationRequested
                    ? 'Retry cancellation'
                    : 'Cancel run'}
              </Button>
            )}
          </>
        )}
      </div>
    </div>
  );
}
