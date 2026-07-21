import { useEffect, useRef, useState, type KeyboardEvent } from 'react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { Ban, Play, RefreshCw } from 'lucide-react';
import type { WorkflowApiClient } from '../../api/client';
import type { RunApiClient } from '../../api/runClient';
import type {
  RunEventResponse,
  RunLogChunkResponse,
  RunRecordResponse,
} from '../../api/runTypes';
import type {
  Issue,
  ValidateWorkflowResponse,
  WorkflowExecutionAvailability,
} from '../../api/types';
import { Button } from '../../components/Button';
import { ArtifactBrowser } from '../run-artifacts/ArtifactBrowser';
import {
  artifactExtractionOutcome,
  type ArtifactExtractionOutcome,
} from '../run-artifacts/artifactState';
import { QcWorkbench } from '../run-qc/QcWorkbench';
import { ExecutionAvailabilityNotice } from '../workflow-detail/WorkflowAvailability';
import {
  qcIndexingOutcome,
  type QcIndexingOutcome,
} from '../run-qc/qcState';
import { RunEventFeed } from './RunEventFeed';
import { RunIssuePanel } from './RunIssuePanel';
import { RunLogPanel } from './RunLogPanel';
import { RunStatusBadge } from './RunStatusBadge';

export type RunDetailView = 'activity' | 'artifacts' | 'qc';

interface RunProgressPanelProps {
  workflowId?: string | null;
  workflowClient?: WorkflowApiClient;
  executionAvailability?: WorkflowExecutionAvailability | null;
  validationResult?: ValidateWorkflowResponse | null;
  runClient: RunApiClient;
  runId?: string | null;
  onRunCreated?: (runId: string) => void;
  beginPreflight?: boolean;
  preflightRequestId?: string | null;
  onPreflightConsumed?: () => void;
  activeView?: RunDetailView;
  selectedArtifactId?: string | null;
  onViewChange?: (view: RunDetailView) => void;
  onArtifactSelect?: (artifactId: string) => void;
}

const pollingStatuses = new Set(['created', 'validating', 'queued', 'running']);
const terminalStatuses = new Set(['succeeded', 'failed', 'cancelled']);
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

type ActionIssueKind = 'create' | 'preflight' | 'start' | 'cancel' | null;

interface CreateRunAttempt {
  authority: string;
  workflowId: string;
  snapshotId: string;
}

export function isRunPollingPaused(
  status: string | undefined,
  startedAt: number,
  now: number,
): boolean {
  return (
    status !== undefined &&
    pollingStatuses.has(status) &&
    now - startedAt >= MAX_POLL_DURATION_MS
  );
}

export interface RunSnapshot {
  run: RunRecordResponse | null;
  events: RunEventResponse[];
  eventsTruncated: boolean;
  stdoutLogs: RunLogChunkResponse[];
  stderrLogs: RunLogChunkResponse[];
  issues: Issue[];
  truncated: boolean;
}

function extractionOutcome(snapshot: RunSnapshot | undefined): ArtifactExtractionOutcome {
  return artifactExtractionOutcome(
    snapshot?.events ?? [],
    snapshot?.eventsTruncated ?? false,
  );
}

function qcOutcome(snapshot: RunSnapshot | undefined): QcIndexingOutcome {
  return qcIndexingOutcome(
    snapshot?.events ?? [],
    snapshot?.eventsTruncated ?? false,
  );
}

export function shouldPollRunSnapshot(
  snapshot: RunSnapshot | undefined,
  activeView: RunDetailView,
): boolean {
  const status = snapshot?.run?.status;
  if (!status) return false;
  if (pollingStatuses.has(status)) return true;
  if (status !== 'succeeded') return false;
  if (activeView === 'artifacts') {
    return extractionOutcome(snapshot).kind === 'pending';
  }
  return activeView === 'qc' && qcOutcome(snapshot).kind === 'pending';
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

function emptySnapshot(run: RunRecordResponse): RunSnapshot {
  return {
    run,
    events: [],
    eventsTruncated: false,
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
      eventsTruncated: false,
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
      eventsTruncated: false,
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
    eventsTruncated: events.truncated,
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
  workflowClient,
  executionAvailability,
  validationResult = null,
  runClient,
  runId = null,
  onRunCreated,
  beginPreflight = false,
  preflightRequestId = null,
  onPreflightConsumed,
  activeView = 'activity',
  selectedArtifactId = null,
  onViewChange = () => undefined,
  onArtifactSelect = () => undefined,
}: RunProgressPanelProps) {
  const [localRunId, setLocalRunId] = useState<string | null>(null);
  const [activeLogStream, setActiveLogStream] = useState<'stdout' | 'stderr'>('stdout');
  const [actionIssues, setActionIssues] = useState<Issue[]>([]);
  const [actionIssueKind, setActionIssueKind] = useState<ActionIssueKind>(null);
  const [cancellationRequested, setCancellationRequested] = useState(false);
  const [pollingPaused, setPollingPaused] = useState(false);
  const [pollRevision, setPollRevision] = useState(0);
  const preflightConsumedFor = useRef<string | null>(null);
  const pollStartedAt = useRef(Date.now());
  const queryClient = useQueryClient();
  const targetRunId = runId ?? localRunId;
  const createAuthority =
    workflowId !== null &&
    validationResult?.ok === true &&
    validationResult.snapshot !== null
      ? `${workflowId}\u0000${validationResult.snapshot.snapshot_id}`
      : null;
  const createAuthorityRef = useRef(createAuthority);
  createAuthorityRef.current = createAuthority;

  const runQuery = useQuery({
    queryKey: runProgressQueryKey(targetRunId),
    queryFn: () => loadRunSnapshot(runClient, targetRunId!),
    enabled: targetRunId !== null,
    refetchInterval: (query) => {
      if (!shouldPollRunSnapshot(query.state.data, activeView)) return false;
      if (pollingPaused) return false;
      return POLL_INTERVAL_MS;
    },
  });
  const snapshot = runQuery.data;
  const run = snapshot?.run ?? null;
  const availabilityWorkflowId = run?.workflow_id ?? workflowId;
  const workflowAvailabilityQuery = useQuery({
    queryKey: ['workflow', availabilityWorkflowId, 'execution-availability'],
    queryFn: () => workflowClient!.getWorkflow(availabilityWorkflowId!),
    enabled:
      executionAvailability === undefined &&
      workflowClient !== undefined &&
      availabilityWorkflowId !== null,
    retry: false,
  });
  const resolvedExecutionAvailability =
    executionAvailability !== undefined
      ? executionAvailability
      : workflowAvailabilityQuery.data?.workflow?.availability.execution ?? null;

  const clearActionIssues = () => {
    setActionIssues([]);
    setActionIssueKind(null);
  };

  const showActionIssues = (kind: Exclude<ActionIssueKind, null>, issues: Issue[]) => {
    setActionIssueKind(kind);
    setActionIssues(issues);
  };

  const resumePolling = () => {
    pollStartedAt.current = Date.now();
    setPollingPaused(false);
    setPollRevision((revision) => revision + 1);
  };

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
        showActionIssues('preflight', response.issues);
        return;
      }
      clearActionIssues();
      if (response.run) updateRun(response.run);
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
    onError: () => {
      showActionIssues('preflight', [
        safeUiIssue('PREFLIGHT_UNAVAILABLE', 'Preflight could not be started. Try again.'),
      ]);
    },
  });

  const startMutation = useMutation({
    mutationFn: (id: string) => runClient.startRun(id),
    onSuccess: (response) => {
      if (!response.ok) {
        showActionIssues('start', response.issues);
        void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
        return;
      }
      clearActionIssues();
      if (response.run) updateRun(response.run);
      resumePolling();
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
    onError: () => {
      showActionIssues('start', [
        safeUiIssue(
          'RUN_START_UNAVAILABLE',
          'Could not confirm whether the run was submitted. Refresh before retrying.',
        ),
      ]);
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
  });

  const cancelMutation = useMutation({
    mutationFn: (id: string) => runClient.cancelRun(id),
    onSuccess: (response) => {
      if (!response.ok) {
        showActionIssues('cancel', response.issues);
        void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
        return;
      }
      clearActionIssues();
      if (response.run) {
        updateRun(response.run);
        setCancellationRequested(response.run.status === 'running');
      }
      resumePolling();
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
    onError: () => {
      showActionIssues('cancel', [
        safeUiIssue(
          'RUN_CANCELLATION_UNAVAILABLE',
          'Could not confirm cancellation. Refresh or retry; the canonical run status remains authoritative.',
        ),
      ]);
      void queryClient.invalidateQueries({ queryKey: runProgressQueryKey(targetRunId) });
    },
  });

  const createMutation = useMutation({
    mutationFn: (attempt: CreateRunAttempt) =>
      runClient.createRun(attempt.workflowId, {
        snapshot_id: attempt.snapshotId,
      }),
  });

  useEffect(() => {
    setLocalRunId(null);
    setActiveLogStream('stdout');
    clearActionIssues();
    setCancellationRequested(false);
    setPollingPaused(false);
    preflightConsumedFor.current = null;
    pollStartedAt.current = Date.now();
    setPollRevision((revision) => revision + 1);
  }, [runId, runClient, workflowId, validationResult?.snapshot?.snapshot_id]);

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

  const artifactOutcome = extractionOutcome(snapshot);
  const qcIndexing = qcOutcome(snapshot);
  const shouldPollSnapshot = shouldPollRunSnapshot(snapshot, activeView);
  const cancellationEventPresent =
    snapshot?.events.some((event) => event.event_type === 'cancellation_requested') ?? false;
  const showCancellationRequested =
    run?.status === 'running' && (cancellationRequested || cancellationEventPresent);

  useEffect(() => {
    if (run && run.status !== 'running') setCancellationRequested(false);
  }, [run]);

  useEffect(() => {
    if (!shouldPollSnapshot) {
      setPollingPaused(false);
      return;
    }
    const remaining = MAX_POLL_DURATION_MS - (Date.now() - pollStartedAt.current);
    if (remaining <= 0) {
      setPollingPaused(true);
      return;
    }
    setPollingPaused(false);
    const timer = window.setTimeout(() => setPollingPaused(true), remaining);
    return () => window.clearTimeout(timer);
  }, [pollRevision, shouldPollSnapshot]);

  useEffect(() => {
    if (!run || actionIssueKind === null) return;
    const startWasConfirmed =
      actionIssueKind === 'start' && run.status !== 'planned';
    const cancellationReachedTerminal =
      actionIssueKind === 'cancel' && terminalStatuses.has(run.status);
    if (startWasConfirmed || cancellationReachedTerminal) clearActionIssues();
  }, [actionIssueKind, run?.status]);

  const canCreateRun =
    createAuthority !== null && resolvedExecutionAvailability === 'available';
  const executionActionPending = startMutation.isPending || cancelMutation.isPending;

  async function handleCreateRun() {
    clearActionIssues();
    const snapshotId = validationResult?.snapshot?.snapshot_id;
    if (
      !workflowId ||
      !snapshotId ||
      createAuthority === null ||
      resolvedExecutionAvailability !== 'available'
    ) {
      return;
    }
    const attempt: CreateRunAttempt = {
      authority: createAuthority,
      workflowId,
      snapshotId,
    };
    try {
      const response = await createMutation.mutateAsync(attempt);
      if (createAuthorityRef.current !== attempt.authority) {
        showActionIssues('create', [
          safeUiIssue(
            'RUN_CREATE_INPUTS_CHANGED',
            'Inputs changed while run creation was running. The earlier request may have created a run; review canonical runs before retrying.',
          ),
        ]);
        return;
      }
      if (!response.ok) {
        showActionIssues('create', response.issues);
        return;
      }
      if (!response.run) {
        showActionIssues('create', [
          safeUiIssue('RUN_RECORD_MISSING', 'The run record was missing.'),
        ]);
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
      if (createAuthorityRef.current !== attempt.authority) {
        showActionIssues('create', [
          safeUiIssue(
            'RUN_CREATE_INPUTS_CHANGED',
            'Inputs changed while run creation was running. The earlier request may have created a run; review canonical runs before retrying.',
          ),
        ]);
        return;
      }
      showActionIssues('create', [
        safeUiIssue('RUN_CREATE_UNAVAILABLE', 'The run record could not be created. Try again.'),
      ]);
    }
  }

  async function handleRefresh(): Promise<void> {
    clearActionIssues();
    resumePolling();
    const result = await runQuery.refetch();
    if (result.isError) {
      throw new Error('canonical run status refresh failed');
    }
  }

  const queryIssues = snapshot?.issues ?? [];
  const unexpectedQueryIssues = runQuery.error
    ? [safeUiIssue('RUN_PROGRESS_UNAVAILABLE', 'Run progress could not be loaded. Try again.')]
    : [];
  const visibleIssues = [...actionIssues, ...unexpectedQueryIssues, ...queryIssues];
  const canCancelRun = run !== null && cancellableStatuses.has(run.status);
  const showStartRun = run?.status === 'planned';
  const canStartRun =
    showStartRun && resolvedExecutionAvailability === 'available';
  const canPreflightRun = run?.status === 'created';
  const showNoRunPlaceholder = !run && !runQuery.isLoading && visibleIssues.length === 0 && !runId;
  const showMissingRunError =
    !run && !runQuery.isLoading && visibleIssues.length > 0 && runId !== null;

  function handleTabKey(
    event: KeyboardEvent<HTMLButtonElement>,
    view: RunDetailView,
  ) {
    const views: RunDetailView[] = ['activity', 'artifacts', 'qc'];
    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    const currentIndex = views.indexOf(view);
    const nextView =
      event.key === 'Home'
        ? views[0]
        : event.key === 'End'
          ? views.at(-1)!
          : event.key === 'ArrowRight'
            ? views[(currentIndex + 1) % views.length]
            : views[(currentIndex - 1 + views.length) % views.length];
    onViewChange(nextView);
    document.getElementById(`run-${nextView}-tab`)?.focus();
  }

  return (
    <div className="space-y-4" data-testid="run-progress-panel">
      {(createAuthority !== null || showStartRun) &&
        resolvedExecutionAvailability !== 'available' && (
          <ExecutionAvailabilityNotice
            availability={resolvedExecutionAvailability}
          />
        )}
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

          {pollingPaused && shouldPollSnapshot && (
            <div
              className="rounded border border-[var(--color-border)] bg-[var(--color-surface)] p-2 text-sm text-[var(--color-text-muted)]"
              role="status"
              data-testid="polling-paused"
            >
              Automatic updates paused. Refresh to resume.
            </div>
          )}

          {run.error && <RunIssuePanel issues={[run.error]} title="Run failure" />}

          <div className="flex flex-wrap gap-2">
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
            {showStartRun && (
              <Button
                variant="primary"
                onClick={() => {
                  if (canStartRun) startMutation.mutate(run.run_id);
                }}
                disabled={!canStartRun || executionActionPending}
                aria-label="Start run"
                data-testid="start-run-button"
              >
                <Play className="mr-1.5 h-4 w-4" aria-hidden="true" />
                {startMutation.isPending ? 'Starting run…' : 'Start run'}
              </Button>
            )}
            <Button
              variant="secondary"
              onClick={() => void handleRefresh().catch(() => undefined)}
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
          </div>

          {snapshot?.truncated && (
            <p className="text-xs text-[var(--color-text-muted)]" role="status">
              Showing the first {MAX_PAGES * PAGE_LIMIT} entries per stream. Additional entries are not loaded automatically.
            </p>
          )}

          {runId && (
            <div
              className="flex border-b border-[var(--color-border)]"
              role="tablist"
              aria-label="Run detail views"
            >
              {(['activity', 'artifacts', 'qc'] as const).map((view) => (
                <button
                  key={view}
                  id={`run-${view}-tab`}
                  type="button"
                  role="tab"
                  aria-selected={activeView === view}
                  aria-controls={`run-${view}-panel`}
                  tabIndex={activeView === view ? 0 : -1}
                  className={`border-b-2 px-3 py-2 text-sm font-medium focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)] focus:ring-inset ${
                    activeView === view
                      ? 'border-[var(--color-accent)] text-[var(--color-accent-hover)]'
                      : 'border-transparent text-[var(--color-text-muted)] hover:text-[var(--color-text)]'
                  }`}
                  onClick={() => onViewChange(view)}
                  onKeyDown={(event) => handleTabKey(event, view)}
                >
                  {view === 'activity'
                    ? 'Activity'
                    : view === 'artifacts'
                      ? 'Artifacts'
                      : 'QC'}
                </button>
              ))}
            </div>
          )}

          {(!runId || activeView === 'activity') && (
            <div
              id="run-activity-panel"
              role={runId ? 'tabpanel' : undefined}
              aria-labelledby={runId ? 'run-activity-tab' : undefined}
              className="space-y-3"
            >
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

          {runId && activeView === 'artifacts' && (
            <div
              id="run-artifacts-panel"
              role="tabpanel"
              aria-labelledby="run-artifacts-tab"
              className="min-w-0"
            >
              <ArtifactBrowser
                runId={run.run_id}
                runStatus={run.status}
                outcome={artifactOutcome}
                selectedArtifactId={selectedArtifactId}
                onSelectArtifact={onArtifactSelect}
                onRefreshStatus={handleRefresh}
              />
            </div>
          )}

          {runId && activeView === 'qc' && (
            <div
              id="run-qc-panel"
              role="tabpanel"
              aria-labelledby="run-qc-tab"
              className="min-w-0"
            >
              <QcWorkbench
                runId={run.run_id}
                runStatus={run.status}
                outcome={qcIndexing}
                onOpenSourceArtifact={onArtifactSelect}
                onRefreshStatus={handleRefresh}
              />
            </div>
          )}
        </div>
      )}

      {!runId && (
        <div className="flex flex-wrap gap-2">
          <Button
            variant="primary"
            onClick={handleCreateRun}
            disabled={!canCreateRun || createMutation.isPending}
            aria-label="Create run"
            data-testid="create-run-button"
          >
            {createMutation.isPending ? 'Creating run…' : 'Create run'}
          </Button>
        </div>
      )}
    </div>
  );
}
