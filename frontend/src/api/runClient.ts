import type { Issue } from './types';
import type {
  RunCreateRequest,
  RunRecordResponse,
  RunResponse,
  RunEventResponse,
  RunEventsResponse,
  RunLogsResponse,
} from './runTypes';

export interface RunApiError extends Error {
  status: number;
  statusText: string;
}

export interface RunApiClient {
  createRun(workflowId: string, request: RunCreateRequest): Promise<RunResponse>;
  preflightRun(runId: string): Promise<RunResponse>;
  getRun(runId: string): Promise<RunResponse>;
  listRunEvents(
    runId: string,
    options?: { after?: string; limit?: number },
  ): Promise<RunEventsResponse>;
  listRunLogs(
    runId: string,
    options?: { streamName?: string; after?: string; limit?: number },
  ): Promise<RunLogsResponse>;
  cancelRun(runId: string): Promise<RunResponse>;
}

function getBaseUrl(): string {
  const envBase = import.meta.env?.VITE_API_BASE_URL;
  return typeof envBase === 'string' ? envBase : '';
}

async function getJson<T>(url: string): Promise<T> {
  const response = await fetch(url, { method: 'GET', credentials: 'omit' });
  return parseResponse<T>(response);
}

async function postJson<T>(url: string, body: unknown): Promise<T> {
  const response = await fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
    credentials: 'omit',
  });
  return parseResponse<T>(response);
}

async function parseResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    const error = new Error(
      `Run API request failed: ${response.status} ${response.statusText}`,
    ) as RunApiError;
    error.status = response.status;
    error.statusText = response.statusText;
    throw error;
  }
  return response.json() as Promise<T>;
}

export function createRunApiClient(baseUrl?: string): RunApiClient {
  const root = baseUrl ?? getBaseUrl();

  return {
    async createRun(workflowId, request) {
      const url = `${root}/api/v1/workflows/${encodeURIComponent(workflowId)}/runs`;
      return postJson<RunResponse>(url, request);
    },

    async getRun(runId) {
      const url = `${root}/api/v1/runs/${encodeURIComponent(runId)}`;
      return getJson<RunResponse>(url);
    },

    async preflightRun(runId) {
      const url = `${root}/api/v1/runs/${encodeURIComponent(runId)}/preflight`;
      return postJson<RunResponse>(url, {});
    },

    async listRunEvents(runId, options = {}) {
      const params = new URLSearchParams();
      if (options.after !== undefined) params.set('after', options.after);
      if (options.limit !== undefined) params.set('limit', String(options.limit));
      const query = params.toString();
      const url = `${root}/api/v1/runs/${encodeURIComponent(runId)}/events${query ? `?${query}` : ''}`;
      return getJson<RunEventsResponse>(url);
    },

    async listRunLogs(runId, options = {}) {
      const params = new URLSearchParams();
      if (options.streamName !== undefined) params.set('stream_name', options.streamName);
      if (options.after !== undefined) params.set('after', options.after);
      if (options.limit !== undefined) params.set('limit', String(options.limit));
      const query = params.toString();
      const url = `${root}/api/v1/runs/${encodeURIComponent(runId)}/logs${query ? `?${query}` : ''}`;
      return getJson<RunLogsResponse>(url);
    },

    async cancelRun(runId) {
      const url = `${root}/api/v1/runs/${encodeURIComponent(runId)}/cancel`;
      return postJson<RunResponse>(url, {});
    },
  };
}

function runNotFoundIssue(runId: string): Issue {
  return {
    code: 'RUN_NOT_FOUND',
    message: 'Run was not found.',
    severity: 'error',
    path: 'run_id',
    source: 'stub',
    technical_message: null,
    hint: null,
    context: { run_id: runId },
  };
}

export function createStubRunApiClient(): RunApiClient {
  let counter = 0;
  const runs = new Map<string, RunRecordResponse>();
  const cancelledRuns = new Set<string>();
  const preflightedRuns = new Set<string>();

  function nextRunId(): string {
    counter += 1;
    return `stub-run-${counter}`;
  }

  function now(): string {
    return new Date().toISOString();
  }

  return {
    async createRun(workflowId, request) {
      const runId = nextRunId();
      const timestamp = now();
      const run: RunRecordResponse = {
        run_id: runId,
        workflow_id: workflowId,
        inputs: {
          config: request.config,
          samples: request.samples,
          options: request.options ?? {},
        },
        status: 'created',
        created_at: timestamp,
        updated_at: timestamp,
        started_at: null,
        ended_at: null,
        current_stage: null,
        cancellation_reason: null,
        error: null,
        tags: request.tags ?? {},
      };
      runs.set(runId, run);
      return { ok: true, run, issues: [] };
    },

    async getRun(runId) {
      const run = runs.get(runId);
      if (!run) {
        return { ok: false, run: null, issues: [runNotFoundIssue(runId)] };
      }
      return { ok: true, run, issues: [] };
    },

    async preflightRun(runId) {
      const run = runs.get(runId);
      if (!run) {
        return { ok: false, run: null, issues: [runNotFoundIssue(runId)] };
      }
      if (run.status !== 'created') {
        return {
          ok: false,
          run: null,
          issues: [
            {
              code: 'PREFLIGHT_ALREADY_TRIGGERED',
              message: 'Preflight has already been triggered for this run.',
              severity: 'error',
              path: 'run_id',
              source: 'stub',
              technical_message: null,
              hint: null,
              context: { current_status: run.status },
            },
          ],
        };
      }
      run.status = 'planned';
      run.current_stage = 'preflight';
      run.updated_at = now();
      preflightedRuns.add(runId);
      return { ok: true, run, issues: [] };
    },

    async listRunEvents(runId) {
      const run = runs.get(runId);
      if (!run) {
        return { ok: false, run_id: runId, events: [], next_cursor: null, issues: [runNotFoundIssue(runId)] };
      }

      const events: RunEventResponse[] = [
        {
          event_id: `${runId}-evt-1`,
          run_id: runId,
          sequence: 1,
          event_type: 'status_changed',
          timestamp: run.created_at,
          status: 'created',
          stage: null,
          message: 'Run created.',
          context: { previous_status: null, new_status: 'created' },
          issue: null,
        },
      ];

      if (cancelledRuns.has(runId)) {
        events.push({
          event_id: `${runId}-evt-2`,
          run_id: runId,
          sequence: 2,
          event_type: 'status_changed',
          timestamp: run.updated_at,
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

      if (preflightedRuns.has(runId)) {
        events.push({
          event_id: `${runId}-evt-preflight`,
          run_id: runId,
          sequence: events.length + 1,
          event_type: 'preflight_completed',
          timestamp: run.updated_at,
          status: 'planned',
          stage: 'preflight',
          message: 'Stub preflight completed.',
          context: {},
          issue: null,
        });
      }

      return { ok: true, run_id: runId, events, next_cursor: null, issues: [] };
    },

    async listRunLogs(runId, options = {}) {
      const run = runs.get(runId);
      if (!run) {
        return {
          ok: false,
          run_id: runId,
          stream_name: options.streamName ?? 'stdout',
          chunks: [],
          next_cursor: null,
          issues: [runNotFoundIssue(runId)],
        };
      }
      return {
        ok: true,
        run_id: runId,
        stream_name: options.streamName ?? 'stdout',
        chunks: [],
        next_cursor: null,
        issues: [],
      };
    },

    async cancelRun(runId) {
      const run = runs.get(runId);
      if (!run) {
        return { ok: false, run: null, issues: [runNotFoundIssue(runId)] };
      }

      const terminalStatuses = ['succeeded', 'failed', 'cancelled'];
      if (!terminalStatuses.includes(run.status)) {
        run.status = 'cancelled';
        run.cancellation_reason = 'User requested cancellation.';
        run.updated_at = now();
        cancelledRuns.add(runId);
      }

      return { ok: true, run, issues: [] };
    },
  };
}
