import type {
  RunHistoryResponse,
  RunStatus,
  RunSummaryResponse,
} from '../../api/generated/models';
import { RunStatus as RunStatusValues } from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';

export const RUN_HISTORY_PAGE_SIZE = 50;
export const RUN_HISTORY_STATUSES = Object.values(RunStatusValues);

const SAFE_WORKFLOW_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,254}$/;
const SAFE_RUN_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$/;
const SAFE_STAGE = /^[A-Za-z][A-Za-z0-9_.-]{0,254}$/;
const SAFE_CURSOR = /^runhist_[A-Za-z0-9_-]+$/;
const ISO_TIMESTAMP = /^\d{4}-\d{2}-\d{2}T.*(?:Z|[+-]\d{2}:\d{2})$/;

export interface RunHistoryFilters {
  workflowId: string | null;
  status: RunStatus | null;
}

export type RunHistoryFilterResult =
  | { ok: true; filters: RunHistoryFilters }
  | { ok: false; filters: null };

export class RunHistoryProtocolError extends Error {
  constructor() {
    super('run history protocol response is invalid');
    this.name = 'RunHistoryProtocolError';
  }
}

export function parseRunHistoryFilters(
  searchParams: URLSearchParams,
): RunHistoryFilterResult {
  const workflows = searchParams.getAll('workflow_id');
  const statuses = searchParams.getAll('status');
  if (workflows.length > 1 || statuses.length > 1) {
    return { ok: false, filters: null };
  }
  const workflowId = workflows[0] ?? null;
  if (
    workflowId !== null &&
    (workflowId.trim() !== workflowId || !SAFE_WORKFLOW_ID.test(workflowId))
  ) {
    return { ok: false, filters: null };
  }
  const status = statuses[0] ?? null;
  if (
    status !== null &&
    !RUN_HISTORY_STATUSES.includes(status as RunStatus)
  ) {
    return { ok: false, filters: null };
  }
  return {
    ok: true,
    filters: { workflowId, status: status as RunStatus | null },
  };
}

export function runHistoryQueryKey(filters: RunHistoryFilters) {
  return ['run-history', filters.workflowId, filters.status] as const;
}

export function validateRunHistoryResponse(
  response: RunHistoryResponse,
  requestedCursor: string | undefined,
): RunHistoryResponse {
  if (
    response.ok !== true ||
    !Array.isArray(response.runs) ||
    !Array.isArray(response.issues) ||
    response.issues.length !== 0 ||
    !isSafeCursor(response.next_cursor) ||
    (requestedCursor !== undefined && response.next_cursor === requestedCursor)
  ) {
    throw new RunHistoryProtocolError();
  }
  const runs = response.runs.map(projectRunSummary);
  if (new Set(runs.map((run) => run.run_id)).size !== runs.length) {
    throw new RunHistoryProtocolError();
  }
  return { ok: true, runs, next_cursor: response.next_cursor, issues: [] };
}

export function flattenRunHistoryPages(
  pages: RunHistoryResponse[],
): RunSummaryResponse[] {
  const seen = new Set<string>();
  const runs: RunSummaryResponse[] = [];
  for (const page of pages) {
    for (const run of page.runs) {
      if (!seen.has(run.run_id)) {
        seen.add(run.run_id);
        runs.push(run);
      }
    }
  }
  return runs;
}

export function runHistoryPaginationAnomaly(
  pages: RunHistoryResponse[],
  pageParams: unknown[],
): boolean {
  const usedCursors = new Set<string>();
  const runIds = new Set<string>();
  for (const pageParam of pageParams) {
    if (typeof pageParam !== 'string') continue;
    if (usedCursors.has(pageParam)) return true;
    usedCursors.add(pageParam);
  }
  for (const [index, page] of pages.entries()) {
    for (const run of page.runs) {
      if (runIds.has(run.run_id)) return true;
      runIds.add(run.run_id);
    }
    if (page.next_cursor !== null) {
      const currentParam = pageParams[index];
      const followingParam = pageParams[index + 1];
      if (page.next_cursor === currentParam) return true;
      if (followingParam !== undefined) {
        if (page.next_cursor !== followingParam) return true;
      } else if (usedCursors.has(page.next_cursor)) {
        return true;
      }
    }
  }
  return false;
}

export function safeRunHistoryNextCursor(
  page: RunHistoryResponse,
  pageParams: unknown[],
): string | undefined {
  const cursor = page.next_cursor;
  if (cursor === null || !isSafeCursor(cursor)) return undefined;
  return pageParams.includes(cursor) ? undefined : cursor;
}

export function isRunHistoryCursorError(error: unknown): boolean {
  return (
    error instanceof RunHistoryProtocolError ||
    (error instanceof ApiError &&
      (error.code === 'RUN_HISTORY_CURSOR_INVALID' ||
        error.code === 'RUN_HISTORY_CURSOR_NOT_FOUND'))
  );
}

function projectRunSummary(value: unknown): RunSummaryResponse {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) {
    throw new RunHistoryProtocolError();
  }
  const run = value as Record<string, unknown>;
  if (
    typeof run.run_id !== 'string' ||
    !SAFE_RUN_ID.test(run.run_id) ||
    typeof run.workflow_id !== 'string' ||
    !SAFE_WORKFLOW_ID.test(run.workflow_id) ||
    typeof run.status !== 'string' ||
    !RUN_HISTORY_STATUSES.includes(run.status as RunStatus) ||
    !isTimestamp(run.created_at) ||
    !isTimestamp(run.updated_at) ||
    !isOptionalTimestamp(run.started_at) ||
    !isOptionalTimestamp(run.ended_at) ||
    !isOptionalStage(run.current_stage)
  ) {
    throw new RunHistoryProtocolError();
  }
  const created = Date.parse(run.created_at);
  const updated = Date.parse(run.updated_at);
  const started = run.started_at === null ? null : Date.parse(run.started_at);
  const ended = run.ended_at === null ? null : Date.parse(run.ended_at);
  if (
    updated < created ||
    (started !== null && (started < created || started > updated)) ||
    (ended !== null && (ended < created || ended > updated)) ||
    (started !== null && ended !== null && ended < started)
  ) {
    throw new RunHistoryProtocolError();
  }
  return {
    run_id: run.run_id,
    workflow_id: run.workflow_id,
    status: run.status as RunStatus,
    created_at: run.created_at,
    updated_at: run.updated_at,
    started_at: run.started_at,
    ended_at: run.ended_at,
    current_stage: run.current_stage,
  };
}

function isSafeCursor(value: unknown): value is string | null {
  return (
    value === null ||
    (typeof value === 'string' &&
      value.length <= 1024 &&
      SAFE_CURSOR.test(value))
  );
}

function isTimestamp(value: unknown): value is string {
  return (
    typeof value === 'string' &&
    value.length <= 64 &&
    ISO_TIMESTAMP.test(value) &&
    !Number.isNaN(Date.parse(value))
  );
}

function isOptionalTimestamp(value: unknown): value is string | null {
  return value === null || isTimestamp(value);
}

function isOptionalStage(value: unknown): value is string | null {
  return value === null || (typeof value === 'string' && SAFE_STAGE.test(value));
}
