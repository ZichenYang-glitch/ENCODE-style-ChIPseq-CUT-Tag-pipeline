import type {
  RunCreateRequest as GeneratedRunCreateRequest,
} from './generated/models';
import type { Issue } from './types';

export interface RunRecordResponse {
  run_id: string;
  workflow_id: string;
  inputs: Record<string, unknown>;
  status: string;
  created_at: string;
  updated_at: string;
  started_at: string | null;
  ended_at: string | null;
  current_stage: string | null;
  cancellation_reason: string | null;
  error: Issue | null;
  tags: Record<string, string>;
}

export interface RunEventResponse {
  event_id: string;
  run_id: string;
  sequence: number;
  event_type: string;
  timestamp: string;
  status: string | null;
  stage: string | null;
  message: string;
  context: Record<string, unknown>;
  issue: Issue | null;
}

export interface RunLogChunkResponse {
  chunk_id: string;
  run_id: string;
  stream_name: string;
  sequence: number;
  timestamp: string;
  lines: string[];
}

export type RunCreateRequest = GeneratedRunCreateRequest;

export interface RunResponse {
  ok: boolean;
  run: RunRecordResponse | null;
  issues: Issue[];
}

export interface RunEventsResponse {
  ok: boolean;
  run_id: string;
  events: RunEventResponse[];
  next_cursor: string | null;
  issues: Issue[];
}

export interface RunLogsResponse {
  ok: boolean;
  run_id: string;
  stream_name: string;
  chunks: RunLogChunkResponse[];
  next_cursor: string | null;
  issues: Issue[];
}
