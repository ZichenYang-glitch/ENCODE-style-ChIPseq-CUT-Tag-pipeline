import type {
  ReferenceProfileRevisionResponse,
  RunCreateRequest as GeneratedRunCreateRequest,
} from './generated/models';
import type { Issue } from './types';

export type ReferenceProfileSummary = ReferenceProfileRevisionResponse;

export function readReferenceProfileSummary(
  value: unknown,
): ReferenceProfileSummary | null {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) {
    return null;
  }
  const item = value as Record<string, unknown>;
  if (
    typeof item.profile_id !== 'string' ||
    !/^refp_[0-9a-f]{32}$/.test(item.profile_id) ||
    typeof item.revision_id !== 'string' ||
    !/^refpr_[0-9a-f]{32}$/.test(item.revision_id) ||
    typeof item.revision_number !== 'number' ||
    !Number.isSafeInteger(item.revision_number) ||
    item.revision_number < 1 ||
    typeof item.display_name !== 'string' ||
    item.display_name.length < 1 ||
    item.display_name.length > 255 ||
    typeof item.organism !== 'string' ||
    item.organism.length < 1 ||
    item.organism.length > 255 ||
    typeof item.assembly !== 'string' ||
    item.assembly.length < 1 ||
    item.assembly.length > 255 ||
    typeof item.identity_sha256 !== 'string' ||
    !/^[0-9a-f]{64}$/.test(item.identity_sha256)
  ) {
    return null;
  }
  return {
    profile_id: item.profile_id,
    revision_id: item.revision_id,
    revision_number: item.revision_number,
    display_name: item.display_name,
    organism: item.organism,
    assembly: item.assembly,
    identity_sha256: item.identity_sha256,
  };
}

export interface RunRecordResponse {
  run_id: string;
  workflow_id: string;
  status: string;
  created_at: string;
  updated_at: string;
  started_at: string | null;
  ended_at: string | null;
  current_stage: string | null;
  cancellation_reason: string | null;
  error: Issue | null;
  tags: Record<string, string>;
  reference_profile?: ReferenceProfileSummary | null;
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
