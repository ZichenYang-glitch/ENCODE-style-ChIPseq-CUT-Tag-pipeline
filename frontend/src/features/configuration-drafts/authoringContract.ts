import { ApiError } from '../../api/fetcher';
import type {
  AuthoringCapabilityResponse,
  AuthoringConnectionResponse,
  AuthoringDraftDocument,
  AuthoringDraftWriteRequest,
  AuthoringListReceiptsParams,
  AuthoringProposalResponse,
  AuthoringQuestionRequest,
  AuthoringSessionResponse,
  AuthoringSessionSummary,
  AuthoringTurnCreateRequest,
  AuthoringTurnResponse,
  AuthoringValidationResponse,
} from '../../api/generated/models';
import { compareJsonKeys, isJsonObject, isPlainObject } from '../input-workbench/jsonSafety';
import type {
  AssistantConnectionView,
  AuthoringCapabilityView,
  CapabilityField,
  CapabilityFieldType,
  DraftValidationStatus,
  DraftValidationView,
  ExplanationView,
  ProposalDecisionState,
  ProposalView,
  ScalarFieldValue,
  StructuredQuestion,
  SuggestionIntent,
} from './model';

export const REASON_COPY: Record<string, string> = {
  DRAFT_EXPIRED: 'This draft has expired.',
  PROPOSAL_EXPIRED: 'This suggestion has expired and was not adopted.',
  PROPOSAL_INVALIDATED: 'This suggestion is no longer valid and was not adopted.',
  STALE_DRAFT: 'The draft changed on the server. Nothing was overwritten.',
  CONNECTION_OWNED_BY_OTHER_LOGIN:
    'Another login owns the assistant connection. Connect again from this login.',
  BRIDGE_UNAVAILABLE:
    'The assistant is unavailable. You can keep editing and validating by hand.',
  ADMISSION_FAILED: 'The assistant is unavailable.',
  AUTHENTICATION_REQUIRED: 'Sign in again before continuing.',
  NOT_FOUND: 'That configuration draft is not available.',
  FORBIDDEN: 'That configuration draft is not available to this login.',
  INVALID_REQUEST: 'The request was not accepted.',
  REQUEST_TOO_LARGE: 'The request is too large.',
  DECISION_CONFLICT: 'This suggestion was already decided. The draft was not changed again.',
  TURN_REQUEST_CONFLICT:
    'This request identity was already used with different content. It was not sent again.',
  TURN_RESULT_CONFLICT: 'This request already has a different result.',
  CONNECTION_NOT_READY: 'The assistant is not ready.',
  CONNECTION_LOST: 'The assistant connection ended before the request finished.',
  CAPABILITY_UNAVAILABLE: 'This workflow has no suggestion fields. You can still edit by hand.',
  WORKFLOW_UNAVAILABLE: 'This workflow is not available.',
  MODEL_OUTPUT_REJECTED: 'The assistant response was not accepted. The draft was not changed.',
  AUTH_SESSION_ENDED: 'The login that owned the assistant connection has ended.',
  CLOSED_BY_OWNER: 'The assistant connection is closed.',
  SAMPLE_SELECTION_REQUIRED: 'Select at least one sample revision before validation.',
  SCOPE_DENIED: 'The assistant request is outside the approved actions.',
  CAPABILITY_CHANGED: 'The suggestion catalog changed. Review the current fields before requesting again.',
};

const GENERIC_REASON = 'The request was not accepted.';

export function reasonCopy(code: string | null | undefined): string {
  if (!code) return GENERIC_REASON;
  return REASON_COPY[code] ?? GENERIC_REASON;
}

export function stableJson(value: unknown): string {
  return JSON.stringify(sortJson(value));
}

function sortJson(value: unknown): unknown {
  if (Array.isArray(value)) return value.map(sortJson);
  if (!isPlainObject(value)) return value;
  const sorted: Record<string, unknown> = {};
  for (const key of Object.keys(value).sort(compareJsonKeys)) {
    sorted[key] = sortJson(value[key]);
  }
  return sorted;
}

export function receiptListParams(
  cursor: string | null,
  limit = 20,
): AuthoringListReceiptsParams {
  return typeof cursor === 'string' ? { cursor, limit } : { limit };
}

export function stableDocumentSnapshot(config: unknown, options: unknown): string {
  return stableJson({ config, options });
}

export interface CanonicalDraftContent {
  config?: unknown;
  options?: unknown;
  samples?: unknown;
  sample_revision_ids?: readonly string[] | null;
  reference_profile_revision_id?: string | null;
  input_selections?: ReadonlyArray<{
    input_file_revision_ids?: readonly string[] | null;
    input_use_key: string;
    occurrence: number;
  }> | null;
}

/**
 * Compare the draft document the server already normalizes.
 * Object key order is ignored. Array order is kept. This is not a success receipt.
 */
export function canonicalDraftSnapshot(content: CanonicalDraftContent): string {
  return stableJson({
    config: content.config ?? {},
    input_selections: [...(content.input_selections ?? [])].map((item) => ({
      input_file_revision_ids: [...(item.input_file_revision_ids ?? [])],
      input_use_key: item.input_use_key,
      occurrence: item.occurrence,
    })),
    options: content.options ?? {},
    reference_profile_revision_id: content.reference_profile_revision_id ?? null,
    sample_revision_ids: [...(content.sample_revision_ids ?? [])],
    samples: content.samples ?? null,
  });
}

export type DraftSaveDisposition = 'confirmed' | 'unconfirmed' | 'conflict';

export function draftSaveDisposition(input: {
  sent: CanonicalDraftContent;
  observed: CanonicalDraftContent | null | undefined;
  previousRevision: number;
  observedRevision: number;
  acknowledged: boolean;
}): DraftSaveDisposition {
  if (!input.observed) return 'unconfirmed';
  const sameDocument =
    canonicalDraftSnapshot(input.sent) === canonicalDraftSnapshot(input.observed);
  if (!sameDocument) {
    if (input.acknowledged || input.observedRevision !== input.previousRevision) {
      return 'conflict';
    }
    return 'unconfirmed';
  }
  if (!input.acknowledged && input.observedRevision === input.previousRevision) {
    return 'unconfirmed';
  }
  return 'confirmed';
}

export const MISSING_WORKFLOW_SAMPLE_ROWS =
  'This draft has no workflow sample rows. Checked project sample revisions are not FASTQ bindings, and this page cannot add file paths. A snapshot needs sample rows the server already holds.';

export function workflowSampleRowsAbsent(samples: unknown): boolean {
  if (samples == null) return true;
  if (Array.isArray(samples)) return samples.length === 0;
  return false;
}

export function withMissingWorkflowSampleContext(
  summary: string,
  issues: readonly { code: string }[],
  samples: unknown,
): string {
  if (!workflowSampleRowsAbsent(samples)) return summary;
  if (!issues.some((issue) => issue.code === 'ENCODE_REFERENCE_BINDING_CONFLICT')) {
    return summary;
  }
  if (summary.includes(MISSING_WORKFLOW_SAMPLE_ROWS)) return summary;
  return `${summary} ${MISSING_WORKFLOW_SAMPLE_ROWS}`.trim();
}

export interface ResponseScope {
  userId: string;
  projectId: string;
  workflowId: string;
  sid: string;
}

export function responseIsCurrent(requested: ResponseScope, current: ResponseScope): boolean {
  return (
    requested.userId === current.userId &&
    requested.projectId === current.projectId &&
    requested.workflowId === current.workflowId &&
    requested.sid === current.sid
  );
}

export interface TurnScope extends ResponseScope {
  schemaDigest: string;
  projectionDigest: string | null;
}

export function turnResponseIsCurrent(requested: TurnScope, current: TurnScope): boolean {
  return (
    responseIsCurrent(requested, current) &&
    requested.schemaDigest === current.schemaDigest &&
    requested.projectionDigest === current.projectionDigest
  );
}

export type AuthoringFailureKind =
  | 'unauthorized'
  | 'not-found'
  | 'expired'
  | 'stale'
  | 'capability-changed'
  | 'decision-conflict'
  | 'bridge-unavailable'
  | 'network'
  | 'rejected';

export interface AuthoringFailure {
  kind: AuthoringFailureKind;
  message: string;
  reasonCode: string | null;
}

export function classifyAuthoringError(error: unknown): AuthoringFailure {
  if (!(error instanceof ApiError)) {
    return {
      kind: 'network',
      message: 'The response was not confirmed. Nothing was marked complete.',
      reasonCode: null,
    };
  }
  const reasonCode = error.code === 'API_ERROR' ? null : error.code;
  if (error.status === 401 || reasonCode === 'AUTHENTICATION_REQUIRED') {
    return { kind: 'unauthorized', message: reasonCopy('AUTHENTICATION_REQUIRED'), reasonCode };
  }
  if (error.status === 404 || reasonCode === 'NOT_FOUND' || reasonCode === 'WORKFLOW_UNAVAILABLE') {
    return { kind: 'not-found', message: reasonCopy(reasonCode ?? 'NOT_FOUND'), reasonCode };
  }
  if (error.status === 410 || reasonCode === 'DRAFT_EXPIRED' || reasonCode === 'PROPOSAL_EXPIRED') {
    return {
      kind: 'expired',
      message: reasonCopy(reasonCode ?? 'DRAFT_EXPIRED'),
      reasonCode,
    };
  }
  if (reasonCode === 'BRIDGE_UNAVAILABLE' || reasonCode === 'ADMISSION_FAILED') {
    return { kind: 'bridge-unavailable', message: reasonCopy(reasonCode), reasonCode };
  }
  if (reasonCode === 'STALE_DRAFT') {
    return { kind: 'stale', message: reasonCopy(reasonCode), reasonCode };
  }
  if (reasonCode === 'CAPABILITY_CHANGED') {
    return { kind: 'capability-changed', message: reasonCopy(reasonCode), reasonCode };
  }
  if (reasonCode === 'DECISION_CONFLICT') {
    return { kind: 'decision-conflict', message: reasonCopy(reasonCode), reasonCode };
  }
  return { kind: 'rejected', message: reasonCopy(reasonCode), reasonCode };
}

export function connectionView(
  connection: AuthoringConnectionResponse,
  options: { draftRestored: boolean; connectPending: boolean },
): AssistantConnectionView {
  if (options.connectPending) {
    return { state: 'connecting', draftRestored: options.draftRestored };
  }
  if (connection.state === 'ready') {
    return { state: 'available', draftRestored: options.draftRestored };
  }
  if (
    connection.state === 'blocked' &&
    (connection.reason_code === 'BRIDGE_UNAVAILABLE' ||
      connection.reason_code === 'ADMISSION_FAILED')
  ) {
    return { state: 'bridge-unavailable', draftRestored: options.draftRestored };
  }
  if (connection.state === 'blocked') {
    return { state: 'disconnected', draftRestored: options.draftRestored };
  }
  if (
    connection.reason_code === 'CLOSED_BY_OWNER' ||
    connection.reason_code === 'AUTH_SESSION_ENDED' ||
    connection.reason_code === 'CONNECTION_LOST'
  ) {
    return { state: 'ended', draftRestored: options.draftRestored };
  }
  return { state: 'disconnected', draftRestored: options.draftRestored };
}

const FIELD_TYPES = new Set<CapabilityFieldType>([
  'integer',
  'number',
  'boolean',
  'string',
]);

export function capabilityView(
  capability: AuthoringCapabilityResponse | null | undefined,
): AuthoringCapabilityView | null {
  if (!capability) return null;
  const intents = capability.allowed_intents.filter(
    (intent): intent is SuggestionIntent =>
      intent === 'suggest' || intent === 'explain_issues',
  );
  const fields: CapabilityField[] = [];
  for (const field of capability.fields) {
    if (!FIELD_TYPES.has(field.type as CapabilityFieldType)) continue;
    const type = field.type as CapabilityFieldType;
    const minimum = numberOrNull(field.constraints.minimum);
    const maximum = numberOrNull(field.constraints.maximum);
    fields.push({
      token: field.token,
      label: field.label,
      type,
      ...(minimum !== undefined ? { minimum } : {}),
      ...(maximum !== undefined ? { maximum } : {}),
    });
  }
  return {
    available: capability.available,
    unavailableReason: capability.available
      ? null
      : reasonCopy(capability.reason_code ?? 'CAPABILITY_UNAVAILABLE'),
    intents,
    fields,
    issues: capability.allowed_issue_codes.map((issue) => ({
      code: issue.code,
      label: issue.label,
    })),
  };
}

function numberOrNull(value: unknown): number | null | undefined {
  if (value === undefined) return undefined;
  if (value === null) return null;
  return typeof value === 'number' && Number.isFinite(value) ? value : undefined;
}

export interface QuestionRequestResult {
  question: AuthoringQuestionRequest;
  omittedNonIntegerPreferences: boolean;
  issueCodeNotSent: boolean;
}

export function questionRequest(
  question: StructuredQuestion,
  capability: AuthoringCapabilityView,
): { ok: true; value: QuestionRequestResult } | { ok: false; message: string } {
  if (!capability.available || !capability.intents.includes(question.intent)) {
    return { ok: false, message: 'Choose a request type from the current list.' };
  }
  const fields = new Map(capability.fields.map((field) => [field.token, field]));
  const tokens = question.fieldTokens.filter((token) => fields.has(token));
  if (question.intent === 'suggest' && (tokens.length === 0 || tokens.length !== question.fieldTokens.length)) {
    return { ok: false, message: 'Choose a field from the current list.' };
  }
  const preferences: { field_token: string; value: number }[] = [];
  let omittedNonIntegerPreferences = false;
  for (const preference of question.preferences) {
    const field = fields.get(preference.fieldToken);
    if (
      field?.type === 'integer' &&
      tokens.includes(preference.fieldToken) &&
      typeof preference.value === 'number' &&
      Number.isSafeInteger(preference.value)
    ) {
      preferences.push({ field_token: preference.fieldToken, value: preference.value });
      continue;
    }
    omittedNonIntegerPreferences = true;
  }
  return {
    ok: true,
    value: {
      question: {
        intent: question.intent,
        field_tokens: question.intent === 'explain_issues' ? [] : tokens,
        preferences: question.intent === 'explain_issues' ? [] : preferences,
      },
      omittedNonIntegerPreferences:
        question.intent === 'suggest' && omittedNonIntegerPreferences,
      issueCodeNotSent: question.intent === 'explain_issues',
    },
  };
}

export function readProjectedScalar(value: unknown): ScalarFieldValue {
  if (!isPlainObject(value) || value.present === false) return { presence: 'missing' };
  if (value.present !== true) return { presence: 'missing' };
  if (value.value === null) return { presence: 'null' };
  if (
    typeof value.value === 'string' ||
    typeof value.value === 'number' ||
    typeof value.value === 'boolean'
  ) {
    return { presence: 'present', value: value.value };
  }
  return { presence: 'missing' };
}

export function proposalView(
  turn: AuthoringTurnResponse,
  capability: AuthoringCapabilityView | null,
): ProposalView | null {
  if (turn.result_kind === 'explanation' || turn.adoptable === false && turn.explanation) {
    return null;
  }
  if (!turn.proposal && turn.result_kind !== 'proposal') {
    if (turn.status === 'generating' || turn.status === 'failed') {
      return {
        state: turnState(turn),
        changes: [],
        digestAvailable: typeof turn.proposal_digest === 'string' && turn.proposal_digest.length > 0,
        adoptable: false,
      };
    }
    return null;
  }
  return {
    state: turnState(turn),
    changes: proposalChanges(turn.proposal, capability),
    digestAvailable: typeof turn.proposal_digest === 'string' && turn.proposal_digest.length > 0,
    adoptable: turn.adoptable === true && turnState(turn) === 'pending',
  };
}

function turnState(turn: AuthoringTurnResponse): ProposalDecisionState {
  if (turn.status === 'generating') return 'unconfirmed';
  if (turn.reason_code === 'PROPOSAL_EXPIRED') return 'expired';
  if (turn.reason_code === 'PROPOSAL_INVALIDATED' || turn.proposal_state === 'invalidated') {
    return 'invalidated';
  }
  if (turn.reason_code === 'STALE_DRAFT') return 'conflict';
  if (turn.proposal_state === 'adopted') return 'adopted';
  if (turn.proposal_state === 'rejected') return 'rejected';
  if (turn.status === 'failed') return 'unconfirmed';
  return 'pending';
}

function proposalChanges(
  proposal: AuthoringProposalResponse | null | undefined,
  capability: AuthoringCapabilityView | null,
): ProposalView['changes'] {
  if (!proposal) return [];
  const fields = new Map((capability?.fields ?? []).map((field) => [field.token, field.label]));
  return proposal.changes.map((change) => {
    const label = fields.get(change.field_token);
    return {
      label: label ?? 'Unlabeled field',
      before: readProjectedScalar(change.before),
      after: readProjectedScalar(change.after),
      recognized: label !== undefined,
    };
  });
}

export function explanationView(
  turn: AuthoringTurnResponse,
  capability: AuthoringCapabilityView | null,
): ExplanationView | null {
  const codes = turn.explanation?.reason_codes;
  if (!codes || turn.result_kind !== 'explanation') return null;
  const labels = new Map(
    (capability?.issues ?? []).map((issue) => [issue.code, issue.label]),
  );
  return {
    contentAvailable: codes.length > 0,
    labels: codes.map((code) => approvedExplanationLabel(labels.get(code))),
  };
}

function approvedExplanationLabel(label: string | undefined): string {
  if (label && label.trim() && !/[/\\]|\n|traceback/i.test(label)) return label;
  return 'This explanation is not in the current catalog.';
}

export function validationView(input: {
  stored: StoredValidation | null;
  currentRevision: number | null;
  schemaDigest: string | null;
  dirty: boolean;
}): DraftValidationView {
  if (!input.stored) return { status: 'idle' };
  if (
    input.dirty ||
    input.currentRevision === null ||
    input.stored.revision !== input.currentRevision ||
    input.stored.schemaDigest !== input.schemaDigest
  ) {
    return { status: 'stale' };
  }
  return input.stored.view;
}

export interface StoredValidation {
  revision: number;
  schemaDigest: string;
  view: DraftValidationView;
}

export function validationFromResponse(
  response: AuthoringValidationResponse,
  schemaDigest: string,
): StoredValidation {
  return {
    revision: response.draft_revision,
    schemaDigest,
    view: validationSuccessView(response),
  };
}

export function validationFromIssues(
  revision: number,
  schemaDigest: string,
  issues: readonly { code: string; message: string }[],
): StoredValidation {
  return {
    revision,
    schemaDigest,
    view: {
      status: 'issues',
      summary: issueSummary(issues),
    },
  };
}

function issueSummary(issues: readonly { code: string; message: string }[]): string {
  const summary = issues.map((issue) => publicIssueText(issue)).join(' ');
  return summary || reasonCopy(null);
}

function publicIssueText(issue: { code: string; message: string }): string {
  if (issue.code in REASON_COPY) return reasonCopy(issue.code);
  const message = issue.message.trim();
  if (message && !/[/\\]|\n|traceback/i.test(message)) return message;
  return reasonCopy(null);
}

function validationSuccessView(response: AuthoringValidationResponse): DraftValidationView {
  if (!response.ok) {
    return { status: 'issues', summary: issueSummary(response.issues ?? []) };
  }
  if (response.snapshot) {
    return {
      status: 'snapshot',
      snapshotId: response.snapshot.snapshot_id,
      revisionLabel: `Revision ${response.draft_revision}`,
      expiresLabel: formatServerTimestamp(response.snapshot.expires_at),
    };
  }
  return { status: 'no-snapshot', revisionLabel: `Revision ${response.draft_revision}` };
}

export function formatServerTimestamp(value: string | null | undefined): string {
  if (!value) return 'Not reported';
  const parsed = Date.parse(value);
  if (Number.isNaN(parsed)) return 'Not reported';
  return `${new Intl.DateTimeFormat('en', {
    dateStyle: 'medium',
    timeStyle: 'short',
    timeZone: 'UTC',
  }).format(new Date(parsed))} UTC`;
}

export function draftWriteRequest(input: {
  expectedRevision: number;
  config: unknown;
  options: unknown;
  sampleRevisionIds: readonly string[];
  referenceProfileRevisionId: string | null;
  preservedSamples: AuthoringDraftDocument['samples'];
  preservedInputSelections: AuthoringDraftDocument['input_selections'];
}): AuthoringDraftWriteRequest | null {
  if (!isJsonObject(input.config) || !isJsonObject(input.options)) return null;
  return {
    expected_revision: input.expectedRevision,
    config: input.config,
    options: input.options,
    sample_revision_ids: [...input.sampleRevisionIds],
    input_selections: input.preservedInputSelections.map((selection) => ({
      input_use_key: selection.input_use_key,
      occurrence: selection.occurrence,
      input_file_revision_ids: [...selection.input_file_revision_ids],
    })),
    ...(input.referenceProfileRevisionId
      ? { reference_profile_revision_id: input.referenceProfileRevisionId }
      : {}),
    ...(input.preservedSamples !== undefined ? { samples: input.preservedSamples } : {}),
  };
}

export function serverDocumentSnapshot(draft: AuthoringDraftDocument | null | undefined): string | null {
  if (!draft) return null;
  return stableDocumentSnapshot(draft.config, draft.options);
}

export function sameMembers(left: readonly string[], right: readonly string[]): boolean {
  if (left.length !== right.length) return false;
  const pending = [...right];
  for (const value of left) {
    const index = pending.indexOf(value);
    if (index < 0) return false;
    pending.splice(index, 1);
  }
  return pending.length === 0;
}

export function sessionScope(
  session: Pick<AuthoringSessionResponse, 'project_id' | 'workflow_id' | 'sid'>,
  userId: string,
): ResponseScope {
  return {
    userId,
    projectId: session.project_id,
    workflowId: session.workflow_id,
    sid: session.sid,
  };
}

export function turnScope(session: AuthoringSessionResponse, userId: string): TurnScope {
  return {
    ...sessionScope(session, userId),
    schemaDigest: session.schema_digest,
    projectionDigest: session.projection_digest ?? null,
  };
}

export function listLabel(summary: AuthoringSessionSummary, projectName: string, workflowName: string) {
  return {
    sid: summary.sid,
    projectName,
    workflowName,
    state: summary.draft_state,
    updatedLabel: `revision ${summary.draft_revision}`,
    expiresLabel: formatServerTimestamp(summary.expires_at),
  };
}

export function turnCreateBody(
  clientTurnId: string,
  expectedRevision: number,
  question: AuthoringQuestionRequest,
): AuthoringTurnCreateRequest {
  return {
    client_turn_id: clientTurnId,
    expected_revision: expectedRevision,
    question,
  };
}

export function validationStatusAllowsCurrent(status: DraftValidationStatus): boolean {
  return status === 'issues' || status === 'no-snapshot' || status === 'snapshot';
}

const VALIDATION_TRANSPORT_REASONS = new Set([
  'AUTHENTICATION_REQUIRED',
  'CAPABILITY_CHANGED',
  'DRAFT_EXPIRED',
  'FORBIDDEN',
  'INVALID_REQUEST',
  'NOT_FOUND',
  'REQUEST_TOO_LARGE',
  'STALE_DRAFT',
  'WORKFLOW_UNAVAILABLE',
]);

export function isValidationBusinessFailure(failure: AuthoringFailure): boolean {
  return (
    failure.kind === 'rejected' &&
    (failure.reasonCode === null || !VALIDATION_TRANSPORT_REASONS.has(failure.reasonCode))
  );
}
