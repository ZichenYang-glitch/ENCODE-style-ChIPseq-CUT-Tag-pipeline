export type ScalarPresence = 'missing' | 'null' | 'present';

export interface ScalarFieldValue {
  presence: ScalarPresence;
  /** Set only when presence is `present`. Empty text, 0, and false are real values. */
  value?: string | number | boolean;
}

export type DraftRecordState = 'active' | 'expired';

export interface ConfigurationDraftSummary {
  sid: string;
  projectName: string;
  workflowName: string;
  state: DraftRecordState;
  updatedLabel: string;
  expiresLabel: string;
}

export type DraftListStatus =
  | 'loading'
  | 'ready'
  | 'empty'
  | 'error'
  | 'unavailable';

export interface DraftListModel {
  status: DraftListStatus;
  drafts: readonly ConfigurationDraftSummary[];
  pageLimit: number;
  hasPrevious: boolean;
  hasNext: boolean;
  errorMessage?: string;
}

export type AssistantConnectionState =
  | 'disconnected'
  | 'connecting'
  | 'available'
  | 'ended'
  | 'bridge-unavailable';

export interface AssistantConnectionView {
  state: AssistantConnectionState;
  draftRestored: boolean;
}

export type SuggestionIntent = 'suggest' | 'explain_issues';

export type CapabilityFieldType = 'integer' | 'number' | 'boolean' | 'string';

export interface CapabilityField {
  token: string;
  label: string;
  type: CapabilityFieldType;
  minimum?: number | null;
  maximum?: number | null;
}

export interface CapabilityIssue {
  code: string;
  label: string;
}

export interface AuthoringCapabilityView {
  available: boolean;
  unavailableReason: string | null;
  intents: readonly SuggestionIntent[];
  fields: readonly CapabilityField[];
  issues: readonly CapabilityIssue[];
}

export interface ProposalChangeView {
  label: string;
  before: ScalarFieldValue;
  after: ScalarFieldValue;
  recognized: boolean;
}

export type ProposalDecisionState =
  | 'pending'
  | 'adopting'
  | 'adopted'
  | 'rejected'
  | 'expired'
  | 'conflict'
  | 'invalidated'
  | 'unconfirmed';

export interface ProposalView {
  state: ProposalDecisionState;
  changes: readonly ProposalChangeView[];
  digestAvailable: boolean;
  adoptable: boolean;
}

export interface ExplanationView {
  contentAvailable: boolean;
  labels: readonly string[];
}

export type DraftValidationStatus =
  | 'unavailable'
  | 'idle'
  | 'unconfirmed'
  | 'request-failed'
  | 'stale'
  | 'issues'
  | 'no-snapshot'
  | 'snapshot';

export interface DraftValidationView {
  status: DraftValidationStatus;
  summary?: string;
  snapshotId?: string;
  revisionLabel?: string;
  expiresLabel?: string;
}

export interface StructuredQuestion {
  intent: SuggestionIntent;
  fieldTokens: readonly string[];
  preferences: readonly {
    fieldToken: string;
    value: string | number | boolean;
  }[];
  issueCode: string | null;
}

export type SuggestionRequestStatus = 'idle' | 'unconfirmed' | 'failed';

export const DEFAULT_DRAFT_PAGE_LIMIT = 20;
export const MAX_DRAFT_PAGE_LIMIT = 50;

export const AUTHORING_WRITES_UNAVAILABLE =
  'This action waits for the reviewed configuration-draft interface. Nothing was saved, adopted, or validated.';

export function boundedPageLimit(value: number): number {
  if (!Number.isInteger(value) || value < 1) return DEFAULT_DRAFT_PAGE_LIMIT;
  return Math.min(value, MAX_DRAFT_PAGE_LIMIT);
}
