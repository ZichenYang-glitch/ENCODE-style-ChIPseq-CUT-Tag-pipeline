import { useEffect, useMemo, useRef, useState } from 'react';
import type { Dispatch, SetStateAction } from 'react';
import { useMutation } from '@tanstack/react-query';
import { Play, ShieldCheck } from 'lucide-react';
import { useNavigate } from 'react-router-dom';
import { createRun } from '../../api/generated/runs/runs';
import { validateWorkflow } from '../../api/generated/workflows/workflows';
import type {
  InputFileRevisionSelectionRequest,
  IssueResponse,
  ValidationRequest,
  ValidatedInputSnapshotResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { readReferenceProfileSummary } from '../../api/runTypes';
import type { WorkflowAvailability } from '../../api/types';
import { Button } from '../../components/Button';
import { readExperimentBinding } from '../experiment-binding/handoff';
import { ExecutionAvailabilityNotice } from '../workflow-detail/WorkflowAvailability';
import type { InputDraftController } from './useInputDraft';

import type { SafeIssue, ValidationFeedback } from './validationFeedback';

interface SnapshotState {
  snapshot: ValidatedInputSnapshotResponse;
  revision: number;
}

interface ValidationAttempt {
  requestId: number;
  payload: ValidationRequest & { reference_profile_revision_id: string };
  revision: number;
  referenceProfileRevisionId: string;
  contextKey: string;
}

interface CreateAttempt {
  snapshotId: string;
  revision: number;
}

export interface ExperimentValidationRequest {
  projectId: string;
  sampleRevisionIds: readonly string[];
  inputSelections: readonly InputFileRevisionSelectionRequest[];
  experimentRevisionId: string | null;
}

type BindingOutcome =
  | { kind: 'idle' }
  | { kind: 'unresolved' }
  | { kind: 'frozen'; experimentRevisionId: string }
  | { kind: 'failed' }
  | { kind: 'unconfirmed' }
  | { kind: 'missing' };

interface ValidatedSubmissionProps {
  workflowId: string;
  draft: InputDraftController;
  availability: WorkflowAvailability | null;
  referenceSelectionAvailable: boolean;
  experimentRequest?: ExperimentValidationRequest | null;
  onValidationFeedback: Dispatch<SetStateAction<ValidationFeedback | null>>;
}

function safeIssues(issues: IssueResponse[] | undefined): SafeIssue[] {
  return (issues ?? []).map((issue) => ({
    code: issue.code,
    message: issue.message,
    severity:
      issue.severity === 'warning' || issue.severity === 'info'
        ? issue.severity
        : 'error',
    path: issue.path,
    hint: issue.hint,
  }));
}

function requestIssues(error: unknown, fallback: SafeIssue): SafeIssue[] {
  if (error instanceof ApiError && error.issues.length > 0) {
    return error.issues.map((issue) => ({
      code: issue.code,
      message: issue.message,
      severity:
        issue.severity === 'warning' || issue.severity === 'info'
          ? issue.severity
          : 'error',
      path: issue.path,
      hint: issue.hint,
    }));
  }
  return [{ ...fallback, severity: fallback.severity ?? 'error' }];
}

export function ValidatedSubmission({
  workflowId,
  draft,
  availability,
  referenceSelectionAvailable,
  experimentRequest = null,
  onValidationFeedback,
}: ValidatedSubmissionProps) {
  const navigate = useNavigate();
  const [snapshotState, setSnapshotState] = useState<SnapshotState | null>(null);
  const [bindingOutcome, setBindingOutcome] = useState<BindingOutcome>({ kind: 'idle' });
  // Creation outcomes survive editing and revalidation, including uncertain outcomes.
  const [issues, setIssues] = useState<SafeIssue[]>([]);
  const latestValidationRequest = useRef(0);
  const activeSnapshot =
    snapshotState?.revision === draft.state.semanticRevision
      ? snapshotState.snapshot
      : null;
  const executionAvailable = availability?.execution === 'available';
  const contextKey = experimentRequest === null
    ? ''
    : JSON.stringify({
        projectId: experimentRequest.projectId,
        sampleRevisionIds: experimentRequest.sampleRevisionIds,
        inputSelections: experimentRequest.inputSelections,
        experimentRevisionId: experimentRequest.experimentRevisionId,
      });
  const contextKeyRef = useRef(contextKey);
  contextKeyRef.current = contextKey;
  const previousContextKey = useRef(contextKey);
  useEffect(() => {
    const contextChanged = previousContextKey.current !== contextKey;
    previousContextKey.current = contextKey;
    setSnapshotState(null);
    setBindingOutcome({ kind: 'idle' });
    if (!contextChanged) return;
    onValidationFeedback((current) =>
      current === null
        ? current
        : { revision: -1, issues: [], changedInFlight: false },
    );
  }, [contextKey, draft.state.semanticRevision, onValidationFeedback]);

  const validationMutation = useMutation({
    mutationFn: (attempt: ValidationAttempt) =>
      validateWorkflow(workflowId, attempt.payload),
    onSuccess: (response, attempt) => {
      if (
        attempt.requestId !== latestValidationRequest.current ||
        attempt.contextKey !== contextKeyRef.current
      ) {
        return;
      }
      const report = (issues: SafeIssue[], notice?: string, snapshotId?: string) =>
        onValidationFeedback({ revision: attempt.revision, issues, notice, snapshotId });
      const failBinding = () => {
        if (experimentRequest !== null) setBindingOutcome({ kind: 'failed' });
      };
      if (attempt.revision !== draft.state.semanticRevision) {
        setSnapshotState(null);
        failBinding();
        onValidationFeedback({ revision: attempt.revision, issues: [], changedInFlight: true });
        return;
      }
      const responseIssues = safeIssues(response.issues);
      if (response.snapshot !== null) {
        const frozenReference = readReferenceProfileSummary(response.snapshot.reference_profile);
        if (frozenReference?.revision_id !== attempt.referenceProfileRevisionId) {
          setSnapshotState(null);
          failBinding();
          report([{
            code: 'REFERENCE_PROFILE_NOT_CONFIRMED',
            message: 'Backend validation did not confirm the selected reference revision.',
            path: 'reference_profile_revision_id',
          }]);
          return;
        }
      }
      if (response.ok && response.snapshot === null && !executionAvailable) {
        setSnapshotState(null);
        failBinding();
        report(responseIssues, 'Backend validation succeeded. No runnable snapshot was issued because execution is unavailable.');
        return;
      }
      if (!response.ok || response.snapshot === null || response.snapshot.workflow_id !== workflowId) {
        setSnapshotState(null);
        failBinding();
        report(responseIssues.length > 0 ? responseIssues : [{
          code: 'VALIDATION_NOT_CONFIRMED',
          message: 'Backend validation did not return a usable snapshot.',
        }]);
        return;
      }
      setSnapshotState({ snapshot: response.snapshot, revision: attempt.revision });
      if (experimentRequest !== null) {
        const binding = readExperimentBinding(response.snapshot.experiment_binding);
        setBindingOutcome(
          binding === null
            ? { kind: 'missing' }
            : binding.status === 'unresolved'
              ? { kind: 'unresolved' }
              : { kind: 'frozen', experimentRevisionId: binding.experimentRevisionId },
        );
      }
      report(responseIssues, executionAvailable
        ? 'Backend validation succeeded. This exact draft can create one run.'
        : 'Backend validation succeeded. This exact draft is saved, but execution is unavailable.',
        response.snapshot.snapshot_id);
    },
    onError: (error, attempt) => {
      if (
        attempt.requestId !== latestValidationRequest.current ||
        attempt.contextKey !== contextKeyRef.current
      ) {
        return;
      }
      setSnapshotState(null);
      if (experimentRequest !== null) {
        setBindingOutcome(
          error instanceof ApiError ? { kind: 'failed' } : { kind: 'unconfirmed' },
        );
      }
      if (attempt.revision !== draft.state.semanticRevision) {
        onValidationFeedback({ revision: attempt.revision, issues: [], changedInFlight: true });
        return;
      }
      onValidationFeedback({
        revision: attempt.revision,
        issues: requestIssues(error, {
          code: 'VALIDATION_UNAVAILABLE',
          message: 'Validation could not be confirmed. Retry when the API is available.',
        }),
      });
    },
  });

  const createMutation = useMutation({
    mutationFn: (attempt: CreateAttempt) =>
      createRun(workflowId, {
        snapshot_id: attempt.snapshotId,
      }),
    onSuccess: (response, attempt) => {
      if (attempt.revision !== draft.state.semanticRevision) {
        setSnapshotState(null);
        setIssues([
          {
            code: 'RUN_CREATE_INPUTS_CHANGED',
            message:
              'Inputs changed while run creation was running. The earlier request may have created a run; review canonical runs before retrying.',
          },
        ]);
        return;
      }
      const run = response.run ?? null;
      if (!response.ok || run === null) {
        setIssues(
          safeIssues(response.issues).concat(
            run === null
              ? [
                  {
                    code: 'RUN_CREATE_UNCONFIRMED',
                    message:
                      'Run creation could not be confirmed. Retry with the same validated snapshot.',
                  },
                ]
              : [],
          ),
        );
        return;
      }
      navigate(`/runs/${run.run_id}`, {
        state: {
          beginPreflight: true,
          preflightRequestId: crypto.randomUUID(),
        },
      });
    },
    onError: (error, attempt) => {
      if (attempt.revision !== draft.state.semanticRevision) {
        setSnapshotState(null);
        setIssues([
          {
            code: 'RUN_CREATE_INPUTS_CHANGED',
            message:
              'Inputs changed while run creation was running. The earlier request may have created a run; review canonical runs before retrying.',
          },
        ]);
        return;
      }
      const apiCode = error instanceof ApiError ? error.code : null;
      if (
        apiCode === 'VALIDATED_SNAPSHOT_EXPIRED' ||
        apiCode === 'VALIDATED_SNAPSHOT_STALE'
      ) {
        setSnapshotState((current) =>
          current?.snapshot.snapshot_id === attempt.snapshotId ? null : current,
        );
        onValidationFeedback((current) =>
          current !== null &&
          current.revision === attempt.revision &&
          current.snapshotId === attempt.snapshotId
            ? {
                ...current,
                notice:
                  'The validated snapshot is no longer valid. Validate the current draft again.',
              }
            : current,
        );
      }
      setIssues(
        requestIssues(error, {
          code: 'RUN_CREATE_UNCONFIRMED',
          message:
            'Run creation could not be confirmed. Retry with the same validated snapshot to read the canonical outcome.',
        }),
      );
    },
  });

  const expiryLabel = useMemo(() => {
    if (activeSnapshot === null) return null;
    return new Intl.DateTimeFormat(undefined, {
      dateStyle: 'medium',
      timeStyle: 'short',
    }).format(new Date(activeSnapshot.expires_at));
  }, [activeSnapshot]);
  const advisoryIssues = issues.filter(
    (issue) => issue.severity === 'warning' || issue.severity === 'info',
  );
  const errorIssues = issues.filter(
    (issue) => issue.severity !== 'warning' && issue.severity !== 'info',
  );

  return (
    <section
      className="mt-4 min-w-0 border-t border-[var(--color-border)] pt-4"
      aria-labelledby="validated-submission-title"
    >
      <ExecutionAvailabilityNotice availability={availability} />
      <div className="flex min-w-0 flex-col gap-3 sm:flex-row sm:items-start sm:justify-between">
        <div className="min-w-0">
          <h3 id="validated-submission-title" className="text-sm font-semibold">
            Backend validation and run creation
          </h3>
          <p className="mt-1 text-xs text-[var(--color-text-muted)]">
            The adapter is authoritative. A run can use only the exact draft saved by a successful validation.
          </p>
        </div>
        <div className="flex flex-wrap gap-2">
          <Button
            type="button"
            variant="secondary"
            className="gap-1.5"
            disabled={
              !draft.reviewReady ||
              !referenceSelectionAvailable ||
              draft.state.referenceProfileRevisionId === null ||
              validationMutation.isPending ||
              createMutation.isPending
            }
            onClick={() => {
              const referenceProfileRevisionId =
                draft.state.referenceProfileRevisionId;
              if (
                draft.reviewReady &&
                draft.review.ok &&
                referenceSelectionAvailable &&
                referenceProfileRevisionId !== null
              ) {
                onValidationFeedback({ revision: draft.state.semanticRevision, issues: [] });
                validationMutation.mutate({
                  requestId: ++latestValidationRequest.current,
                  payload: {
                    ...draft.review.payload,
                    reference_profile_revision_id: referenceProfileRevisionId,
                    ...(experimentRequest === null
                      ? {}
                      : {
                          project_id: experimentRequest.projectId,
                          sample_revision_ids: [...experimentRequest.sampleRevisionIds],
                          ...(experimentRequest.inputSelections.length === 0
                            ? {}
                            : {
                                input_selections: experimentRequest.inputSelections.map(
                                  (selection) => ({
                                    input_use_key: selection.input_use_key,
                                    occurrence: selection.occurrence,
                                    input_file_revision_ids: [
                                      ...selection.input_file_revision_ids,
                                    ],
                                  }),
                                ),
                              }),
                          ...(experimentRequest.experimentRevisionId === null
                            ? {}
                            : {
                                experiment_revision_id:
                                  experimentRequest.experimentRevisionId,
                              }),
                        }),
                  },
                  revision: draft.state.semanticRevision,
                  referenceProfileRevisionId,
                  contextKey,
                });
              }
            }}
            aria-label="Validate current inputs"
            data-testid="validate-draft-button"
          >
            <ShieldCheck aria-hidden="true" size={16} />
            {validationMutation.isPending ? 'Validating…' : 'Validate inputs'}
          </Button>
          <Button
            type="button"
            variant="primary"
            className="gap-1.5"
            disabled={
              !executionAvailable ||
              activeSnapshot === null ||
              validationMutation.isPending ||
              createMutation.isPending
            }
            onClick={() => {
              if (
                executionAvailable &&
                activeSnapshot !== null &&
                snapshotState !== null
              ) {
                createMutation.mutate({
                  snapshotId: activeSnapshot.snapshot_id,
                  revision: snapshotState.revision,
                });
              }
            }}
            aria-label="Create run from validated inputs"
            data-testid="create-validated-run-button"
          >
            <Play aria-hidden="true" size={16} />
            {createMutation.isPending ? 'Creating run…' : 'Create run'}
          </Button>
        </div>
      </div>

      {experimentRequest !== null && (
        <p className="mt-3 text-sm" role="status" data-testid="experiment-binding-status">
          {bindingStatusCopy(bindingOutcome)}
        </p>
      )}

      <div className="mt-3 min-h-12" aria-live="polite">
        {activeSnapshot !== null && expiryLabel !== null && (
          <p className="mb-2 break-words text-xs">
            First use expires {expiryLabel}. Snapshot{' '}
            <code className="break-all">{activeSnapshot.snapshot_id}</code>
          </p>
        )}
        {advisoryIssues.length > 0 && (
          <div
            className="rounded border border-amber-200 bg-[var(--color-warning-bg)] px-3 py-2 text-sm text-[var(--color-warning)]"
            role="status"
            data-testid="creation-advisories"
          >
            <ul className="space-y-1">
              {advisoryIssues.map((issue, index) => (
                <li key={`${issue.code}:${issue.path ?? ''}:${index}`} className="min-w-0">
                  <span className="font-mono text-xs">{issue.code}</span>{' '}
                  <span className="break-words [overflow-wrap:anywhere]">{issue.message}</span>
                  {issue.path && (
                    <span className="ml-1 text-xs text-[var(--color-text-muted)]">
                      ({issue.path})
                    </span>
                  )}
                  {issue.hint && (
                    <p className="mt-0.5 break-words text-xs">{issue.hint}</p>
                  )}
                </li>
              ))}
            </ul>
          </div>
        )}
        {errorIssues.length > 0 && (
          <div
            className="rounded border border-red-200 bg-red-50 px-3 py-2 text-sm text-[var(--color-error)]"
            role="alert"
          >
            <ul className="space-y-1">
              {errorIssues.map((issue, index) => (
                <li key={`${issue.code}:${issue.path ?? ''}:${index}`} className="min-w-0">
                  <span className="font-mono text-xs">{issue.code}</span>{' '}
                  <span className="break-words [overflow-wrap:anywhere]">{issue.message}</span>
                  {issue.path && (
                    <span className="ml-1 text-xs text-[var(--color-text-muted)]">
                      ({issue.path})
                    </span>
                  )}
                  {issue.hint && (
                    <p className="mt-0.5 break-words text-xs">{issue.hint}</p>
                  )}
                </li>
              ))}
            </ul>
          </div>
        )}
      </div>
    </section>
  );
}

function bindingStatusCopy(outcome: BindingOutcome): string {
  switch (outcome.kind) {
    case 'idle':
      return 'No experiment binding. Choosing a revision or opening this page does not create one.';
    case 'unresolved':
      return 'Not bound. This validation did not select an experiment revision.';
    case 'frozen':
      return `Bound to experiment revision ${outcome.experimentRevisionId} for this validation snapshot.`;
    case 'failed':
      return 'No experiment binding was saved.';
    case 'unconfirmed':
      return 'The validation result is unconfirmed. Retry to see whether an experiment binding was saved.';
    case 'missing':
      return 'The validation response did not include an experiment binding status.';
    default: {
      const unreachable: never = outcome;
      return unreachable;
    }
  }
}
