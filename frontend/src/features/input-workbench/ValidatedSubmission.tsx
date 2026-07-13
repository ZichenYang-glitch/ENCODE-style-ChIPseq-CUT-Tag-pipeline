import { useEffect, useMemo, useState } from 'react';
import { useMutation } from '@tanstack/react-query';
import { CheckCircle2, CircleAlert, Play, ShieldCheck } from 'lucide-react';
import { useNavigate } from 'react-router-dom';
import { createRun } from '../../api/generated/runs/runs';
import { validateWorkflow } from '../../api/generated/workflows/workflows';
import type {
  IssueResponse,
  ValidationRequest,
  ValidatedInputSnapshotResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import type { InputDraftController } from './useInputDraft';

interface SafeIssue {
  code: string;
  message: string;
  path?: string | null;
  hint?: string | null;
}

interface SnapshotState {
  snapshot: ValidatedInputSnapshotResponse;
  revision: number;
}

interface ValidationAttempt {
  payload: ValidationRequest;
  revision: number;
}

interface ValidatedSubmissionProps {
  workflowId: string;
  draft: InputDraftController;
}

function safeIssues(issues: IssueResponse[] | undefined): SafeIssue[] {
  return (issues ?? []).map((issue) => ({
    code: issue.code,
    message: issue.message,
    path: issue.path,
    hint: issue.hint,
  }));
}

function requestIssues(error: unknown, fallback: SafeIssue): SafeIssue[] {
  if (error instanceof ApiError && error.issues.length > 0) {
    return error.issues.map((issue) => ({
      code: issue.code,
      message: issue.message,
      path: issue.path,
      hint: issue.hint,
    }));
  }
  return [fallback];
}

export function ValidatedSubmission({
  workflowId,
  draft,
}: ValidatedSubmissionProps) {
  const navigate = useNavigate();
  const [snapshotState, setSnapshotState] = useState<SnapshotState | null>(null);
  const [issues, setIssues] = useState<SafeIssue[]>([]);
  const [notice, setNotice] = useState<string | null>(null);
  const activeSnapshot =
    snapshotState?.revision === draft.state.semanticRevision
      ? snapshotState.snapshot
      : null;
  const snapshotInvalidated = snapshotState !== null && activeSnapshot === null;

  useEffect(() => {
    if (snapshotInvalidated) {
      setNotice('Inputs changed after validation. Validate the current draft again.');
      setIssues([]);
    }
  }, [snapshotInvalidated]);

  const validationMutation = useMutation({
    mutationFn: (attempt: ValidationAttempt) =>
      validateWorkflow(workflowId, attempt.payload),
    onSuccess: (response, attempt) => {
      if (attempt.revision !== draft.state.semanticRevision) {
        setSnapshotState(null);
        setIssues([]);
        setNotice(
          'Inputs changed while validation was running. Validate the current draft again.',
        );
        return;
      }
      const responseIssues = safeIssues(response.issues);
      if (
        !response.ok ||
        response.snapshot === null ||
        response.snapshot.workflow_id !== workflowId
      ) {
        setSnapshotState(null);
        setIssues(
          responseIssues.length > 0
            ? responseIssues
            : [
                {
                  code: 'VALIDATION_NOT_CONFIRMED',
                  message: 'Backend validation did not return a usable snapshot.',
                },
              ],
        );
        setNotice(null);
        return;
      }
      setSnapshotState({
        snapshot: response.snapshot,
        revision: draft.state.semanticRevision,
      });
      setIssues(responseIssues);
      setNotice('Backend validation succeeded. This exact draft can create one run.');
    },
    onError: (error, attempt) => {
      setSnapshotState(null);
      if (attempt.revision !== draft.state.semanticRevision) {
        setIssues([]);
        setNotice(
          'Inputs changed while validation was running. Validate the current draft again.',
        );
        return;
      }
      setIssues(
        requestIssues(error, {
          code: 'VALIDATION_UNAVAILABLE',
          message: 'Validation could not be confirmed. Retry when the API is available.',
        }),
      );
      setNotice(null);
    },
  });

  const createMutation = useMutation({
    mutationFn: async () => {
      if (activeSnapshot === null) {
        throw new Error('validated snapshot is unavailable');
      }
      return createRun(workflowId, {
        snapshot_id: activeSnapshot.snapshot_id,
      });
    },
    onSuccess: (response) => {
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
    onError: (error) => {
      const apiCode = error instanceof ApiError ? error.code : null;
      if (
        apiCode === 'VALIDATED_SNAPSHOT_EXPIRED' ||
        apiCode === 'VALIDATED_SNAPSHOT_STALE'
      ) {
        setSnapshotState(null);
      }
      setIssues(
        requestIssues(error, {
          code: 'RUN_CREATE_UNCONFIRMED',
          message:
            'Run creation could not be confirmed. Retry with the same validated snapshot to read the canonical outcome.',
        }),
      );
      setNotice(null);
    },
  });

  const expiryLabel = useMemo(() => {
    if (activeSnapshot === null) return null;
    return new Intl.DateTimeFormat(undefined, {
      dateStyle: 'medium',
      timeStyle: 'short',
    }).format(new Date(activeSnapshot.expires_at));
  }, [activeSnapshot]);

  return (
    <section
      className="mt-4 min-w-0 border-t border-[var(--color-border)] pt-4"
      aria-labelledby="validated-submission-title"
    >
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
            disabled={!draft.reviewReady || validationMutation.isPending || createMutation.isPending}
            onClick={() => {
              if (draft.reviewReady && draft.review.ok) {
                validationMutation.mutate({
                  payload: draft.review.payload,
                  revision: draft.state.semanticRevision,
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
            disabled={activeSnapshot === null || validationMutation.isPending || createMutation.isPending}
            onClick={() => createMutation.mutate()}
            aria-label="Create run from validated inputs"
            data-testid="create-validated-run-button"
          >
            <Play aria-hidden="true" size={16} />
            {createMutation.isPending ? 'Creating run…' : 'Create run'}
          </Button>
        </div>
      </div>

      <div className="mt-3 min-h-12" aria-live="polite">
        {notice !== null && (
          <div
            className={`flex items-start gap-2 rounded border px-3 py-2 text-sm ${
              activeSnapshot !== null
                ? 'border-emerald-200 bg-emerald-50 text-emerald-800'
                : 'border-amber-200 bg-[var(--color-warning-bg)] text-[var(--color-warning)]'
            }`}
            role="status"
          >
            {activeSnapshot !== null ? (
              <CheckCircle2 className="mt-0.5 shrink-0" aria-hidden="true" size={16} />
            ) : (
              <CircleAlert className="mt-0.5 shrink-0" aria-hidden="true" size={16} />
            )}
            <div className="min-w-0">
              <p>{notice}</p>
              {activeSnapshot !== null && expiryLabel !== null && (
                <p className="mt-1 break-words text-xs">
                  First use expires {expiryLabel}. Snapshot{' '}
                  <code className="break-all">{activeSnapshot.snapshot_id}</code>
                </p>
              )}
            </div>
          </div>
        )}
        {issues.length > 0 && (
          <div
            className="rounded border border-red-200 bg-red-50 px-3 py-2 text-sm text-[var(--color-error)]"
            role="alert"
          >
            <ul className="space-y-1">
              {issues.map((issue, index) => (
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
