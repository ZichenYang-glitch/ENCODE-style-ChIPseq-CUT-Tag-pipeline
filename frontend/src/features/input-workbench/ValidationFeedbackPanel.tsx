import type { RJSFSchema } from '@rjsf/utils';
import { Button } from '../../components/Button';
import { configIssueLocation, type ValidationFeedback } from './validationFeedback';

interface Props {
  feedback: ValidationFeedback | null;
  revision: number;
  schema: RJSFSchema;
  onLocate: (path: string) => void;
  onReview: () => void;
}

export function ValidationFeedbackPanel({ feedback, revision, schema, onLocate, onReview }: Props) {
  if (feedback === null) return null;
  const stale = feedback.revision !== revision;
  const issues = stale ? [] : feedback.issues;
  const notice = stale
    ? feedback.changedInFlight
      ? 'Inputs changed while validation was running. Validate the current draft again.'
      : 'Inputs changed after validation. Validate the current draft again.'
    : feedback.notice;
  return (
    <div className="my-3 min-w-0 space-y-2" aria-live="polite" data-testid="validation-feedback">
      {notice && (
        <div role="status" className="rounded border border-[var(--color-border)] px-3 py-2 text-sm">
          <p>{notice}</p>
          {stale && <Button type="button" variant="secondary" className="mt-2" onClick={onReview}>Review current inputs</Button>}
        </div>
      )}
      {(['error', 'advisory'] as const).map((kind) => {
        const selected = issues.filter((issue) =>
          (issue.severity === 'warning' || issue.severity === 'info') === (kind === 'advisory'));
        return selected.length > 0 && (
          <div key={kind} role={kind === 'error' ? 'alert' : 'status'}
            data-testid={kind === 'error' ? 'validation-errors' : 'validation-advisories'}
            className="min-w-0 rounded border border-[var(--color-border)] bg-[var(--color-surface)] px-3 py-2 text-sm">
            <p className="font-semibold">{kind === 'error' ? 'Inputs need attention' : 'Validation notes'}</p>
            <ul className="mt-2 space-y-3">
              {selected.map((issue, index) => {
                const location = configIssueLocation(schema, issue.path);
                return (
                  <li key={`${issue.code}:${index}`} className="min-w-0 break-words [overflow-wrap:anywhere]">
                    <span className="font-mono text-xs">{issue.code}</span>{' '}{issue.message}
                    {issue.path && <p className="text-xs">{issue.path}</p>}
                    {issue.hint && <p className="mt-1">{issue.hint}</p>}
                    {location && (
                      <Button type="button" variant="secondary" className="mt-2 max-w-full whitespace-normal text-left"
                        onClick={() => onLocate(location.path)}>Go to {location.path}</Button>
                    )}
                  </li>
                );
              })}
            </ul>
          </div>
        );
      })}
    </div>
  );
}
