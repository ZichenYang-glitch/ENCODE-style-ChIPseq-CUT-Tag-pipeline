import type { Issue } from '../../api/types';

interface RunIssuePanelProps {
  issues: Issue[];
  title?: string;
}

export function RunIssuePanel({ issues, title = 'Run issue' }: RunIssuePanelProps) {
  if (issues.length === 0) return null;

  return (
    <section
      className="space-y-2 rounded border border-[var(--color-error)] bg-[var(--color-error-bg)] p-3"
      aria-label={title}
      role="alert"
      data-testid="run-progress-error"
    >
      {issues.map((issue, index) => (
        <div key={`${issue.code}-${index}`} className="text-sm text-[var(--color-error)]">
          <div className="font-semibold">
            {issue.code}: {issue.message}
          </div>
          {(issue.source || issue.path) && (
            <div className="mt-0.5 text-xs">
              {[issue.source, issue.path].filter(Boolean).join(' · ')}
            </div>
          )}
          {issue.hint && <div className="mt-1 text-xs">{issue.hint}</div>}
        </div>
      ))}
    </section>
  );
}
