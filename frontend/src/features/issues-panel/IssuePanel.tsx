import type { Issue, Severity } from '../../api/types';
import { Badge } from '../../components/Badge';

interface IssuePanelProps {
  issues: Issue[];
  onAskAgent?: (issue: Issue) => void;
}

const severityRank: Record<Severity, number> = {
  error: 0,
  warning: 1,
  info: 2,
};

function groupIssues(issues: Issue[]): Record<string, Issue[]> {
  const grouped: Record<string, Issue[]> = {};
  for (const issue of issues) {
    const source = issue.source ?? 'unknown';
    if (!grouped[source]) {
      grouped[source] = [];
    }
    grouped[source].push(issue);
  }
  for (const source of Object.keys(grouped)) {
    grouped[source].sort((a, b) => severityRank[a.severity] - severityRank[b.severity]);
  }
  return grouped;
}

function countBySeverity(issues: Issue[]) {
  return issues.reduce(
    (acc, issue) => {
      acc[issue.severity] += 1;
      return acc;
    },
    { error: 0, warning: 0, info: 0 },
  );
}

export function IssuePanel({ issues, onAskAgent }: IssuePanelProps) {
  const grouped = groupIssues(issues);
  const counts = countBySeverity(issues);

  if (issues.length === 0) {
    return (
      <div className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-3 text-sm text-[var(--color-text-muted)]">
        No validation issues yet. Run validation to review structured feedback.
      </div>
    );
  }

  return (
    <div className="space-y-3">
      <div className="flex flex-wrap gap-3 text-xs">
        <span className="font-medium text-[var(--color-text)]">Issues:</span>
        {counts.error > 0 && (
          <span className="font-medium text-[var(--color-error)]">{counts.error} errors</span>
        )}
        {counts.warning > 0 && (
          <span className="font-medium text-[var(--color-warning)]">{counts.warning} warnings</span>
        )}
        {counts.info > 0 && (
          <span className="font-medium text-[var(--color-info)]">{counts.info} info</span>
        )}
      </div>
      {Object.entries(grouped).map(([source, sourceIssues]) => (
        <div key={source}>
          <h4 className="mb-2 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
            {source}
          </h4>
          <ul className="space-y-2">
            {sourceIssues.map((issue, index) => (
              <li
                key={`${issue.code}-${issue.path ?? 'null'}-${index}`}
                className={`rounded border border-[var(--color-border)] bg-[var(--color-surface)] p-3 shadow-sm ${severityBorder(issue.severity)}`}
              >
                <div className="flex flex-wrap items-center gap-2">
                  <Badge severity={issue.severity} />
                  {issue.path && (
                    <code className="rounded bg-[var(--color-bg)] px-1.5 py-0.5 text-xs text-[var(--color-text-muted)]">
                      {issue.path}
                    </code>
                  )}
                  {issue.code && (
                    <code className="rounded bg-[var(--color-bg)] px-1.5 py-0.5 text-xs text-[var(--color-text-muted)]">
                      {issue.code}
                    </code>
                  )}
                </div>
                <p className="mt-2 text-sm font-medium text-[var(--color-text)]">{issue.message}</p>
                {issue.hint && (
                  <p className="mt-1 text-sm text-[var(--color-text-muted)]">
                    {issue.hint}
                  </p>
                )}
                {onAskAgent && (
                  <div className="mt-2">
                    <button
                      type="button"
                      onClick={() => onAskAgent(issue)}
                      className="text-xs text-[var(--color-accent)] hover:underline focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
                      aria-label={`Ask Agent about ${issue.code}`}
                    >
                      Ask Agent
                    </button>
                  </div>
                )}
                {issue.technical_message && (
                  <details className="mt-2">
                    <summary className="cursor-pointer text-xs text-[var(--color-accent)] hover:underline focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]">Details</summary>
                    <pre className="mt-2 overflow-auto rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-2 text-xs text-[var(--color-text)]">
                      {issue.technical_message}
                    </pre>
                  </details>
                )}
              </li>
            ))}
          </ul>
        </div>
      ))}
    </div>
  );
}

function severityBorder(severity: Severity): string {
  const borderColors = {
    error: 'border-l-4 border-l-[var(--color-error)]',
    warning: 'border-l-4 border-l-[var(--color-warning)]',
    info: 'border-l-4 border-l-[var(--color-info)]',
  };
  return borderColors[severity];
}
