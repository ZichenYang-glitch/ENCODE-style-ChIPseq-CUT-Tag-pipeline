import type { Issue } from '../../api/types';

interface AgentSidebarProps {
  workflowId: string | null;
  issues: Issue[];
}

export function AgentSidebar({ workflowId, issues }: AgentSidebarProps) {
  return (
    <div className="rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] p-3 shadow-sm">
      <header className="mb-3 border-b border-[var(--color-border)] pb-2">
        <h3 className="text-sm font-semibold text-[var(--color-text)]">
          Validation Assistant — Read Only
        </h3>
        <p className="mt-1 text-xs text-[var(--color-text-muted)]">
          I can help you understand validation errors and schema requirements. I
          cannot run, modify, or delete workflows.
        </p>
      </header>
      <div className="space-y-2 text-xs text-[var(--color-text-muted)]">
        <div>Workflow: {workflowId ?? 'None selected'}</div>
        <div>Issues: {issues.length}</div>
      </div>
    </div>
  );
}
