import type { WorkflowSummary } from '../../api/types';
import { Link } from 'react-router-dom';
import { ExecutionAvailabilityBadge } from '../workflow-detail/WorkflowAvailability';

interface WorkflowCatalogProps {
  workflows: WorkflowSummary[];
  selectedWorkflowId: string | null;
}

export function WorkflowCatalog({
  workflows,
  selectedWorkflowId,
}: WorkflowCatalogProps) {
  return (
    <div className="space-y-1">
      {workflows.map((workflow) => {
        const id = workflow.metadata.workflow_id;
        const isSelected = id === selectedWorkflowId;
        return (
          <Link
            key={id}
            to={`/workflows/${encodeURIComponent(id)}`}
            className={`block w-full rounded-[4px] border p-2 text-left transition-colors ${
              isSelected
                ? 'border-[var(--color-accent)] bg-[var(--color-surface-selected)]'
                : 'border-[var(--color-border)] bg-[var(--color-surface)] hover:bg-[var(--color-surface-subtle)]'
            }`}
          >
            <div className="flex min-w-0 flex-wrap items-start justify-between gap-2">
              <div className="min-w-0 text-sm font-medium text-[var(--color-text)]">
                {workflow.metadata.name}
              </div>
              <ExecutionAvailabilityBadge
                availability={workflow.availability.execution}
              />
            </div>
            <div className="truncate text-xs text-[var(--color-text-muted)]" title={id}>
              {id}
            </div>
          </Link>
        );
      })}
    </div>
  );
}
