import type { WorkflowSummary } from '../../api/types';
import { ExecutionAvailabilityBadge } from '../workflow-detail/WorkflowAvailability';

interface WorkflowCatalogProps {
  workflows: WorkflowSummary[];
  selectedWorkflowId: string | null;
  onSelect: (workflowId: string) => void;
}

export function WorkflowCatalog({
  workflows,
  selectedWorkflowId,
  onSelect,
}: WorkflowCatalogProps) {
  return (
    <div className="space-y-1">
      {workflows.map((workflow) => {
        const id = workflow.metadata.workflow_id;
        const isSelected = id === selectedWorkflowId;
        return (
          <button
            key={id}
            onClick={() => onSelect(id)}
            className={`w-full rounded border p-2 text-left transition-colors focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)] ${
              isSelected
                ? 'border-[var(--color-accent)] bg-[var(--color-info-bg)]'
                : 'border-[var(--color-border)] bg-[var(--color-surface)] hover:bg-[var(--color-bg)]'
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
            <div className="break-all text-xs text-[var(--color-text-muted)]">
              {id}
            </div>
          </button>
        );
      })}
    </div>
  );
}
