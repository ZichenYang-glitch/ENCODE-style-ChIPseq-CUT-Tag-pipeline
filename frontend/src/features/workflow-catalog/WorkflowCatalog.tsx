import type { WorkflowSummary } from '../../api/types';

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
    <div className="space-y-2">
      {workflows.map((workflow) => {
        const id = workflow.metadata.workflow_id;
        const isSelected = id === selectedWorkflowId;
        return (
          <button
            key={id}
            onClick={() => onSelect(id)}
            className={`w-full rounded border p-3 text-left transition-colors ${
              isSelected
                ? 'border-[var(--color-accent)] bg-[var(--color-info-bg)]'
                : 'border-[var(--color-border)] bg-[var(--color-bg)] hover:bg-[var(--color-surface)]'
            }`}
          >
            <div className="font-medium text-[var(--color-text)]">
              {workflow.metadata.name}
            </div>
            <div className="text-xs text-[var(--color-text-muted)]">{id}</div>
          </button>
        );
      })}
    </div>
  );
}
