import type { WorkflowSchema } from '../../api/types';

interface WorkflowDetailProps {
  workflowId: string;
  workflowName: string;
  schemaHints: WorkflowSchema | null;
}

export function WorkflowDetail({
  workflowId,
  workflowName,
  schemaHints,
}: WorkflowDetailProps) {
  if (!schemaHints) {
    return (
      <p className="text-sm text-[var(--color-text-muted)]">
        No schema hints available for {workflowId}.
      </p>
    );
  }

  return (
    <div className="space-y-3">
      <div className="text-sm text-[var(--color-text-muted)]">
        <span className="font-medium text-[var(--color-text)]">{workflowName}</span>
        <span className="mx-1">·</span>
        <code className="text-xs">{workflowId}</code>
      </div>
      <div className="grid grid-cols-1 gap-3 md:grid-cols-3">
        <SchemaSection title="Config schema" hints={schemaHints.config_schema} />
        <SchemaSection title="Sample schema" hints={schemaHints.sample_schema} />
        <SchemaSection title="Options schema" hints={schemaHints.option_schema} />
      </div>
    </div>
  );
}

function SchemaSection({
  title,
  hints,
}: {
  title: string;
  hints: Record<string, unknown>;
}) {
  return (
    <div className="flex flex-col">
      <h4 className="mb-1 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
        {title}
      </h4>
      <pre className="overflow-auto rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-2 text-xs text-[var(--color-text)]">
        {JSON.stringify(hints, null, 2)}
      </pre>
    </div>
  );
}
