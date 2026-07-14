import { ChevronRight, Code2 } from 'lucide-react';
import type { WorkflowSchema } from '../../api/types';

interface WorkflowDetailProps {
  workflowId: string;
  workflowName: string;
}

export function WorkflowDetail({
  workflowId,
  workflowName,
}: WorkflowDetailProps) {
  return (
    <div className="min-w-0 text-sm text-[var(--color-text-muted)]">
      <span className="font-medium text-[var(--color-text)]">{workflowName}</span>
      <span className="mx-1">·</span>
      <code className="break-all text-xs">{workflowId}</code>
    </div>
  );
}

interface DeveloperSchemaDetailsProps {
  workflowId: string;
  schemaHints: WorkflowSchema | null;
}

export function DeveloperSchemaDetails({
  workflowId,
  schemaHints,
}: DeveloperSchemaDetailsProps) {
  if (!schemaHints) {
    return (
      <p className="text-sm text-[var(--color-text-muted)]">
        No developer schema is available for {workflowId}.
      </p>
    );
  }

  return (
    <details
      className="group min-w-0 rounded border border-[var(--color-border)] bg-[var(--color-bg)]"
      data-testid="developer-schema-details"
    >
      <summary className="flex cursor-pointer list-none items-center gap-2 rounded px-3 py-2 text-sm font-medium text-[var(--color-text)] focus:outline-none focus:ring-2 focus:ring-inset focus:ring-[var(--color-accent)] [&::-webkit-details-marker]:hidden">
        <ChevronRight
          aria-hidden="true"
          className="shrink-0 transition-transform group-open:rotate-90"
          size={16}
        />
        <Code2 aria-hidden="true" className="shrink-0" size={16} />
        Developer schema
      </summary>
      <div className="grid min-w-0 grid-cols-1 gap-3 border-t border-[var(--color-border)] p-3 md:grid-cols-3">
        <SchemaSection title="Config schema" hints={schemaHints.config_schema} />
        <SchemaSection title="Sample schema" hints={schemaHints.sample_schema} />
        <SchemaSection title="Options schema" hints={schemaHints.option_schema} />
      </div>
    </details>
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
    <div className="flex min-w-0 flex-col">
      <h4 className="mb-1 text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
        {title}
      </h4>
      <pre className="max-h-96 w-full overflow-auto rounded border border-[var(--color-border)] bg-[var(--color-surface)] p-2 text-xs text-[var(--color-text)]">
        {JSON.stringify(hints, null, 2)}
      </pre>
    </div>
  );
}
