import { ChevronRight, Code2 } from 'lucide-react';
import type { WorkflowSchema, WorkflowSummary } from '../../api/types';
import { ExecutionAvailabilityBadge } from './WorkflowAvailability';

interface WorkflowDetailProps {
  workflow: WorkflowSummary;
}

export function WorkflowDetail({ workflow }: WorkflowDetailProps) {
  const { metadata } = workflow;
  return (
    <div className="min-w-0 space-y-3 text-sm text-[var(--color-text-muted)]">
      <div className="min-w-0">
        <div className="flex min-w-0 flex-wrap items-center gap-2">
          <span className="font-medium text-[var(--color-text)]">
            {metadata.name}
          </span>
          <ExecutionAvailabilityBadge
            availability={workflow.availability.execution}
          />
        </div>
        <code className="mt-1 block break-all text-xs">
          {metadata.workflow_id}
        </code>
        {metadata.description && (
          <p className="mt-2 break-words">{metadata.description}</p>
        )}
      </div>
      <div className="flex min-w-0 flex-wrap gap-2 text-xs">
        <span className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2 py-1">
          Adapter {metadata.version}
        </span>
        <span className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2 py-1">
          Schema {workflow.schema_version}
        </span>
        {metadata.engines.map((engine) => (
          <span
            key={engine}
            className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2 py-1"
          >
            {engine}
          </span>
        ))}
      </div>
      {workflow.upstream_identity && (
        <p className="break-words text-xs">
          Upstream: {workflow.upstream_identity.name}{' '}
          {workflow.upstream_identity.version}
          <span className="mx-1">·</span>
          <code className="break-all">
            {workflow.upstream_identity.revision}
          </code>
        </p>
      )}
      <p className="min-w-0 break-words text-xs">
        Input authoring is available. Execution status:{' '}
        <code className="break-all">{workflow.availability.reason_code}</code>
      </p>
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
