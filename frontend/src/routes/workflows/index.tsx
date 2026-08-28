import { useQuery } from '@tanstack/react-query';
import { FilePenLine, Workflow } from 'lucide-react';
import { Link } from 'react-router-dom';
import { useClients } from '../../api/client-context';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { ExecutionAvailabilityBadge } from '../../features/workflow-detail/WorkflowAvailability';

export function WorkflowsIndexPage() {
  const { workflowClient } = useClients();
  const { data, isLoading, error } = useQuery({
    queryKey: ['workflows'],
    queryFn: () => workflowClient.listWorkflows(),
  });

  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Catalog">
        <div className="min-w-0">
          <header>
            <h2 className="flex items-center gap-2 text-xl font-semibold leading-7">
              <Workflow aria-hidden="true" size={20} />
              Workflows
            </h2>
            <p className="mt-1 text-sm text-[var(--color-text-muted)]">
              Choose an adapter-owned workflow, review its contract, or author inputs.
            </p>
          </header>
          {isLoading ? (
            <div className="mt-4 space-y-2" aria-label="Loading workflows" role="status">
              {[0, 1, 2].map((row) => (
                <div key={row} className="h-14 animate-pulse rounded-[4px] bg-[var(--color-skeleton)]" />
              ))}
            </div>
          ) : error || !data?.ok ? (
            <div
              className="mt-4 border-l-[3px] border-[var(--color-error)] bg-[var(--color-error-bg)] px-3 py-2 text-sm text-[var(--color-error)]"
              role="alert"
            >
              {error instanceof Error
                ? error.message
                : data?.issues[0]?.message ?? 'Unable to load workflows.'}
            </div>
          ) : data.workflows.length === 0 ? (
            <p className="mt-4 border-y border-[var(--color-border)] py-4 text-sm text-[var(--color-text-muted)]" role="status">
              No workflows are registered. Ask the local operator to review the installation.
            </p>
          ) : (
            <ul className="mt-4 divide-y divide-[var(--color-border)] border-y border-[var(--color-border)]">
              {data.workflows.map((workflow) => {
                const id = workflow.metadata.workflow_id;
                return (
                  <li key={id} className="grid min-w-0 gap-3 py-3 sm:grid-cols-[minmax(0,1fr)_auto] sm:items-center">
                    <div className="min-w-0">
                      <div className="flex min-w-0 flex-wrap items-center gap-2">
                        <Link
                          className="min-w-0 font-semibold text-[var(--color-link)] hover:underline"
                          to={`/workflows/${encodeURIComponent(id)}`}
                        >
                          {workflow.metadata.name}
                        </Link>
                        <ExecutionAvailabilityBadge availability={workflow.availability.execution} />
                      </div>
                      <p className="mt-1 truncate font-mono text-xs text-[var(--color-text-faint)]" title={id}>
                        {id}
                      </p>
                      <p className="mt-1 text-xs text-[var(--color-text-muted)]">
                        {workflow.metadata.engines.join(', ')} · Schema {workflow.schema_version}
                      </p>
                    </div>
                    <Button asChild variant="primary" className="w-full gap-1.5 sm:w-auto">
                      <Link to={`/workflows/${encodeURIComponent(id)}/new-run`}>
                        <FilePenLine aria-hidden="true" size={16} />
                        Author inputs
                      </Link>
                    </Button>
                  </li>
                );
              })}
            </ul>
          )}
        </div>
      </Panel>
    </section>
  );
}
