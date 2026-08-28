import { Outlet, useLocation, useParams } from 'react-router-dom';
import { useQuery } from '@tanstack/react-query';
import { Panel } from '../../components/Panel';
import { WorkflowCatalog } from '../../features/workflow-catalog/WorkflowCatalog';
import { useClients } from '../../api/client-context';

export function WorkflowsLayout() {
  const { pathname } = useLocation();
  const { workflowId } = useParams<{ workflowId: string }>();
  const showCatalog =
    workflowId !== undefined &&
    pathname === `/workflows/${encodeURIComponent(workflowId)}`;
  const { workflowClient } = useClients();
  const { data, isLoading, error } = useQuery({
    queryKey: ['workflows'],
    queryFn: () => workflowClient.listWorkflows(),
    enabled: showCatalog,
  });

  return (
    <>
      {showCatalog && (
        <aside className="w-full shrink-0 lg:w-56">
          <Panel title="Workflows">
            {isLoading ? (
              <p className="text-sm text-[var(--color-text-muted)]">Loading…</p>
            ) : error || !data?.ok ? (
              <p className="text-sm text-[var(--color-error)]" role="alert">
                {error instanceof Error
                  ? error.message
                  : data?.issues[0]?.message ?? 'Unable to load workflows.'}
              </p>
            ) : (
              <WorkflowCatalog
                workflows={data.workflows}
                selectedWorkflowId={workflowId ?? null}
              />
            )}
          </Panel>
        </aside>
      )}
      <div className="flex min-w-0 flex-1 flex-col gap-3 lg:flex-row">
        <Outlet />
      </div>
    </>
  );
}
