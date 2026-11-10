import { useQuery } from '@tanstack/react-query';
import { getWorkflowSchema } from '../../api/generated/workflows/workflows';
import { useClients } from '../../api/client-context';
import { ConfigurationDraftWorkbench } from '../../features/configuration-drafts/ConfigurationDraftWorkbench';
import { readWorkbenchSchema } from '../../features/input-workbench/schemaContract';

interface ConfigurationDraftsPageProps {
  workflowId: string;
}

export function ConfigurationDraftsPage({
  workflowId,
}: ConfigurationDraftsPageProps) {
  const { workflowClient } = useClients();
  const detailQuery = useQuery({
    queryKey: ['workflow', workflowId],
    queryFn: () => workflowClient.getWorkflow(workflowId),
    retry: false,
  });
  const schemaQuery = useQuery({
    queryKey: ['workflow-authoring-schema', workflowId],
    queryFn: () => getWorkflowSchema(workflowId),
    retry: false,
  });
  const workflowName =
    detailQuery.data?.ok === true
      ? detailQuery.data.workflow?.metadata.name ?? 'Workflow name unavailable'
      : detailQuery.isError || detailQuery.data?.ok === false
        ? 'Workflow name unavailable'
        : 'Loading workflow';
  const parsed = schemaQuery.data?.schema
    ? readWorkbenchSchema(schemaQuery.data.schema)
    : null;
  const schema = parsed?.ok ? parsed.value : null;
  const schemaMessage = schemaQuery.isLoading
    ? null
    : schemaQuery.isError
      ? 'The workflow editor could not be loaded.'
      : schemaQuery.data && !schemaQuery.data.schema
        ? 'This workflow does not publish an editor contract.'
        : parsed && !parsed.ok
          ? parsed.message
          : null;

  if (schemaQuery.isLoading) {
    return (
      <section
        className="min-h-64 min-w-0 flex-1 animate-pulse rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
        aria-label="Loading configuration drafts"
      >
        <div className="h-5 w-48 rounded-[4px] bg-[var(--color-skeleton)]" />
        <div className="mt-4 h-10 rounded-[4px] bg-[var(--color-skeleton)]" />
        <div className="mt-4 h-64 rounded-[4px] bg-[var(--color-skeleton)]" />
      </section>
    );
  }

  return (
    <ConfigurationDraftWorkbench
      workflowName={workflowName}
      projectName={null}
      schema={schema}
      schemaIdentity={
        schema
          ? `${workflowId}:${schema.contract.schema_version}`
          : workflowId
      }
      schemaMessage={schemaMessage}
      list={{
        status: 'unavailable',
        drafts: [],
        pageLimit: 20,
        hasPrevious: false,
        hasNext: false,
      }}
      selectedSid={null}
      draftExpired={false}
      connection={{ state: 'disconnected', draftRestored: false }}
      capability={null}
      proposal={null}
      explanation={null}
      validation={{ status: 'unavailable' }}
      requestStatus="idle"
      writesEnabled={false}
    />
  );
}
