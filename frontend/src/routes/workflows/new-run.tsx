import { useQuery } from '@tanstack/react-query';
import { RotateCcw } from 'lucide-react';
import { getWorkflowSchema } from '../../api/generated/workflows/workflows';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { InputWorkbench } from '../../features/input-workbench/InputWorkbench';
import { readWorkbenchSchema } from '../../features/input-workbench/schemaContract';

interface NewRunWorkbenchPageProps {
  workflowId: string;
}

export function NewRunWorkbenchPage({ workflowId }: NewRunWorkbenchPageProps) {
  const query = useQuery({
    queryKey: ['workflow-authoring-schema', workflowId],
    queryFn: () => getWorkflowSchema(workflowId),
    retry: false,
  });

  if (query.isLoading) {
    return (
      <section
        data-testid="input-workbench-loading"
        className="min-h-[38rem] min-w-0 flex-1 animate-pulse rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
        aria-label="Loading input workbench"
      >
        <div className="h-5 w-48 rounded bg-slate-200" />
        <div className="mt-4 h-10 rounded bg-slate-100" />
        <div className="mt-4 h-96 rounded bg-slate-100" />
      </section>
    );
  }

  if (query.isError) {
    return (
      <WorkbenchError
        message="Unable to load the workflow authoring contract."
        onRetry={() => void query.refetch()}
      />
    );
  }

  const response = query.data;
  if (!response?.schema) {
    return (
      <WorkbenchError
        message="This workflow does not publish an authoring contract."
        onRetry={() => void query.refetch()}
      />
    );
  }
  const parsed = readWorkbenchSchema(response.schema);
  if (!parsed.ok) {
    return (
      <WorkbenchError
        message={parsed.message}
        onRetry={() => void query.refetch()}
      />
    );
  }

  return (
    <InputWorkbench
      key={`${workflowId}:${parsed.value.contract.schema_version}`}
      workflowId={workflowId}
      schema={parsed.value}
    />
  );
}

function WorkbenchError({
  message,
  onRetry,
}: {
  message: string;
  onRetry: () => void;
}) {
  return (
    <section className="min-h-[24rem] min-w-0 flex-1">
      <Panel title="Input workbench unavailable">
        <div className="space-y-3">
          <p className="text-sm text-[var(--color-error)]" role="alert">
            {message}
          </p>
          <Button type="button" className="gap-1.5" onClick={onRetry}>
            <RotateCcw aria-hidden="true" size={15} />
            <span aria-label="Retry schema">Retry</span>
          </Button>
        </div>
      </Panel>
    </section>
  );
}
