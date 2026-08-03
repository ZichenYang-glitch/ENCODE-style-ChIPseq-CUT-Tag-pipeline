import { useQuery } from '@tanstack/react-query';
import { RotateCcw } from 'lucide-react';
import {
  getWorkflowSchema,
  listCompatibleReferenceProfiles,
} from '../../api/generated/workflows/workflows';
import type { ReferenceProfileListResponse } from '../../api/generated/models';
import { useClients } from '../../api/client-context';
import {
  readReferenceProfileSummary,
  type ReferenceProfileSummary,
} from '../../api/runTypes';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { InputWorkbench } from '../../features/input-workbench/InputWorkbench';
import { readWorkbenchSchema } from '../../features/input-workbench/schemaContract';

interface NewRunWorkbenchPageProps {
  workflowId: string;
}

export function NewRunWorkbenchPage({ workflowId }: NewRunWorkbenchPageProps) {
  const { workflowClient } = useClients();
  const query = useQuery({
    queryKey: ['workflow-authoring-schema', workflowId],
    queryFn: () => getWorkflowSchema(workflowId),
    retry: false,
  });
  const detailQuery = useQuery({
    queryKey: ['workflow', workflowId],
    queryFn: () => workflowClient.getWorkflow(workflowId),
    retry: false,
  });
  const referenceQuery = useQuery({
    queryKey: ['workflow-reference-profiles', workflowId],
    queryFn: () => listCompatibleReferenceProfiles(workflowId),
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
  const references = readCompatibleProfiles(referenceQuery.data, workflowId);

  return (
    <InputWorkbench
      key={`${workflowId}:${parsed.value.contract.schema_version}`}
      workflowId={workflowId}
      schema={parsed.value}
      availability={detailQuery.data?.workflow?.availability ?? null}
      referenceProfiles={references.profiles}
      referenceProfilesLoading={referenceQuery.isLoading}
      referenceProfilesRefreshing={referenceQuery.isFetching}
      referenceProfilesError={
        referenceQuery.isError
          ? 'Unable to load enabled reference profiles.'
          : references.error
      }
      onRefreshReferenceProfiles={() => void referenceQuery.refetch()}
    />
  );
}

function readCompatibleProfiles(
  response: ReferenceProfileListResponse | undefined,
  workflowId: string,
): { profiles: ReferenceProfileSummary[] | null; error: string | null } {
  if (response === undefined) return { profiles: null, error: null };
  if (!response.ok) {
    return {
      profiles: [],
      error: 'Enabled reference profiles are unavailable for this workflow.',
    };
  }
  if (response.workflow_id !== workflowId || !Array.isArray(response.profiles)) {
    return {
      profiles: [],
      error: 'The reference profile response could not be verified.',
    };
  }
  const profiles: ReferenceProfileSummary[] = [];
  for (const candidate of response.profiles) {
    const profile = readReferenceProfileSummary(candidate);
    if (profile === null) {
      return {
        profiles: [],
        error: 'The reference profile response could not be verified.',
      };
    }
    profiles.push(profile);
  }
  if (
    new Set(profiles.map((profile) => profile.revision_id)).size !==
    profiles.length
  ) {
    return {
      profiles: [],
      error: 'The reference profile response could not be verified.',
    };
  }
  return { profiles, error: null };
}

function WorkbenchError({
  message,
  onRetry,
}: {
  message: string;
  onRetry: () => void;
}) {
  return (
    <section className="min-h-[38rem] min-w-0 flex-1">
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
