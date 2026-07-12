import { useEffect, useMemo, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useQuery } from '@tanstack/react-query';
import { useClients } from '../../api/client-context';
import type {
  Issue,
  ValidateWorkflowResponse,
  WorkflowInputs,
  WorkflowSchema,
} from '../../api/types';
import { Panel } from '../../components/Panel';
import { AgentSidebar } from '../../features/agent-sidebar/AgentSidebar';
import { IssuePanel } from '../../features/issues-panel/IssuePanel';
import { RunProgressPanel } from '../../features/run-progress/RunProgressPanel';
import { ValidationWorkspace } from '../../features/validation-workspace/ValidationWorkspace';
import { WorkflowDetail } from '../../features/workflow-detail/WorkflowDetail';

interface WorkflowDetailPageProps {
  workflowId: string;
}

function createParseErrorIssue(path: 'config' | 'options', message: string): Issue {
  return {
    code: 'FRONTEND_INPUT_INVALID',
    message: `Invalid JSON for ${path}.`,
    severity: 'error',
    path,
    source: 'ui',
    technical_message: message,
    hint: 'Use valid JSON.',
    context: {},
  };
}

function createNotFoundIssue(workflowId: string): Issue {
  return {
    code: 'WORKFLOW_NOT_FOUND',
    message: 'Workflow was not found.',
    severity: 'error',
    path: 'workflow_id',
    source: 'registry',
    technical_message: null,
    hint: 'Select a workflow from the catalog.',
    context: { workflow_id: workflowId },
  };
}

export function WorkflowDetailPage({ workflowId }: WorkflowDetailPageProps) {
  const navigate = useNavigate();
  const { workflowClient, agentClient, runClient } = useClients();

  const {
    data: workflowsData,
    isLoading: workflowsLoading,
    error: workflowsError,
  } = useQuery({
    queryKey: ['workflows'],
    queryFn: () => workflowClient.listWorkflows(),
  });

  const {
    data: schemaData,
    isLoading: schemaLoading,
    error: schemaError,
  } = useQuery({
    queryKey: ['workflow', workflowId, 'schema'],
    queryFn: () => workflowClient.getWorkflowSchema(workflowId),
  });

  const workflow = useMemo(
    () =>
      workflowsData?.workflows.find(
        (w) => w.metadata.workflow_id === workflowId,
      ),
    [workflowsData, workflowId],
  );

  const schemaHints: WorkflowSchema | null = schemaData?.schema_hints ?? null;

  const [configText, setConfigText] = useState('{}');
  const [samplesText, setSamplesText] = useState('');
  const [optionsText, setOptionsText] = useState('{}');
  const [validationResult, setValidationResult] =
    useState<ValidateWorkflowResponse | null>(null);
  const [validatedInputs, setValidatedInputs] = useState<WorkflowInputs | null>(
    null,
  );
  const [loading, setLoading] = useState(false);
  const [agentDraftMessage, setAgentDraftMessage] = useState('');
  const [agentFocusedIssue, setAgentFocusedIssue] = useState<Issue | null>(null);

  useEffect(() => {
    setConfigText('{}');
    setSamplesText('');
    setOptionsText('{}');
    setValidationResult(null);
    setValidatedInputs(null);
    setAgentDraftMessage('');
    setAgentFocusedIssue(null);
  }, [workflowId]);

  const isNotFound =
    !workflowsLoading &&
    !schemaLoading &&
    workflowsData !== undefined &&
    schemaData !== undefined &&
    ((workflowsData.ok && workflow === undefined) ||
      schemaData.issues.some((issue) => issue.code === 'WORKFLOW_NOT_FOUND'));

  const loadError = workflowsError ?? schemaError;
  const loadIssue =
    workflowsData?.ok === false
      ? workflowsData.issues[0]
      : schemaData?.ok === false && !isNotFound
        ? schemaData.issues[0]
        : null;

  const notFoundIssues = useMemo(
    () => (isNotFound ? [createNotFoundIssue(workflowId)] : []),
    [isNotFound, workflowId],
  );

  const displayedIssues = useMemo(
    () =>
      isNotFound
        ? notFoundIssues
        : validationResult?.issues ?? [],
    [isNotFound, notFoundIssues, validationResult],
  );

  async function handleValidate() {
    setLoading(true);
    let config: Record<string, unknown> = {};
    let options: Record<string, unknown> = {};
    const issues: Issue[] = [];

    try {
      config = JSON.parse(configText);
    } catch (error) {
      issues.push(
        createParseErrorIssue(
          'config',
          error instanceof Error ? error.message : 'Unknown parse error',
        ),
      );
    }

    try {
      options = JSON.parse(optionsText);
    } catch (error) {
      issues.push(
        createParseErrorIssue(
          'options',
          error instanceof Error ? error.message : 'Unknown parse error',
        ),
      );
    }

    if (issues.length > 0) {
      setValidationResult({
        ok: false,
        workflow_id: workflowId,
        value: null,
        issues,
      });
      setValidatedInputs(null);
      setLoading(false);
      return;
    }

    const inputs: WorkflowInputs = {
      config,
      samples: samplesText || null,
      options,
    };

    try {
      const response = await workflowClient.validateWorkflow(workflowId, inputs);
      setValidationResult(response);
      if (response.ok) {
        setValidatedInputs(inputs);
      } else {
        setValidatedInputs(null);
      }
    } catch (error) {
      setValidationResult({
        ok: false,
        workflow_id: workflowId,
        value: null,
        issues: [
          {
            code: 'API_UNAVAILABLE',
            message:
              error instanceof Error ? error.message : 'Validation request failed.',
            severity: 'error',
            path: null,
            source: 'ui',
            technical_message: null,
            hint: null,
            context: {},
          },
        ],
      });
      setValidatedInputs(null);
    } finally {
      setLoading(false);
    }
  }

  function handleConfigChange(value: string) {
    setConfigText(value);
    setValidationResult(null);
    setValidatedInputs(null);
  }

  function handleSamplesChange(value: string) {
    setSamplesText(value);
    setValidationResult(null);
    setValidatedInputs(null);
  }

  function handleOptionsChange(value: string) {
    setOptionsText(value);
    setValidationResult(null);
    setValidatedInputs(null);
  }

  function handleAskAgent(issue: Issue) {
    setAgentDraftMessage(`Explain issue ${issue.code}.`);
    setAgentFocusedIssue(issue);
  }

  function handleAgentDraftConsumed() {
    setAgentDraftMessage('');
  }

  function handleAgentFocusedIssueConsumed() {
    setAgentFocusedIssue(null);
  }

  function handleRunCreated(runId: string) {
    navigate(`/runs/${runId}`, {
      state: {
        beginPreflight: true,
        preflightRequestId: crypto.randomUUID(),
      },
    });
  }

  if (isNotFound) {
    return (
      <section className="flex min-w-0 flex-1 flex-col gap-3">
        <Panel title="Workflow not found">
          <p className="text-sm text-[var(--color-error)]">
            No workflow with ID <code>{workflowId}</code> was found.
          </p>
        </Panel>
        <Panel title="Validation results">
          <IssuePanel issues={displayedIssues} onAskAgent={handleAskAgent} />
        </Panel>
      </section>
    );
  }

  if (loadError || loadIssue) {
    return (
      <section className="flex min-w-0 flex-1 flex-col gap-3">
        <Panel title="Workflow unavailable">
          <p className="text-sm text-[var(--color-error)]" role="alert">
            {loadError instanceof Error
              ? loadError.message
              : loadIssue?.message ?? 'Unable to load workflow details.'}
          </p>
        </Panel>
      </section>
    );
  }

  return (
    <>
      <section className="flex min-w-0 flex-1 flex-col gap-3">
        {workflow && (
          <Panel title="Workflow detail">
            <WorkflowDetail
              workflowId={workflow.metadata.workflow_id}
              workflowName={workflow.metadata.name}
              schemaHints={schemaHints}
            />
          </Panel>
        )}
        {workflow && (
          <Panel title="Validation workspace">
            <ValidationWorkspace
              workflowId={workflow.metadata.workflow_id}
              configText={configText}
              samplesText={samplesText}
              optionsText={optionsText}
              loading={loading}
              onConfigChange={handleConfigChange}
              onSamplesChange={handleSamplesChange}
              onOptionsChange={handleOptionsChange}
              onValidate={handleValidate}
            />
          </Panel>
        )}
        {workflow && (
          <Panel title="Run progress">
            <RunProgressPanel
              workflowId={workflow.metadata.workflow_id}
              validationResult={validationResult}
              validatedInputs={validatedInputs}
              runClient={runClient}
              onRunCreated={handleRunCreated}
            />
          </Panel>
        )}
        <Panel title="Validation results">
          <IssuePanel issues={displayedIssues} onAskAgent={handleAskAgent} />
        </Panel>
      </section>

      <aside className="flex w-full flex-col gap-3 lg:w-72">
        <AgentSidebar
          workflowId={workflowId}
          issues={displayedIssues}
          agentClient={agentClient}
          draftMessage={agentDraftMessage || undefined}
          onDraftConsumed={handleAgentDraftConsumed}
          focusedIssue={agentFocusedIssue}
          onFocusedIssueConsumed={handleAgentFocusedIssueConsumed}
        />
      </aside>
    </>
  );
}
