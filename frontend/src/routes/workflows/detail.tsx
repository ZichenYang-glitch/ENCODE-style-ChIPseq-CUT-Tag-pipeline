import { useEffect, useMemo, useRef, useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { FilePenLine } from 'lucide-react';
import { useQuery } from '@tanstack/react-query';
import { useClients } from '../../api/client-context';
import type {
  Issue,
  ValidateWorkflowResponse,
  WorkflowInputs,
  WorkflowSchema,
} from '../../api/types';
import { Panel } from '../../components/Panel';
import { Button } from '../../components/Button';
import { AgentSidebar } from '../../features/agent-sidebar/AgentSidebar';
import { IssuePanel } from '../../features/issues-panel/IssuePanel';
import { RunProgressPanel } from '../../features/run-progress/RunProgressPanel';
import { ValidationWorkspace } from '../../features/validation-workspace/ValidationWorkspace';
import {
  DeveloperSchemaDetails,
  WorkflowDetail,
} from '../../features/workflow-detail/WorkflowDetail';

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
    data: workflowData,
    isLoading: workflowLoading,
    error: workflowError,
  } = useQuery({
    queryKey: ['workflow', workflowId],
    queryFn: () => workflowClient.getWorkflow(workflowId),
  });

  const {
    data: schemaData,
    isLoading: schemaLoading,
    error: schemaError,
  } = useQuery({
    queryKey: ['workflow', workflowId, 'schema'],
    queryFn: () => workflowClient.getWorkflowSchema(workflowId),
  });

  const workflow = workflowData?.workflow ?? null;

  const schemaHints: WorkflowSchema | null = schemaData?.schema_hints ?? null;
  const supportsServerPathSamples =
    schemaHints?.input_modes?.samples.includes('server_path') === true;

  const [configText, setConfigText] = useState('{}');
  const [samplesText, setSamplesText] = useState('');
  const [optionsText, setOptionsText] = useState('{}');
  const [validationResult, setValidationResult] =
    useState<ValidateWorkflowResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const validationAttemptRef = useRef(0);
  const workflowIdRef = useRef(workflowId);
  workflowIdRef.current = workflowId;
  const [agentDraftMessage, setAgentDraftMessage] = useState('');
  const [agentFocusedIssue, setAgentFocusedIssue] = useState<Issue | null>(null);

  useEffect(() => {
    validationAttemptRef.current += 1;
    setConfigText('{}');
    setSamplesText('');
    setOptionsText('{}');
    setValidationResult(null);
    setLoading(false);
    setAgentDraftMessage('');
    setAgentFocusedIssue(null);
  }, [workflowId]);

  const isNotFound =
    !workflowLoading &&
    !schemaLoading &&
    workflowData !== undefined &&
    schemaData !== undefined &&
    ((!workflowData.ok &&
      workflowData.issues.some((issue) => issue.code === 'WORKFLOW_NOT_FOUND')) ||
      schemaData.issues.some((issue) => issue.code === 'WORKFLOW_NOT_FOUND'));

  const loadError = workflowError ?? schemaError;
  const loadIssue =
    workflowData?.ok === false && !isNotFound
      ? workflowData.issues[0]
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
    const attempt = validationAttemptRef.current + 1;
    validationAttemptRef.current = attempt;
    const isCurrentAttempt = () =>
      validationAttemptRef.current === attempt &&
      workflowIdRef.current === workflowId;
    setLoading(true);
    let config: WorkflowInputs['config'] = {};
    let options: NonNullable<WorkflowInputs['options']> = {};
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
      if (isCurrentAttempt()) {
        setValidationResult({
          ok: false,
          workflow_id: workflowId,
          value: null,
          snapshot: null,
          issues,
        });
        setLoading(false);
      }
      return;
    }

    const inputs: WorkflowInputs = {
      config,
      samples: samplesText || null,
      options,
    };

    try {
      const response = await workflowClient.validateWorkflow(workflowId, inputs);
      if (isCurrentAttempt()) setValidationResult(response);
    } catch {
      if (isCurrentAttempt()) {
        setValidationResult({
          ok: false,
          workflow_id: workflowId,
          value: null,
          snapshot: null,
          issues: [
            {
              code: 'API_UNAVAILABLE',
              message: 'Validation request could not be confirmed. Try again.',
              severity: 'error',
              path: null,
              source: 'ui',
              technical_message: null,
              hint: null,
              context: {},
            },
          ],
        });
      }
    } finally {
      if (isCurrentAttempt()) setLoading(false);
    }
  }

  function handleConfigChange(value: string) {
    validationAttemptRef.current += 1;
    setConfigText(value);
    setValidationResult(null);
    setLoading(false);
  }

  function handleSamplesChange(value: string) {
    validationAttemptRef.current += 1;
    setSamplesText(value);
    setValidationResult(null);
    setLoading(false);
  }

  function handleOptionsChange(value: string) {
    validationAttemptRef.current += 1;
    setOptionsText(value);
    setValidationResult(null);
    setLoading(false);
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
            <div className="space-y-3">
              <WorkflowDetail
                workflow={workflow}
              />
              <Button asChild variant="primary" className="gap-1.5">
                <Link to={`/workflows/${workflow.metadata.workflow_id}/new-run`}>
                  <FilePenLine aria-hidden="true" size={16} />
                  Author inputs
                </Link>
              </Button>
              <DeveloperSchemaDetails
                workflowId={workflow.metadata.workflow_id}
                schemaHints={schemaHints}
              />
            </div>
          </Panel>
        )}
        {workflow && supportsServerPathSamples && (
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
        {workflow && supportsServerPathSamples && (
          <Panel title="Run progress">
            <RunProgressPanel
              workflowId={workflow.metadata.workflow_id}
              validationResult={validationResult}
              runClient={runClient}
              executionAvailability={workflow.availability.execution}
              onRunCreated={handleRunCreated}
            />
          </Panel>
        )}
        {supportsServerPathSamples && (
          <Panel title="Validation results">
            <IssuePanel issues={displayedIssues} onAskAgent={handleAskAgent} />
          </Panel>
        )}
      </section>

      {supportsServerPathSamples && (
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
      )}
    </>
  );
}
