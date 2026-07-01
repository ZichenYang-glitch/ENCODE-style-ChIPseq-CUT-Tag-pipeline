import { useEffect, useMemo, useState } from 'react';
import { createStubWorkflowClient } from './api/client';
import type {
  Issue,
  ValidateWorkflowResponse,
  WorkflowInputs,
  WorkflowSchema,
  WorkflowSummary,
} from './api/types';
import { Panel } from './components/Panel';
import { AgentSidebar } from './features/agent-sidebar/AgentSidebar';
import { IssuePanel } from './features/issues-panel/IssuePanel';
import { ValidationWorkspace } from './features/validation-workspace/ValidationWorkspace';
import { WorkflowCatalog } from './features/workflow-catalog/WorkflowCatalog';
import { WorkflowDetail } from './features/workflow-detail/WorkflowDetail';

const client = createStubWorkflowClient();

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

export function App() {
  const [workflows, setWorkflows] = useState<WorkflowSummary[]>([]);
  const [selectedWorkflowId, setSelectedWorkflowId] = useState<string | null>(
    null,
  );
  const [schemaHints, setSchemaHints] = useState<WorkflowSchema | null>(null);
  const [configText, setConfigText] = useState('{}');
  const [samplesText, setSamplesText] = useState('');
  const [optionsText, setOptionsText] = useState('{}');
  const [validationResult, setValidationResult] =
    useState<ValidateWorkflowResponse | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    let cancelled = false;
    client.listWorkflows().then((response) => {
      if (!cancelled) {
        setWorkflows(response.workflows);
      }
    });
    return () => {
      cancelled = true;
    };
  }, []);

  useEffect(() => {
    if (!selectedWorkflowId) {
      setSchemaHints(null);
      return;
    }
    let cancelled = false;
    setSchemaHints(null);
    client.getWorkflowSchema(selectedWorkflowId).then((response) => {
      if (!cancelled) {
        setSchemaHints(response.schema_hints);
      }
    });
    return () => {
      cancelled = true;
    };
  }, [selectedWorkflowId]);

  const selectedWorkflow = useMemo(
    () => workflows.find((w) => w.metadata.workflow_id === selectedWorkflowId),
    [workflows, selectedWorkflowId],
  );

  const displayedIssues = useMemo(() => {
    return validationResult?.issues ?? [];
  }, [validationResult]);

  async function handleValidate() {
    if (!selectedWorkflowId) return;

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
        workflow_id: selectedWorkflowId,
        value: null,
        issues,
      });
      setLoading(false);
      return;
    }

    const inputs: WorkflowInputs = {
      config,
      samples: samplesText || null,
      options,
    };

    const response = await client.validateWorkflow(selectedWorkflowId, inputs);
    setValidationResult(response);
    setLoading(false);
  }

  return (
    <div className="flex min-h-screen flex-col bg-[var(--color-bg)] text-[var(--color-text)]">
      <header className="border-b border-[var(--color-border)] bg-[var(--color-surface)] px-4 py-3">
        <h1 className="text-base font-semibold tracking-wide text-[var(--color-accent)]">Workflow Platform</h1>
      </header>
      <main className="mx-auto flex w-full max-w-screen-2xl flex-1 flex-col gap-3 p-3 lg:flex-row">
        <aside className="w-full shrink-0 lg:w-56">
          <Panel title="Workflows">
            {workflows.length === 0 ? (
              <p className="text-sm text-[var(--color-text-muted)]">Loading…</p>
            ) : (
              <WorkflowCatalog
                workflows={workflows}
                selectedWorkflowId={selectedWorkflowId}
                onSelect={setSelectedWorkflowId}
              />
            )}
          </Panel>
        </aside>

        <section className="flex min-w-0 flex-1 flex-col gap-3">
          {selectedWorkflow && (
            <Panel title="Workflow detail">
              <WorkflowDetail
                workflowId={selectedWorkflow.metadata.workflow_id}
                workflowName={selectedWorkflow.metadata.name}
                schemaHints={schemaHints}
              />
            </Panel>
          )}
          {selectedWorkflow && (
            <Panel title="Validation workspace">
              <ValidationWorkspace
                workflowId={selectedWorkflow.metadata.workflow_id}
                configText={configText}
                samplesText={samplesText}
                optionsText={optionsText}
                loading={loading}
                onConfigChange={setConfigText}
                onSamplesChange={setSamplesText}
                onOptionsChange={setOptionsText}
                onValidate={handleValidate}
              />
            </Panel>
          )}
          {!selectedWorkflow && (
            <Panel title="Validation workspace">
              <p className="text-sm text-[var(--color-text-muted)]">
                Select a workflow from the catalog to begin validation.
              </p>
            </Panel>
          )}
          <Panel title="Validation results">
            <IssuePanel issues={displayedIssues} />
          </Panel>
        </section>

        <aside className="flex w-full flex-col gap-3 lg:w-72">
          <AgentSidebar
            workflowId={selectedWorkflowId}
            issues={displayedIssues}
          />
        </aside>
      </main>
    </div>
  );
}
