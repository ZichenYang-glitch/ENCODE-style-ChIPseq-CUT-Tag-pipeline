import { useEffect, useState } from 'react';
import { createStubWorkflowClient } from './api/client';
import type { WorkflowSummary } from './api/types';

const client = createStubWorkflowClient();

export function App() {
  const [workflows, setWorkflows] = useState<WorkflowSummary[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    client.listWorkflows().then((response) => {
      if (!cancelled && response.ok) {
        setWorkflows(response.workflows);
      }
      setLoading(false);
    });
    return () => {
      cancelled = true;
    };
  }, []);

  return (
    <div className="min-h-screen bg-[var(--color-bg)] text-[var(--color-text)]">
      <header className="border-b border-[var(--color-border)] px-6 py-4">
        <h1 className="text-xl font-semibold">Workflow Platform</h1>
      </header>
      <main className="p-6">
        {loading ? (
          <p>Loading workflows…</p>
        ) : (
          <section>
            <h2 className="mb-2 text-lg font-medium">Available workflows</h2>
            <ul className="space-y-2">
              {workflows.map((workflow) => (
                <li
                  key={workflow.metadata.workflow_id}
                  className="rounded border border-[var(--color-border)] p-3"
                >
                  <div className="font-medium">{workflow.metadata.name}</div>
                  <div className="text-sm text-[var(--color-text-muted)]">
                    {workflow.metadata.workflow_id}
                  </div>
                </li>
              ))}
            </ul>
          </section>
        )}
      </main>
    </div>
  );
}
