import { Panel } from '../../components/Panel';

export function WorkflowsIndexPage() {
  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Workflow detail">
        <p className="text-sm text-[var(--color-text-muted)]">
          Select a workflow from the catalog to view details.
        </p>
      </Panel>
      <Panel title="Validation workspace">
        <p className="text-sm text-[var(--color-text-muted)]">
          Validation requires a selected workflow.
        </p>
      </Panel>
      <Panel title="Run progress">
        <p className="text-sm text-[var(--color-text-muted)]">
          Create a run after selecting and validating a workflow.
        </p>
      </Panel>
    </section>
  );
}
