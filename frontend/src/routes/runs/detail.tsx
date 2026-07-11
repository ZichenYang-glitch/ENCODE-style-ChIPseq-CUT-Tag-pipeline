import { useParams } from 'react-router-dom';
import { useClients } from '../../api/client-context';
import { Panel } from '../../components/Panel';
import { RunProgressPanel } from '../../features/run-progress/RunProgressPanel';

export function RunDetailPage() {
  const { runId } = useParams<{ runId: string }>();
  const { runClient } = useClients();

  if (!runId) {
    return (
      <Panel title="Run">
        <p className="text-sm text-[var(--color-error)]">Missing run ID.</p>
      </Panel>
    );
  }

  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Run progress">
        <RunProgressPanel runId={runId} runClient={runClient} />
      </Panel>
    </section>
  );
}
