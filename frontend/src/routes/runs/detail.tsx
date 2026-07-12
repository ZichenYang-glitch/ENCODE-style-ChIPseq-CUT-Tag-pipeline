import { useCallback } from 'react';
import { useLocation, useNavigate, useParams } from 'react-router-dom';
import { useClients } from '../../api/client-context';
import { Panel } from '../../components/Panel';
import { RunProgressPanel } from '../../features/run-progress/RunProgressPanel';

export function RunDetailPage() {
  const { runId } = useParams<{ runId: string }>();
  const location = useLocation();
  const navigate = useNavigate();
  const { runClient } = useClients();
  const beginPreflight =
    typeof location.state === 'object' &&
    location.state !== null &&
    'beginPreflight' in location.state &&
    location.state.beginPreflight === true;
  const preflightRequestId =
    typeof location.state === 'object' &&
    location.state !== null &&
    'preflightRequestId' in location.state &&
    typeof location.state.preflightRequestId === 'string'
      ? location.state.preflightRequestId
      : null;
  const handlePreflightConsumed = useCallback(() => {
    navigate(location.pathname, { replace: true, state: null });
  }, [location.pathname, navigate]);

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
        <RunProgressPanel
          runId={runId}
          runClient={runClient}
          beginPreflight={beginPreflight}
          preflightRequestId={preflightRequestId}
          onPreflightConsumed={handlePreflightConsumed}
        />
      </Panel>
    </section>
  );
}
