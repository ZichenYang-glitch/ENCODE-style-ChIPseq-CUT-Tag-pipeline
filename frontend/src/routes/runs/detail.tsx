import { useCallback } from 'react';
import {
  useLocation,
  useNavigate,
  useParams,
  useSearchParams,
} from 'react-router-dom';
import { useClients } from '../../api/client-context';
import { Panel } from '../../components/Panel';
import {
  RunProgressPanel,
  type RunDetailView,
} from '../../features/run-progress/RunProgressPanel';

export function RunDetailPage() {
  const { runId } = useParams<{ runId: string }>();
  const location = useLocation();
  const navigate = useNavigate();
  const [searchParams, setSearchParams] = useSearchParams();
  const { runClient, workflowClient } = useClients();
  const selectedArtifactId = searchParams.get('artifact');
  const requestedView = searchParams.get('view');
  const activeView: RunDetailView =
    requestedView === 'artifacts' || selectedArtifactId !== null
      ? 'artifacts'
      : requestedView === 'qc'
        ? 'qc'
        : 'activity';
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
    navigate(`${location.pathname}${location.search}`, {
      replace: true,
      state: null,
    });
  }, [location.pathname, location.search, navigate]);

  const handleViewChange = useCallback(
    (view: RunDetailView) => {
      const next = new URLSearchParams(searchParams);
      if (view === 'artifacts') {
        next.set('view', 'artifacts');
      } else if (view === 'qc') {
        next.set('view', 'qc');
        next.delete('artifact');
      } else {
        next.delete('view');
        next.delete('artifact');
      }
      setSearchParams(next);
    },
    [searchParams, setSearchParams],
  );

  const handleArtifactSelect = useCallback(
    (artifactId: string) => {
      const next = new URLSearchParams(searchParams);
      next.set('view', 'artifacts');
      next.set('artifact', artifactId);
      setSearchParams(next);
    },
    [searchParams, setSearchParams],
  );

  if (!runId) {
    return (
      <Panel title="Run">
        <p className="text-sm text-[var(--color-error)]">Missing run ID.</p>
      </Panel>
    );
  }

  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Run">
        <RunProgressPanel
          runId={runId}
          runClient={runClient}
          workflowClient={workflowClient}
          beginPreflight={beginPreflight}
          preflightRequestId={preflightRequestId}
          onPreflightConsumed={handlePreflightConsumed}
          activeView={activeView}
          selectedArtifactId={selectedArtifactId}
          onViewChange={handleViewChange}
          onArtifactSelect={handleArtifactSelect}
        />
      </Panel>
    </section>
  );
}
