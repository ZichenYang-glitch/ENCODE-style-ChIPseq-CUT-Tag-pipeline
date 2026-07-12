import { useCallback } from 'react';
import {
  useLocation,
  useNavigate,
  useParams,
  useSearchParams,
} from 'react-router-dom';
import { useClients } from '../../api/client-context';
import { Panel } from '../../components/Panel';
import { RunProgressPanel } from '../../features/run-progress/RunProgressPanel';

export function RunDetailPage() {
  const { runId } = useParams<{ runId: string }>();
  const location = useLocation();
  const navigate = useNavigate();
  const [searchParams, setSearchParams] = useSearchParams();
  const { runClient } = useClients();
  const selectedArtifactId = searchParams.get('artifact');
  const activeView =
    searchParams.get('view') === 'artifacts' || selectedArtifactId !== null
      ? 'artifacts'
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
    (view: 'activity' | 'artifacts') => {
      const next = new URLSearchParams(searchParams);
      if (view === 'artifacts') {
        next.set('view', 'artifacts');
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
