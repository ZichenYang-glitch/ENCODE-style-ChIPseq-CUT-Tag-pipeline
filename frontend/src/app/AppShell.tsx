import {
  Dna,
  FileArchive,
  FilePenLine,
  History,
  ListTree,
} from 'lucide-react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { Link, matchPath, Outlet, useLocation, useNavigate } from 'react-router-dom';
import {
  getTerminalEmailPreference,
  setTerminalEmailPreference,
} from '../api/generated/auth/auth';
import { useAuth } from './auth';
import { AccountMenu } from './AccountMenu';

export function AppShell() {
  const { pathname } = useLocation();
  const navigate = useNavigate();
  const { principal, logout } = useAuth();
  const queryClient = useQueryClient();
  const isMember = principal?.role === 'member';
  const terminalEmailPreference = useQuery({
    queryKey: ['auth', 'terminal-email-preference'],
    queryFn: () => getTerminalEmailPreference(),
    enabled: isMember,
  });
  const updateTerminalEmailPreference = useMutation({
    mutationFn: (enabled: boolean) =>
      setTerminalEmailPreference({ terminal_email_enabled: enabled }),
    onSuccess: (preference) => {
      queryClient.setQueryData(
        ['auth', 'terminal-email-preference'],
        preference,
      );
    },
  });
  const workflowRoute = matchPath(
    { path: '/workflows/:workflowId/*', end: false },
    pathname,
  );
  const workflowDetailRoute = matchPath(
    { path: '/workflows/:workflowId', end: true },
    pathname,
  );
  const workflowCatalogRoute = matchPath(
    { path: '/workflows', end: true },
    pathname,
  );
  const authoringRoute = matchPath(
    { path: '/workflows/:workflowId/new-run', end: true },
    pathname,
  );
  const runHistoryRoute = matchPath({ path: '/runs', end: true }, pathname);
  const runDetailRoute = matchPath(
    { path: '/runs/:runId', end: true },
    pathname,
  );
  const artifactIndexRoute = matchPath(
    { path: '/artifacts', end: true },
    pathname,
  );
  const artifactDetailRoute = matchPath(
    { path: '/artifacts/:runId/:artifactId', end: true },
    pathname,
  );
  const workflowId = workflowRoute?.params.workflowId;
  const workflowsCurrent =
    workflowCatalogRoute !== null || workflowDetailRoute !== null;
  const runsCurrent = runHistoryRoute !== null || runDetailRoute !== null;
  const artifactsCurrent =
    artifactIndexRoute !== null || artifactDetailRoute !== null;

  return (
    <div className="flex min-h-screen flex-col bg-[var(--color-bg)] text-[var(--color-text)]">
      <header className="border-b border-[var(--color-border)] bg-[var(--color-surface)] px-3 py-2 sm:px-4">
        <div className="mx-auto grid w-full max-w-screen-2xl grid-cols-[minmax(0,1fr)_auto] items-center gap-x-3 md:grid-cols-[auto_minmax(0,1fr)_auto] md:gap-x-6">
          <h1 className="min-w-0 text-base font-semibold text-[var(--color-accent)]">
            <Link
              className="inline-flex min-h-11 items-center gap-2 rounded-[4px] sm:min-h-9"
              to="/workflows"
            >
              <Dna aria-hidden="true" size={18} />
              HelixWeave
            </Link>
          </h1>
          <nav
            aria-label="Primary"
            className="col-span-2 row-start-2 grid min-w-0 grid-flow-col auto-cols-fr items-center gap-1 md:col-span-1 md:col-start-2 md:row-start-1 md:flex"
          >
            <Link
              className="primary-nav-link"
              aria-current={workflowsCurrent ? 'page' : undefined}
              to="/workflows"
            >
              <ListTree className="hidden sm:block" aria-hidden="true" size={16} />
              Workflows
            </Link>
            <Link
              className="primary-nav-link"
              aria-current={runsCurrent ? 'page' : undefined}
              to="/runs"
            >
              <History className="hidden sm:block" aria-hidden="true" size={16} />
              Runs
            </Link>
            <Link
              className="primary-nav-link"
              aria-current={artifactsCurrent ? 'page' : undefined}
              to="/artifacts"
            >
              <FileArchive className="hidden sm:block" aria-hidden="true" size={16} />
              Artifacts
            </Link>
            {workflowId && (
              <Link
                className="primary-nav-link"
                aria-label="New analysis"
                title="New analysis"
                aria-current={authoringRoute ? 'page' : undefined}
                to={`/workflows/${encodeURIComponent(workflowId)}/new-run`}
              >
                <FilePenLine className="hidden sm:block" aria-hidden="true" size={16} />
                <span className="sm:hidden">New</span>
                <span className="hidden sm:inline">New analysis</span>
              </Link>
            )}
          </nav>
          {principal !== null && (
            <div className="col-start-2 row-start-1 justify-self-end md:col-start-3">
              <AccountMenu
                principal={principal}
                terminalEmailPreference={terminalEmailPreference.data}
                preferencePending={updateTerminalEmailPreference.isPending}
                onPreferenceChange={(enabled) =>
                  updateTerminalEmailPreference.mutate(enabled)
                }
                onSignOut={() => {
                  void logout().then(() => navigate('/login'));
                }}
              />
            </div>
          )}
        </div>
      </header>
      <main className="mx-auto flex w-full max-w-screen-2xl flex-1 flex-col gap-3 p-3 lg:flex-row lg:p-4">
        <Outlet />
      </main>
    </div>
  );
}
