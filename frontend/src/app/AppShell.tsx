import {
  Dna,
  FileArchive,
  FilePenLine,
  History,
  ListTree,
  LogOut,
} from 'lucide-react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { Link, matchPath, Outlet, useLocation, useNavigate } from 'react-router-dom';
import {
  getTerminalEmailPreference,
  setTerminalEmailPreference,
} from '../api/generated/auth/auth';
import { Button } from '../components/Button';
import { useAuth } from './auth';

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
      <header className="border-b border-[var(--color-border)] bg-[var(--color-surface)] px-4 py-3">
        <div className="mx-auto flex w-full max-w-screen-2xl flex-wrap items-center justify-between gap-x-4 gap-y-2">
          <h1 className="text-base font-semibold tracking-wide text-[var(--color-accent)]">
            <Link
              className="inline-flex items-center gap-2 rounded focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)] focus:ring-offset-2"
              to="/workflows"
            >
              <Dna aria-hidden="true" size={18} />
              HelixWeave
            </Link>
          </h1>
          <nav
            aria-label="Primary"
            className="flex min-w-0 flex-wrap items-center gap-2"
          >
            <Button
              asChild
              className="gap-1.5"
              variant={workflowsCurrent ? 'primary' : 'secondary'}
            >
              <Link
                aria-current={workflowsCurrent ? 'page' : undefined}
                to="/workflows"
              >
                <ListTree aria-hidden="true" size={16} />
                Workflows
              </Link>
            </Button>
            <Button
              asChild
              className="gap-1.5"
              variant={runsCurrent ? 'primary' : 'secondary'}
            >
              <Link
                aria-current={runsCurrent ? 'page' : undefined}
                to="/runs"
              >
                <History aria-hidden="true" size={16} />
                Runs
              </Link>
            </Button>
            <Button
              asChild
              className="gap-1.5"
              variant={artifactsCurrent ? 'primary' : 'secondary'}
            >
              <Link
                aria-current={artifactsCurrent ? 'page' : undefined}
                to="/artifacts"
              >
                <FileArchive aria-hidden="true" size={16} />
                Artifacts
              </Link>
            </Button>
            {workflowId && (
              <Button
                asChild
                className="gap-1.5"
                variant={authoringRoute ? 'primary' : 'secondary'}
              >
                <Link
                  aria-current={authoringRoute ? 'page' : undefined}
                  to={`/workflows/${encodeURIComponent(workflowId)}/new-run`}
                >
                  <FilePenLine aria-hidden="true" size={16} />
                  New analysis
                </Link>
              </Button>
            )}
            {principal !== null && (
              <span className="ml-2 flex items-center gap-2 text-sm text-[var(--color-text-muted)]">
                <span data-testid="current-username">{principal.username}</span>
                {isMember && terminalEmailPreference.data && (
                  <label className="inline-flex items-center gap-1.5">
                    <input
                      aria-label="Email me when my runs finish"
                      checked={
                        terminalEmailPreference.data.terminal_email_enabled
                      }
                      disabled={
                        !terminalEmailPreference.data.address_configured ||
                        updateTerminalEmailPreference.isPending
                      }
                      type="checkbox"
                      onChange={(event) => {
                        updateTerminalEmailPreference.mutate(
                          event.currentTarget.checked,
                        );
                      }}
                    />
                    Terminal emails
                    {!terminalEmailPreference.data.address_configured && (
                      <span title="Ask the local operator to configure your notification address.">
                        (address not configured)
                      </span>
                    )}
                  </label>
                )}
                {principal.role === 'administrator' && (
                  <span title="Administrator terminal emails are managed by the local operator.">
                    Notifications: operator-managed
                  </span>
                )}
                <Button
                  className="gap-1.5"
                  variant="secondary"
                  onClick={() => {
                    void logout().then(() => navigate('/login'));
                  }}
                >
                  <LogOut aria-hidden="true" size={16} />
                  Sign out
                </Button>
              </span>
            )}
          </nav>
        </div>
      </header>
      <main className="mx-auto flex w-full max-w-screen-2xl flex-1 flex-col gap-3 p-3 lg:flex-row">
        <Outlet />
      </main>
    </div>
  );
}
