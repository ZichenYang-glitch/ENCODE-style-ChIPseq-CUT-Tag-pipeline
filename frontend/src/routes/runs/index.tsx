import { useInfiniteQuery, useQuery, useQueryClient } from '@tanstack/react-query';
import { History, RefreshCw, RotateCcw, Workflow } from 'lucide-react';
import type { ReactNode } from 'react';
import { Link, useSearchParams } from 'react-router-dom';
import { listRuns } from '../../api/generated/runs/runs';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { useClients } from '../../api/client-context';
import { RunHistoryList } from '../../features/run-history/RunHistoryList';
import {
  flattenRunHistoryPages,
  isRunHistoryCursorError,
  parseRunHistoryFilters,
  RUN_HISTORY_PAGE_SIZE,
  RUN_HISTORY_STATUSES,
  runHistoryPaginationAnomaly,
  runHistoryQueryKey,
  safeRunHistoryNextCursor,
  validateRunHistoryResponse,
  type RunHistoryFilters,
} from '../../features/run-history/runHistoryState';

export function RunHistoryPage() {
  const [searchParams, setSearchParams] = useSearchParams();
  const parsed = parseRunHistoryFilters(searchParams);
  const filters: RunHistoryFilters = parsed.ok
    ? parsed.filters
    : { workflowId: null, status: null };
  const queryClient = useQueryClient();
  const { workflowClient } = useClients();
  const workflowsQuery = useQuery({
    queryKey: ['workflows'],
    queryFn: () => workflowClient.listWorkflows(),
  });
  const queryKey = runHistoryQueryKey(filters);
  const historyQuery = useInfiniteQuery({
    queryKey,
    enabled: parsed.ok,
    initialPageParam: undefined as string | undefined,
    queryFn: async ({ pageParam }) =>
      validateRunHistoryResponse(
        await listRuns({
          after: pageParam,
          limit: RUN_HISTORY_PAGE_SIZE,
          workflow_id: filters.workflowId ?? undefined,
          status: filters.status ?? undefined,
        }),
        pageParam,
      ),
    getNextPageParam: (lastPage, _allPages, _lastPageParam, allPageParams) =>
      safeRunHistoryNextCursor(lastPage, allPageParams),
    retry: 1,
    retryDelay: 100,
  });

  const pages = historyQuery.data?.pages ?? [];
  const pageParams = historyQuery.data?.pageParams ?? [];
  const runs = flattenRunHistoryPages(pages);
  const protocolRecovery =
    runHistoryPaginationAnomaly(pages, pageParams) ||
    isRunHistoryCursorError(historyQuery.error);
  const hasCachedRuns = runs.length > 0;
  const hasError = historyQuery.error !== null;

  function setFilter(name: 'workflow_id' | 'status', value: string) {
    const next = new URLSearchParams(searchParams);
    if (value) next.set(name, value);
    else next.delete(name);
    setSearchParams(next);
  }

  function resetFilters() {
    const next = new URLSearchParams(searchParams);
    next.delete('workflow_id');
    next.delete('status');
    setSearchParams(next, { replace: true });
  }

  function reloadFromFirstPage() {
    void queryClient.resetQueries({ queryKey, exact: true });
  }

  const workflows = workflowsQuery.data?.ok
    ? workflowsQuery.data.workflows
    : [];
  const knownWorkflowIds = new Set(
    workflows.map((workflow) => workflow.metadata.workflow_id),
  );

  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Run history">
        <div className="min-w-0 space-y-3" data-testid="run-history-page">
          <div className="flex min-w-0 flex-wrap items-start justify-between gap-3">
            <div className="min-w-0">
              <h2 className="flex items-center gap-2 text-sm font-semibold">
                <History aria-hidden="true" size={17} />
                Canonical runs
              </h2>
              <p className="mt-1 text-xs text-[var(--color-text-muted)]">
                Persistent lifecycle summaries from the local platform database.
              </p>
            </div>
            <Button
              variant="secondary"
              onClick={() => void historyQuery.refetch()}
              disabled={!parsed.ok || historyQuery.isFetching}
              aria-label="Refresh run history"
              title="Refresh run history"
            >
              <RefreshCw
                className={`mr-1.5 h-4 w-4 ${historyQuery.isFetching ? 'animate-spin' : ''}`}
                aria-hidden="true"
              />
              Refresh
            </Button>
          </div>

          <div className="grid min-w-0 gap-2 border-y border-[var(--color-border)] py-3 sm:grid-cols-2">
            <label className="min-w-0 text-xs font-medium">
              Workflow
              <select
                className="mt-1 block h-9 w-full min-w-0 rounded border border-[var(--color-border)] bg-[var(--color-surface)] px-2 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
                aria-label="Filter runs by workflow"
                value={parsed.ok ? filters.workflowId ?? '' : ''}
                disabled={!parsed.ok}
                onChange={(event) => setFilter('workflow_id', event.target.value)}
              >
                <option value="">All workflows</option>
                {filters.workflowId && !knownWorkflowIds.has(filters.workflowId) && (
                  <option value={filters.workflowId}>{filters.workflowId}</option>
                )}
                {workflows.map((workflow) => (
                  <option
                    key={workflow.metadata.workflow_id}
                    value={workflow.metadata.workflow_id}
                  >
                    {workflow.metadata.name}
                  </option>
                ))}
              </select>
            </label>
            <label className="min-w-0 text-xs font-medium">
              Status
              <select
                className="mt-1 block h-9 w-full min-w-0 rounded border border-[var(--color-border)] bg-[var(--color-surface)] px-2 text-sm capitalize focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
                aria-label="Filter runs by status"
                value={parsed.ok ? filters.status ?? '' : ''}
                disabled={!parsed.ok}
                onChange={(event) => setFilter('status', event.target.value)}
              >
                <option value="">All statuses</option>
                {RUN_HISTORY_STATUSES.map((status) => (
                  <option key={status} value={status}>{status}</option>
                ))}
              </select>
            </label>
          </div>

          {!parsed.ok ? (
            <StateMessage title="Run filters could not be used" tone="error">
              The URL contains an unsupported or repeated filter. Reset it before loading run history.
              <Button className="mt-3" variant="secondary" onClick={resetFilters}>
                <RotateCcw className="mr-1.5 h-4 w-4" aria-hidden="true" />
                Reset filters
              </Button>
            </StateMessage>
          ) : historyQuery.isPending ? (
            <RunHistorySkeleton />
          ) : hasError && !hasCachedRuns ? (
            <StateMessage title="Run history could not be loaded" tone="error">
              {protocolRecovery
                ? 'The pagination state could not be confirmed. Reload from the first page.'
                : 'The canonical run history is temporarily unavailable.'}
              <Button
                className="mt-3"
                variant="secondary"
                onClick={protocolRecovery ? reloadFromFirstPage : () => void historyQuery.refetch()}
              >
                {protocolRecovery ? 'Reload from first page' : 'Retry run history'}
              </Button>
            </StateMessage>
          ) : runs.length === 0 ? (
            <StateMessage title={filters.workflowId || filters.status ? 'No runs match these filters' : 'No runs yet'}>
              {filters.workflowId || filters.status ? (
                <Button className="mt-3" variant="secondary" onClick={resetFilters}>
                  Clear filters
                </Button>
              ) : (
                <Button asChild className="mt-3" variant="secondary">
                  <Link to="/workflows">
                    <Workflow className="mr-1.5 h-4 w-4" aria-hidden="true" />
                    Browse workflows
                  </Link>
                </Button>
              )}
            </StateMessage>
          ) : (
            <>
              {(hasError || protocolRecovery) && (
                <div className="flex min-w-0 flex-wrap items-center justify-between gap-2 border border-[var(--color-error)] bg-[var(--color-error-bg)] p-2 text-xs" role="alert">
                  <span>
                    {protocolRecovery
                      ? 'Pagination could not be confirmed. Loaded rows are retained.'
                      : 'Refresh failed. Previously loaded canonical rows are retained.'}
                  </span>
                  <Button
                    className="px-2 py-1 text-xs"
                    variant="secondary"
                    onClick={protocolRecovery ? reloadFromFirstPage : () => void historyQuery.refetch()}
                  >
                    {protocolRecovery ? 'Reload from first page' : 'Retry'}
                  </Button>
                </div>
              )}
              <p className="text-xs text-[var(--color-text-muted)]" role="status">
                {runs.length} run{runs.length === 1 ? '' : 's'} loaded
              </p>
              <RunHistoryList
                runs={runs}
                hasNextPage={Boolean(historyQuery.hasNextPage) && !protocolRecovery}
                isFetchingNextPage={historyQuery.isFetchingNextPage}
                onLoadMore={() => void historyQuery.fetchNextPage()}
              />
            </>
          )}
        </div>
      </Panel>
    </section>
  );
}

function StateMessage({
  title,
  tone = 'neutral',
  children,
}: {
  title: string;
  tone?: 'neutral' | 'error';
  children: ReactNode;
}) {
  return (
    <div
      className={`min-h-36 border-y border-[var(--color-border)] py-6 text-sm ${tone === 'error' ? 'text-[var(--color-error)]' : ''}`}
      role={tone === 'error' ? 'alert' : 'status'}
    >
      <h3 className="font-semibold">{title}</h3>
      <div className="mt-1 text-[var(--color-text-muted)]">{children}</div>
    </div>
  );
}

function RunHistorySkeleton() {
  return (
    <div className="min-h-64 animate-pulse space-y-2" aria-label="Loading run history">
      {[0, 1, 2, 3].map((row) => (
        <div key={row} className="h-12 rounded bg-slate-100" />
      ))}
    </div>
  );
}
