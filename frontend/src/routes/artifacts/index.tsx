import { useInfiniteQuery, useQueryClient } from '@tanstack/react-query';
import { FileArchive, RefreshCw, RotateCcw } from 'lucide-react';
import { useState, type FormEvent, type ReactNode } from 'react';
import { Link, useParams, useSearchParams } from 'react-router-dom';
import { listArtifactPublications } from '../../api/generated/artifact-publications/artifact-publications';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { ArtifactPublicationDetail } from '../../features/artifact-publications/ArtifactPublicationDetail';
import { ArtifactPublicationList } from '../../features/artifact-publications/ArtifactPublicationList';
import {
  ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE,
  artifactPublicationPaginationAnomaly,
  artifactPublicationQueryKey,
  artifactPublicationSearchParams,
  flattenArtifactPublicationPages,
  hasActiveArtifactPublicationFilters,
  isArtifactPublicationRecoveryError,
  parseArtifactPublicationDetailIdentity,
  parseArtifactPublicationFilters,
  safeArtifactPublicationNextCursor,
  validateArtifactPublicationListResponse,
  type ArtifactPublicationFilters,
} from '../../features/artifact-publications/artifactPublicationState';

const EMPTY_FILTERS: ArtifactPublicationFilters = {
  projectId: null,
  runId: null,
  workflowId: null,
  outputType: null,
  associatedRunSampleRevisionId: null,
  publishedFrom: null,
  publishedBefore: null,
  currentOnly: true,
  after: null,
  limit: ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE,
};

export function ArtifactPublicationsPage() {
  const [searchParams, setSearchParams] = useSearchParams();
  const parsed = parseArtifactPublicationFilters(searchParams);
  const filters = parsed.ok ? parsed.filters : EMPTY_FILTERS;
  const queryClient = useQueryClient();
  const queryKey = artifactPublicationQueryKey(filters);
  const publicationsQuery = useInfiniteQuery({
    queryKey,
    enabled: parsed.ok,
    initialPageParam: filters.after ?? (undefined as string | undefined),
    queryFn: async ({ pageParam }) =>
      validateArtifactPublicationListResponse(
        await listArtifactPublications({
          project_id: filters.projectId ?? undefined,
          run_id: filters.runId ?? undefined,
          workflow_id: filters.workflowId ?? undefined,
          output_type: filters.outputType ?? undefined,
          associated_run_sample_revision_id:
            filters.associatedRunSampleRevisionId ?? undefined,
          published_from: filters.publishedFrom ?? undefined,
          published_before: filters.publishedBefore ?? undefined,
          current_only: filters.currentOnly,
          after: pageParam,
          limit: filters.limit,
        }),
        pageParam,
      ),
    getNextPageParam: (lastPage, allPages, _lastPageParam, allPageParams) =>
      safeArtifactPublicationNextCursor(lastPage, allPages, allPageParams),
    retry: (failureCount, error) =>
      !(error instanceof ApiError && error.status >= 400 && error.status < 500) &&
      !isArtifactPublicationRecoveryError(error) &&
      failureCount < 1,
    retryDelay: 100,
  });

  const pages = publicationsQuery.data?.pages ?? [];
  const pageParams = publicationsQuery.data?.pageParams ?? [];
  const publications = flattenArtifactPublicationPages(pages);
  const protocolRecovery =
    artifactPublicationPaginationAnomaly(pages, pageParams) ||
    isArtifactPublicationRecoveryError(publicationsQuery.error);
  const hasCachedPublications = publications.length > 0;
  const hasError = publicationsQuery.error !== null;
  const hasFilters = hasActiveArtifactPublicationFilters(filters);
  const activeFilters = parsed.ok ? activeFilterEntries(filters) : [];
  const [filtersOpen, setFiltersOpen] = useState(!parsed.ok || hasFilters);

  function resetFilters() {
    setFiltersOpen(false);
    setSearchParams(new URLSearchParams(), { replace: true });
  }

  function reloadFromFirstPage() {
    if (filters.after !== null) {
      const next = artifactPublicationSearchParams({ ...filters, after: null });
      setSearchParams(next, { replace: true });
      return;
    }
    void queryClient.resetQueries({ queryKey, exact: true });
  }

  function applyFilters(nextFilters: ArtifactPublicationFilters) {
    setFiltersOpen(hasActiveArtifactPublicationFilters(nextFilters));
    setSearchParams(artifactPublicationSearchParams({ ...nextFilters, after: null }));
  }

  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Artifacts">
        <div className="min-w-0 space-y-3" data-testid="artifact-publications-page">
          <div className="flex min-w-0 flex-wrap items-start justify-between gap-3">
            <div className="min-w-0">
              <h2 className="flex items-center gap-2 text-sm font-semibold">
                <FileArchive aria-hidden="true" size={17} />
                Artifact publications
              </h2>
              <p className="mt-1 text-xs text-[var(--color-text-muted)]">
                Cross-run publication metadata from the canonical local index.
              </p>
            </div>
            <Button
              variant="secondary"
              onClick={() => void publicationsQuery.refetch()}
              disabled={!parsed.ok || publicationsQuery.isFetching}
              aria-label="Refresh artifact publications"
            >
              <RefreshCw
                className={`mr-1.5 h-4 w-4 ${publicationsQuery.isFetching ? 'animate-spin' : ''}`}
                aria-hidden="true"
              />
              Refresh
            </Button>
          </div>

          {parsed.ok && (
            <div className="border-y border-[var(--color-border)]">
              {activeFilters.length > 0 && (
                <ul className="flex min-w-0 flex-wrap gap-1.5 border-b border-[var(--color-border)] py-2" aria-label="Active artifact filters">
                  {activeFilters.map(([label, value]) => (
                    <li
                      key={label}
                      className="max-w-full truncate rounded-[4px] border border-[var(--color-info-border)] bg-[var(--color-info-bg)] px-2 py-1 text-xs text-[var(--color-info)]"
                      title={`${label}: ${value}`}
                    >
                      {label}: {value}
                    </li>
                  ))}
                </ul>
              )}
              <details
                key={searchParams.toString()}
                open={filtersOpen}
                onToggle={(event) => setFiltersOpen(event.currentTarget.open)}
                data-testid="artifact-filter-disclosure"
              >
                <summary className="flex min-h-11 cursor-pointer items-center justify-between gap-3 py-2 text-sm font-semibold">
                  <span>Filters ({activeFilters.length} active)</span>
                  <span className="text-xs font-normal text-[var(--color-text-muted)]">
                    Apply and Reset keep URL filter semantics.
                  </span>
                </summary>
                <ArtifactPublicationFilterForm
                  filters={filters}
                  onApply={applyFilters}
                  onReset={resetFilters}
                />
              </details>
            </div>
          )}

          {!parsed.ok && (
            <details open className="border-y border-[var(--color-error-border)]" data-testid="artifact-filter-disclosure">
              <summary className="min-h-11 cursor-pointer py-2 text-sm font-semibold text-[var(--color-error)]">
                Filters (invalid URL)
              </summary>
              <p className="pb-2 text-xs text-[var(--color-text-muted)]">
                The invalid URL filters must be reset before the form can be used.
              </p>
            </details>
          )}

          {!parsed.ok ? (
            <StateMessage title="Artifact filters could not be used" tone="error">
              The URL contains an unsupported, repeated, or unsafe filter. Reset it before loading artifact publications.
              <Button className="mt-3" variant="secondary" onClick={resetFilters}>
                <RotateCcw className="mr-1.5 h-4 w-4" aria-hidden="true" />
                Reset filters
              </Button>
            </StateMessage>
          ) : publicationsQuery.isPending ? (
            <PublicationListSkeleton />
          ) : hasError && !hasCachedPublications ? (
            <StateMessage title="Artifact publications could not be loaded" tone="error">
              {protocolRecovery
                ? 'The publication data or pagination state could not be confirmed. Reload from the first page.'
                : 'The artifact publication index is temporarily unavailable.'}
              <Button
                className="mt-3"
                variant="secondary"
                onClick={protocolRecovery ? reloadFromFirstPage : () => void publicationsQuery.refetch()}
              >
                {protocolRecovery ? 'Reload from first page' : 'Retry publications'}
              </Button>
            </StateMessage>
          ) : publications.length === 0 ? (
            <StateMessage title={hasFilters ? 'No artifact publications match these filters' : 'No artifact publications yet'}>
              {hasFilters ? (
                <Button className="mt-3" variant="secondary" onClick={resetFilters}>
                  Clear filters
                </Button>
              ) : (
                'Successful artifact generations will appear here after publication indexing completes.'
              )}
            </StateMessage>
          ) : (
            <>
              {(hasError || protocolRecovery) && (
                <div
                  className="flex min-w-0 flex-wrap items-center justify-between gap-2 border border-[var(--color-error)] bg-[var(--color-error-bg)] p-2 text-xs"
                  role="alert"
                >
                  <span>
                    {protocolRecovery
                      ? 'Pagination could not be confirmed. Previously confirmed rows are retained.'
                      : publicationsQuery.isFetchNextPageError
                        ? 'The next page could not be loaded. Previously loaded rows are retained.'
                        : 'Refresh failed. Previously loaded publication rows are retained.'}
                  </span>
                  <Button
                    className="px-2 py-1 text-xs"
                    variant="secondary"
                    onClick={
                      protocolRecovery
                        ? reloadFromFirstPage
                        : publicationsQuery.isFetchNextPageError
                          ? () => void publicationsQuery.fetchNextPage()
                          : () => void publicationsQuery.refetch()
                    }
                  >
                    {protocolRecovery
                      ? 'Reload from first page'
                      : publicationsQuery.isFetchNextPageError
                        ? 'Retry page'
                        : 'Retry'}
                  </Button>
                </div>
              )}
              <p className="text-xs text-[var(--color-text-muted)]" role="status">
                {publications.length} publication{publications.length === 1 ? '' : 's'} loaded
              </p>
              <ArtifactPublicationList
                publications={publications}
                hasNextPage={Boolean(publicationsQuery.hasNextPage) && !protocolRecovery}
                isFetchingNextPage={publicationsQuery.isFetchingNextPage}
                onLoadMore={() => void publicationsQuery.fetchNextPage()}
              />
            </>
          )}
        </div>
      </Panel>
    </section>
  );
}

export function ArtifactPublicationDetailPage() {
  const { runId, artifactId } = useParams<{
    runId: string;
    artifactId: string;
  }>();
  const [searchParams] = useSearchParams();
  const identity = parseArtifactPublicationDetailIdentity(
    runId,
    artifactId,
    searchParams,
  );
  if (!identity.ok) {
    return (
      <section className="flex min-w-0 flex-1 flex-col gap-3">
        <Panel title="Artifact publication">
          <StateMessage title="Artifact publication link could not be used" tone="error">
            The URL does not contain one valid run, artifact, and generation identity.
            <Button asChild className="mt-3" variant="secondary">
              <Link to="/artifacts">Return to artifacts</Link>
            </Button>
          </StateMessage>
        </Panel>
      </section>
    );
  }
  return (
    <ArtifactPublicationDetail
      key={`${identity.runId}:${identity.artifactId}:${identity.generation}`}
      runId={identity.runId}
      artifactId={identity.artifactId}
      generation={identity.generation}
    />
  );
}

interface FilterDraft {
  projectId: string;
  runId: string;
  workflowId: string;
  outputType: string;
  associatedRunSampleRevisionId: string;
  publishedFrom: string;
  publishedBefore: string;
  currentOnly: boolean;
  limit: string;
}

function ArtifactPublicationFilterForm({
  filters,
  onApply,
  onReset,
}: {
  filters: ArtifactPublicationFilters;
  onApply: (filters: ArtifactPublicationFilters) => void;
  onReset: () => void;
}) {
  const [draft, setDraft] = useState<FilterDraft>({
    projectId: filters.projectId ?? '',
    runId: filters.runId ?? '',
    workflowId: filters.workflowId ?? '',
    outputType: filters.outputType ?? '',
    associatedRunSampleRevisionId: filters.associatedRunSampleRevisionId ?? '',
    publishedFrom: filters.publishedFrom ?? '',
    publishedBefore: filters.publishedBefore ?? '',
    currentOnly: filters.currentOnly,
    limit: String(filters.limit),
  });
  const [formError, setFormError] = useState(false);

  function update<K extends keyof FilterDraft>(key: K, value: FilterDraft[K]) {
    setDraft((current) => ({ ...current, [key]: value }));
    setFormError(false);
  }

  function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const candidate = new URLSearchParams();
    setDraftValue(candidate, 'project_id', draft.projectId);
    setDraftValue(candidate, 'run_id', draft.runId);
    setDraftValue(candidate, 'workflow_id', draft.workflowId);
    setDraftValue(candidate, 'output_type', draft.outputType);
    setDraftValue(
      candidate,
      'associated_run_sample_revision_id',
      draft.associatedRunSampleRevisionId,
    );
    setDraftValue(candidate, 'published_from', draft.publishedFrom);
    setDraftValue(candidate, 'published_before', draft.publishedBefore);
    if (!draft.currentOnly) candidate.set('current_only', 'false');
    if (draft.limit !== String(ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE)) {
      candidate.set('limit', draft.limit);
    }
    const parsed = parseArtifactPublicationFilters(candidate);
    if (!parsed.ok) {
      setFormError(true);
      return;
    }
    setFormError(false);
    onApply(parsed.filters);
  }

  return (
    <form
      className="border-t border-[var(--color-border)] py-3"
      aria-label="Artifact publication filters"
      onSubmit={submit}
    >
      <div className="grid min-w-0 gap-2 sm:grid-cols-2 xl:grid-cols-4">
        <FilterField
          label="Project"
          value={draft.projectId}
          onChange={(value) => update('projectId', value)}
          placeholder="prj_…"
        />
        <FilterField
          label="Run"
          value={draft.runId}
          onChange={(value) => update('runId', value)}
        />
        <FilterField
          label="Workflow"
          value={draft.workflowId}
          onChange={(value) => update('workflowId', value)}
        />
        <FilterField
          label="Artifact type"
          value={draft.outputType}
          onChange={(value) => update('outputType', value)}
        />
        <FilterField
          label="Associated run sample revision"
          value={draft.associatedRunSampleRevisionId}
          onChange={(value) => update('associatedRunSampleRevisionId', value)}
          placeholder="smpr_…"
        />
        <FilterField
          label="Published from"
          value={draft.publishedFrom}
          onChange={(value) => update('publishedFrom', value)}
          placeholder="2026-08-07T08:00:00Z"
        />
        <FilterField
          label="Published before"
          value={draft.publishedBefore}
          onChange={(value) => update('publishedBefore', value)}
          placeholder="2026-08-08T08:00:00Z"
        />
        <label className="min-w-0 text-xs font-medium">
          Rows per page
          <select
            className="mt-1 block h-9 w-full min-w-0 rounded border border-[var(--color-border)] bg-[var(--color-surface)] px-2 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
            aria-label="Rows per page"
            value={draft.limit}
            onChange={(event) => update('limit', event.target.value)}
          >
            {[25, 50, 100].map((value) => (
              <option key={value} value={value}>{value}</option>
            ))}
          </select>
        </label>
      </div>
      <div className="mt-3 flex min-w-0 flex-wrap items-center justify-between gap-2">
        <label className="inline-flex items-center gap-2 text-xs font-medium">
          <input
            type="checkbox"
            checked={draft.currentOnly}
            onChange={(event) => update('currentOnly', event.target.checked)}
          />
          Current publications only
        </label>
        <div className="flex flex-wrap gap-2">
          <Button type="button" variant="secondary" onClick={onReset}>
            Reset filters
          </Button>
          <Button type="submit" variant="primary">Apply filters</Button>
        </div>
      </div>
      <p className="mt-2 text-xs text-[var(--color-text-muted)]">
        Publication times use timezone-aware RFC 3339 values. Identity filters are
        exact matches; the time interval includes Published from and excludes Published
        before.
      </p>
      {formError && (
        <p className="mt-2 text-xs text-[var(--color-error)]" role="alert">
          Filters could not be applied. Check opaque identities, timestamps, and the publication interval.
        </p>
      )}
    </form>
  );
}

function FilterField({
  label,
  value,
  onChange,
  placeholder,
}: {
  label: string;
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
}) {
  return (
    <label className="min-w-0 text-xs font-medium">
      {label}
      <input
        className="mt-1 block h-9 w-full min-w-0 rounded border border-[var(--color-border)] bg-[var(--color-surface)] px-2 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
        aria-label={label}
        value={value}
        placeholder={placeholder}
        autoComplete="off"
        onChange={(event) => onChange(event.target.value)}
      />
    </label>
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

function PublicationListSkeleton() {
  return (
    <div className="min-h-64 animate-pulse space-y-2" aria-label="Loading artifact publications">
      {[0, 1, 2, 3].map((row) => (
        <div key={row} className="h-14 rounded bg-slate-100" />
      ))}
    </div>
  );
}

function setDraftValue(
  params: URLSearchParams,
  key: string,
  value: string,
): void {
  if (value !== '') params.set(key, value);
}

function activeFilterEntries(
  filters: ArtifactPublicationFilters,
): Array<[string, string]> {
  const entries: Array<[string, string | null]> = [
    ['Project', filters.projectId],
    ['Run', filters.runId],
    ['Workflow', filters.workflowId],
    ['Artifact type', filters.outputType],
    ['Sample revision', filters.associatedRunSampleRevisionId],
    ['Published from', filters.publishedFrom],
    ['Published before', filters.publishedBefore],
  ];
  const active = entries.filter(
    (entry): entry is [string, string] => entry[1] !== null,
  );
  if (!filters.currentOnly) active.push(['Publication state', 'all']);
  if (filters.limit !== ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE) {
    active.push(['Rows per page', String(filters.limit)]);
  }
  return active;
}
