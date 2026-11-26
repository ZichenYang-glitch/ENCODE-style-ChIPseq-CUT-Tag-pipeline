import { useEffect, useState, type FormEvent } from 'react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { authoringListProjects } from '../../api/generated/authoring/authoring';
import {
  catalogCreateExperiment,
  catalogCreateTerm,
  catalogListExperimentRevisionRunResults,
  catalogListExperimentRevisions,
  catalogListExperiments,
  catalogListSampleRevisionRunResults,
  catalogListSamples,
  catalogListTerms,
  catalogReviseExperiment,
} from '../../api/generated/catalog/catalog';
import type {
  CatalogExperimentResponse,
  CatalogExperimentRunListResponse,
  CatalogSampleListResponseSamplesItem,
  CatalogSampleRunListResponseResultsItem,
  CatalogTermSummaryResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { useAuth } from '../../app/auth';
import { Button } from '../../components/Button';
import { PeriTargetPanel } from './PeriTargetPanel';
import {
  CATALOG_PAGE_SIZE,
  RUN_INPUT_SCOPE_TEXT,
  STABLE_KEY,
  collectPages,
  experimentCellTypeLabel,
  editingHead,
  experimentListParams,
  experimentTraceText,
  knownTermSlot,
  revisionChoices,
  newIdempotencyHeaders,
  samplesForWrite,
  singleSlotLabel,
  writeSlot,
  type WriteSlot,
} from './catalogModel';

const FIELD_CLASS =
  'mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-3 text-sm focus:border-[var(--color-focus)] focus:outline-none sm:h-9';

const TERM_TYPES = [
  'assay',
  'target',
  'batch',
  'replicate_role',
  'cell_type',
  'species',
  'tissue',
  'cell_line',
  'treatment',
  'tag',
] as const;

interface AppliedSearch {
  projectId: string;
  assayTermIds: string[];
  cellTypeTermIds: string[];
}

function errorMessage(caught: unknown): string {
  if (caught instanceof ApiError) {
    if (caught.code === 'catalog_revision_conflict') {
      return '冲突。这次提交没有写入。';
    }
    return caught.message;
  }
  return 'The catalog request failed.';
}

function termsOf(
  terms: readonly CatalogTermSummaryResponse[],
  termType: string,
): CatalogTermSummaryResponse[] {
  return terms.filter(
    (term) => term.term_type === termType && term.archived_at == null,
  );
}

function isPublishedRun(
  row: CatalogSampleRunListResponseResultsItem,
): row is Extract<
  CatalogSampleRunListResponseResultsItem,
  { publication_evidence: 'present' }
> {
  return row.publication_evidence === 'present';
}

export function ExperimentCatalogPage() {
  const { loading, principal } = useAuth();
  const queryClient = useQueryClient();
  const [projectId, setProjectId] = useState('');
  const [assayTermIds, setAssayTermIds] = useState<string[]>([]);
  const [cellTypeTermIds, setCellTypeTermIds] = useState<string[]>([]);
  const [applied, setApplied] = useState<AppliedSearch | null>(null);
  const [pageIndex, setPageIndex] = useState(0);
  const [cursors, setCursors] = useState<Array<string | undefined>>([undefined]);
  const [selectedExperimentId, setSelectedExperimentId] = useState<string | null>(
    null,
  );
  const [selectedRevisionId, setSelectedRevisionId] = useState<string | null>(null);
  const [viewedRevision, setViewedRevision] = useState<CatalogExperimentResponse | null>(
    null,
  );
  const [pinnedHead, setPinnedHead] = useState<CatalogExperimentResponse | null>(null);
  const [followHeadAfter, setFollowHeadAfter] = useState<number | null>(null);
  const [historyIndex, setHistoryIndex] = useState(0);
  const [historyCursors, setHistoryCursors] = useState<Array<string | undefined>>([
    undefined,
  ]);
  const [notice, setNotice] = useState<string | null>(null);
  const [formError, setFormError] = useState<string | null>(null);

  const projectsQuery = useQuery({
    queryKey: ['catalog', 'projects'],
    queryFn: () =>
      collectPages(async (cursor) => {
        const page = await authoringListProjects(
          cursor === undefined ? { limit: 50 } : { limit: 50, cursor },
        );
        return { items: page.projects, nextCursor: page.next_cursor ?? null };
      }),
  });
  const termsQuery = useQuery({
    queryKey: ['catalog', 'terms'],
    queryFn: () =>
      collectPages(async (cursor) => {
        const page = await catalogListTerms(
          cursor === undefined ? { limit: 50 } : { limit: 50, cursor },
        );
        return { items: page.terms, nextCursor: page.next_cursor ?? null };
      }),
  });

  const projects = projectsQuery.data?.items ?? [];
  const terms = termsQuery.data?.items ?? [];
  const assayTerms = termsOf(terms, 'assay');
  const cellTypeTerms = termsOf(terms, 'cell_type');

  useEffect(() => {
    const loaded = projectsQuery.data?.items ?? [];
    if (projectId === '' && loaded.length === 1) {
      setProjectId(loaded[0]?.project_id ?? '');
    }
  }, [projectId, projectsQuery.data]);

  const pageCursor = cursors[pageIndex];
  const experimentsQuery = useQuery({
    queryKey: ['catalog', 'experiments', applied, pageCursor ?? null],
    enabled: applied !== null,
    queryFn: () => {
      if (applied === null) throw new Error('Search is not applied.');
      return catalogListExperiments(
        applied.projectId,
        experimentListParams({
          cursor: pageCursor,
          assayTermIds: applied.assayTermIds,
          cellTypeTermIds: applied.cellTypeTermIds,
        }),
      );
    },
  });
  const experiments = experimentsQuery.data?.experiments ?? [];
  const listHead =
    experiments.find((item) => item.experiment_id === selectedExperimentId) ?? null;
  const head = editingHead(listHead, pinnedHead);
  const experimentsUpdatedAt = experimentsQuery.dataUpdatedAt;
  const historyCursor = historyCursors[historyIndex];
  const historyQuery = useQuery({
    queryKey: ['catalog', 'revisions', selectedExperimentId, historyCursor ?? null],
    enabled: selectedExperimentId !== null,
    queryFn: () => {
      if (selectedExperimentId === null) throw new Error('No experiment.');
      return catalogListExperimentRevisions(
        selectedExperimentId,
        historyCursor === undefined
          ? { limit: CATALOG_PAGE_SIZE }
          : { limit: CATALOG_PAGE_SIZE, cursor: historyCursor },
      );
    },
  });
  const historyPage = historyQuery.data?.experiment_revisions ?? [];
  const historyNext = historyQuery.data?.next_cursor ?? null;
  const revisions = revisionChoices(
    [head, viewedRevision?.experiment_id === selectedExperimentId ? viewedRevision : null],
    historyPage,
  );
  const selected =
    revisions.find((item) => item.experiment_revision_id === selectedRevisionId) ??
    head;

  useEffect(() => {
    if (head === null) return;
    if (followHeadAfter !== null && experimentsUpdatedAt !== followHeadAfter) {
      setFollowHeadAfter(null);
      setPinnedHead(null);
      setViewedRevision(head);
      if (selectedRevisionId !== head.experiment_revision_id) {
        setSelectedRevisionId(head.experiment_revision_id);
      }
      return;
    }
    if (selectedRevisionId === null) {
      setViewedRevision(head);
      setSelectedRevisionId(head.experiment_revision_id);
    }
  }, [experimentsUpdatedAt, followHeadAfter, head, selectedRevisionId]);

  const traceQuery = useQuery({
    queryKey: [
      'catalog',
      'experiment-trace',
      applied?.projectId ?? null,
      selected?.experiment_revision_id ?? null,
    ],
    enabled: applied !== null && selected !== null,
    queryFn: () => {
      if (applied === null || selected === null) throw new Error('No revision.');
      return catalogListExperimentRevisionRunResults(
        applied.projectId,
        selected.experiment_revision_id,
        { limit: CATALOG_PAGE_SIZE },
      );
    },
  });
  const sampleRevisionIds = (selected?.samples ?? []).map(
    (sample) => sample.sample_revision_id,
  );
  const sampleRunsQuery = useQuery({
    queryKey: [
      'catalog',
      'sample-runs',
      applied?.projectId ?? null,
      selected?.experiment_revision_id ?? null,
      sampleRevisionIds,
    ],
    enabled: applied !== null && selected !== null && sampleRevisionIds.length > 0,
    queryFn: async () => {
      if (applied === null || selected === null) throw new Error('No revision.');
      const pages = await Promise.all(
        selected.samples.map((sample) =>
          collectPages(async (cursor) => {
            const page = await catalogListSampleRevisionRunResults(
              applied.projectId,
              sample.sample_revision_id,
              cursor === undefined
                ? { limit: CATALOG_PAGE_SIZE }
                : { limit: CATALOG_PAGE_SIZE, cursor },
            );
            return { items: page.results, nextCursor: page.next_cursor ?? null };
          }),
        ),
      );
      return {
        results: uniqueRunResults(pages.flatMap((page) => page.items)),
        truncated: pages.some((page) => page.truncated),
      };
    },
  });

  function refreshCatalog() {
    return queryClient.invalidateQueries({ queryKey: ['catalog'] });
  }

  if (loading || principal === null) {
    return (
      <p role="status" className="text-sm text-[var(--color-text-muted)]">
        Checking session…
      </p>
    );
  }

  const isAdmin = principal.role === 'administrator';
  const nextCursor = experimentsQuery.data?.next_cursor ?? null;

  function resetHistory() {
    setHistoryIndex(0);
    setHistoryCursors([undefined]);
  }

  function clearRevisionSelection() {
    setSelectedRevisionId(null);
    setViewedRevision(null);
    setPinnedHead(null);
    setFollowHeadAfter(null);
  }

  function applySearch(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (projectId === '') return;
    setApplied({
      projectId,
      assayTermIds: [...assayTermIds],
      cellTypeTermIds: [...cellTypeTermIds],
    });
    setPageIndex(0);
    setCursors([undefined]);
    setSelectedExperimentId(null);
    resetHistory();
    clearRevisionSelection();
    setNotice(null);
    setFormError(null);
  }

  function openExperiment(experimentId: string) {
    const opened =
      experiments.find((item) => item.experiment_id === experimentId) ?? null;
    setSelectedExperimentId(experimentId);
    resetHistory();
    setSelectedRevisionId(null);
    setViewedRevision(null);
    setPinnedHead(opened);
    setFollowHeadAfter(null);
    setFormError(null);
  }

  function showNewerRevisions() {
    if (historyNext == null || historyNext === '') return;
    setHistoryCursors((current) => {
      const next = current.slice(0, historyIndex + 1);
      next[historyIndex + 1] = historyNext;
      return next;
    });
    setHistoryIndex((index) => index + 1);
  }

  return (
    <section
      data-testid="experiment-catalog-page"
      className="min-w-0 flex-1 space-y-4"
    >
      <header className="space-y-1">
        <h2 className="text-xl font-semibold">Experiment catalog</h2>
        <p className="max-w-3xl text-sm text-[var(--color-text-muted)]">
          {isAdmin
            ? 'You can record and revise experiment classification.'
            : 'You can search, page through experiments, and open older revisions. Recording and revising classification is limited to administrators.'}
        </p>
      </header>

      <form
        aria-label="Find experiments"
        className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
        onSubmit={applySearch}
      >
        <div className="grid gap-3 md:grid-cols-3">
          <label className="block text-sm font-medium">
            Project
            <select
              className={FIELD_CLASS}
              value={projectId}
              onChange={(event) => setProjectId(event.currentTarget.value)}
            >
              <option value="">Select a project</option>
              {projects.map((project) => (
                <option key={project.project_id} value={project.project_id}>
                  {project.display_name}
                </option>
              ))}
            </select>
          </label>
          <label className="block text-sm font-medium">
            Assay
            <select
              multiple
              className={`${FIELD_CLASS} h-28`}
              value={assayTermIds}
              onChange={(event) =>
                setAssayTermIds(
                  Array.from(event.currentTarget.selectedOptions, (option) => option.value),
                )
              }
            >
              {assayTerms.map((term) => (
                <option key={term.term_id} value={term.term_id}>
                  {term.display_name}
                </option>
              ))}
            </select>
          </label>
          <label className="block text-sm font-medium">
            Cell type
            <select
              multiple
              className={`${FIELD_CLASS} h-28`}
              value={cellTypeTermIds}
              onChange={(event) =>
                setCellTypeTermIds(
                  Array.from(event.currentTarget.selectedOptions, (option) => option.value),
                )
              }
            >
              {cellTypeTerms.map((term) => (
                <option key={term.term_id} value={term.term_id}>
                  {term.display_name}
                </option>
              ))}
            </select>
          </label>
        </div>
        <p className="text-xs text-[var(--color-text-muted)]">
          Values in one field are OR. Different fields are AND. Leaving a field
          empty is not a search for unknown.
        </p>
        {projectsQuery.isError && (
          <p role="alert" className="text-sm text-[var(--color-error)]">
            {errorMessage(projectsQuery.error)}
          </p>
        )}
        {projectsQuery.data?.truncated && (
          <p className="text-xs text-[var(--color-text-muted)]">
            The project list stopped after 20 pages.
          </p>
        )}
        <Button type="submit" variant="primary" disabled={projectId === ''}>
          Find experiments
        </Button>
      </form>

      {applied !== null && (
        <section aria-label="Experiment results" className="space-y-3">
          <div className="flex flex-wrap items-center gap-2">
            <p className="text-sm text-[var(--color-text-muted)]">
              Page {pageIndex + 1}
            </p>
            <Button
              type="button"
              disabled={pageIndex === 0 || experimentsQuery.isFetching}
              onClick={() => setPageIndex((index) => Math.max(0, index - 1))}
            >
              Previous page
            </Button>
            <Button
              type="button"
              disabled={nextCursor == null || experimentsQuery.isFetching}
              onClick={() => {
                if (nextCursor == null) return;
                setCursors((current) => {
                  const next = current.slice(0, pageIndex + 1);
                  next[pageIndex + 1] = nextCursor;
                  return next;
                });
                setPageIndex((index) => index + 1);
              }}
            >
              Next page
            </Button>
          </div>
          {experimentsQuery.isError && (
            <p role="alert" className="text-sm text-[var(--color-error)]">
              {errorMessage(experimentsQuery.error)}
            </p>
          )}
          {experimentsQuery.isSuccess && experiments.length === 0 && (
            <p className="text-sm text-[var(--color-text-muted)]">
              No experiments match this search.
            </p>
          )}
          <ul aria-label="Experiments" className="divide-y divide-[var(--color-border)] rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)]">
            {experiments.map((experiment) => (
              <li key={experiment.experiment_id}>
                <button
                  type="button"
                  className="flex w-full flex-col items-start gap-1 px-4 py-3 text-left hover:bg-[var(--color-surface-subtle)]"
                  aria-current={
                    selectedExperimentId === experiment.experiment_id
                      ? 'true'
                      : undefined
                  }
                  onClick={() => openExperiment(experiment.experiment_id)}
                >
                  <span className="font-medium">{experiment.display_name}</span>
                  <span className="text-xs text-[var(--color-text-muted)]">
                    {singleSlotLabel(experiment.assay)} ·{' '}
                    {experimentCellTypeLabel(experiment.samples)} · 修订{' '}
                    {experiment.revision_number}
                  </span>
                </button>
              </li>
            ))}
          </ul>
        </section>
      )}

      {selected !== null && head !== null && (
        <RevisionDetail
          head={head}
          selected={selected}
          revisions={revisions}
          historyPage={historyIndex + 1}
          canShowOlderRevisions={historyIndex > 0}
          canShowNewerRevisions={
            historyNext != null && historyNext !== '' && !historyQuery.isFetching
          }
          onOlderRevisions={() => setHistoryIndex((index) => Math.max(0, index - 1))}
          onNewerRevisions={showNewerRevisions}
          onSelectRevision={(revisionId) => {
            const revision = revisions.find(
              (item) => item.experiment_revision_id === revisionId,
            );
            setFollowHeadAfter(null);
            setSelectedRevisionId(revisionId);
            if (revision !== undefined) setViewedRevision(revision);
          }}
          trace={traceQuery.data ?? null}
          traceError={traceQuery.isError ? errorMessage(traceQuery.error) : null}
          sampleRuns={sampleRunsQuery.data?.results ?? null}
          sampleRunsTruncated={sampleRunsQuery.data?.truncated ?? false}
          sampleRunsError={
            sampleRunsQuery.isError ? errorMessage(sampleRunsQuery.error) : null
          }
        />
      )}

      {isAdmin && applied && projectId === applied.projectId && head && selected &&
        selected.experiment_revision_id === head.experiment_revision_id && (
          <PeriTargetPanel
            key={head.experiment_revision_id}
            projectId={applied.projectId}
            experiment={head}
            terms={terms}
            onAdopted={(receipt) => {
              setNotice(`Peri · ${receipt.receipt_id}`);
              setPinnedHead(null);
              setFollowHeadAfter(experimentsUpdatedAt);
              void refreshCatalog();
            }}
          />
        )}

      {notice !== null && (
        <p role="status" className="text-sm text-[var(--color-text)]">
          {notice}
        </p>
      )}
      {formError !== null && (
        <p role="alert" className="text-sm text-[var(--color-error)]">
          {formError}
        </p>
      )}

      {isAdmin && (
        <AdminMaintenance
          projectId={projectId}
          terms={terms}
          head={head}
          onNotice={setNotice}
          onError={setFormError}
          onConflict={async () => {
            setPinnedHead(null);
            setFollowHeadAfter(experimentsUpdatedAt);
            await refreshCatalog();
            setFormError('冲突。这次提交没有写入。');
          }}
          onRevised={(revision) => {
            setPinnedHead(revision);
            setViewedRevision(revision);
            setSelectedRevisionId(revision.experiment_revision_id);
            setSelectedExperimentId(revision.experiment_id);
          }}
          onCreated={(revision) => {
            setSelectedExperimentId(revision.experiment_id);
            setPinnedHead(revision);
            setViewedRevision(revision);
            setSelectedRevisionId(revision.experiment_revision_id);
            resetHistory();
          }}
          refresh={refreshCatalog}
        />
      )}
    </section>
  );
}

function RevisionDetail({
  head,
  selected,
  revisions,
  historyPage,
  canShowOlderRevisions,
  canShowNewerRevisions,
  onOlderRevisions,
  onNewerRevisions,
  onSelectRevision,
  trace,
  traceError,
  sampleRuns,
  sampleRunsTruncated,
  sampleRunsError,
}: {
  head: CatalogExperimentResponse | null;
  selected: CatalogExperimentResponse;
  revisions: CatalogExperimentResponse[];
  historyPage: number;
  canShowOlderRevisions: boolean;
  canShowNewerRevisions: boolean;
  onOlderRevisions: () => void;
  onNewerRevisions: () => void;
  onSelectRevision: (revisionId: string) => void;
  trace: CatalogExperimentRunListResponse | null;
  traceError: string | null;
  sampleRuns: CatalogSampleRunListResponseResultsItem[] | null;
  sampleRunsTruncated: boolean;
  sampleRunsError: string | null;
}) {
  const traceText =
    trace === null
      ? null
      : experimentTraceText({
          trace_status: trace.trace_status,
          reason: trace.reason,
        });

  return (
    <section aria-label="Experiment revision" className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4">
      <div className="flex flex-wrap items-end justify-between gap-3">
        <div>
          <h3 className="text-lg font-semibold">{selected.display_name}</h3>
          <p className="text-xs text-[var(--color-text-muted)]">
            {revisionMark(selected, head)} {selected.revision_number}
            {' · '}
            {singleSlotLabel(selected.assay)}
            {' · '}
            {singleSlotLabel(selected.target)}
            {' · '}
            {singleSlotLabel(selected.batch)}
          </p>
        </div>
        <div className="block min-w-64 text-sm font-medium">
          <label htmlFor="experiment-revision">Revision</label>
          <select
            id="experiment-revision"
            className={FIELD_CLASS}
            value={selected.experiment_revision_id}
            onChange={(event) => onSelectRevision(event.currentTarget.value)}
          >
            {revisions.map((revision) => (
              <option
                key={revision.experiment_revision_id}
                value={revision.experiment_revision_id}
              >
                {revisionMark(revision, head)}{' '}
                {revision.revision_number}: {revision.display_name}
              </option>
            ))}
          </select>
        </div>
      </div>
      <div className="flex flex-wrap items-center gap-2">
        <p className="text-sm text-[var(--color-text-muted)]">
          Revision page {historyPage}
        </p>
        <Button
          type="button"
          disabled={!canShowOlderRevisions}
          onClick={onOlderRevisions}
        >
          Previous revision page
        </Button>
        <Button
          type="button"
          disabled={!canShowNewerRevisions}
          onClick={onNewerRevisions}
        >
          Next revision page
        </Button>
      </div>
      <ul aria-label="Pinned samples" className="space-y-1 text-sm">
        {selected.samples.map((sample) => (
          <li key={sample.sample_revision_id}>
            {sample.sample_id} · {experimentCellTypeLabel([sample])} ·{' '}
            {singleSlotLabel(sample.replicate_role)}
          </li>
        ))}
      </ul>
      <section
        aria-label="Experiment trace"
        className="rounded-[4px] border border-[var(--color-border)] bg-[var(--color-surface-subtle)] p-3"
      >
        {traceError !== null && (
          <p role="alert" className="text-sm text-[var(--color-error)]">
            {traceError}
          </p>
        )}
        {traceText !== null && (
          <>
            <h3 className="text-base font-semibold">{traceText.title}</h3>
            <p className="mt-1 text-sm">{traceText.detail}</p>
          </>
        )}
      </section>
      <section
        aria-label="Run input scope"
        className="rounded-[4px] border border-[var(--color-border)] p-3"
      >
        <h3 className="text-base font-semibold">运行输入范围</h3>
        <p className="mt-1 text-sm">{RUN_INPUT_SCOPE_TEXT}</p>
        {sampleRunsError !== null && (
          <p role="alert" className="mt-2 text-sm text-[var(--color-error)]">
            {sampleRunsError}
          </p>
        )}
        {sampleRuns !== null && sampleRuns.length === 0 && (
          <p className="mt-2 text-sm text-[var(--color-text-muted)]">
            这条样本修订没有运行结果。
          </p>
        )}
        {sampleRunsTruncated && (
          <p role="status" className="mt-2 text-sm text-[var(--color-text-muted)]">
            Run results stopped after 20 pages.
          </p>
        )}
        {sampleRuns !== null && sampleRuns.length > 0 && (
          <ul aria-label="Run results" className="mt-2 space-y-2 text-sm">
            {sampleRuns.map((row) => (
              <li key={runResultKey(row)}>
                {isPublishedRun(row) ? (
                  <>
                    <span>{row.artifact_association_kind}</span>
                    {' · '}
                    <span>{row.artifact_id}</span>
                    <ul className="mt-1 text-xs text-[var(--color-text-muted)]">
                      {row.run_sample_binding.associated_run_samples.map((sample) => (
                        <li key={sample.sample_revision_id}>
                          运行样本 {sample.sample_id} 修订 {sample.revision_number}
                        </li>
                      ))}
                    </ul>
                  </>
                ) : (
                  <span>
                    {row.run_id} · publication_evidence {row.publication_evidence}
                  </span>
                )}
              </li>
            ))}
          </ul>
        )}
      </section>
    </section>
  );
}

function runResultKey(row: CatalogSampleRunListResponseResultsItem): string {
  if (isPublishedRun(row)) {
    return `published:${row.run_id}:${row.artifact_id}:${row.artifact_generation}`;
  }
  return `unpublished:${row.run_id}:${row.provenance}`;
}

function uniqueRunResults(
  rows: readonly CatalogSampleRunListResponseResultsItem[],
): CatalogSampleRunListResponseResultsItem[] {
  const seen = new Set<string>();
  const unique: CatalogSampleRunListResponseResultsItem[] = [];
  for (const row of rows) {
    const key = runResultKey(row);
    if (seen.has(key)) continue;
    seen.add(key);
    unique.push(row);
  }
  return unique;
}

function isCurrentRevision(
  revision: { experiment_revision_id: string },
  head: CatalogExperimentResponse | null,
): boolean {
  return head !== null && revision.experiment_revision_id === head.experiment_revision_id;
}

function revisionMark(
  revision: { experiment_revision_id: string },
  head: CatalogExperimentResponse | null,
): string {
  return isCurrentRevision(revision, head) ? '当前' : '旧版';
}

function AdminMaintenance({
  projectId,
  terms,
  head,
  onNotice,
  onError,
  onConflict,
  onRevised,
  onCreated,
  refresh,
}: {
  projectId: string;
  terms: readonly CatalogTermSummaryResponse[];
  head: CatalogExperimentResponse | null;
  onNotice: (notice: string) => void;
  onError: (message: string | null) => void;
  onConflict: () => Promise<void>;
  onRevised: (revision: CatalogExperimentResponse) => void;
  onCreated: (revision: CatalogExperimentResponse) => void;
  refresh: () => Promise<void>;
}) {
  const samplesQuery = useQuery({
    queryKey: ['catalog', 'samples', projectId],
    enabled: projectId !== '',
    queryFn: () =>
      collectPages(async (cursor) => {
        const page = await catalogListSamples(
          projectId,
          cursor === undefined ? { limit: 50 } : { limit: 50, cursor },
        );
        return { items: page.samples, nextCursor: page.next_cursor ?? null };
      }),
  });
  const samples = samplesQuery.data?.items ?? [];
  const samplesStatus = samplesQuery.isError
    ? 'error'
    : samplesQuery.isSuccess
      ? 'ready'
      : 'loading';
  const assayTerms = termsOf(terms, 'assay');
  const targetTerms = termsOf(terms, 'target');
  const batchTerms = termsOf(terms, 'batch');
  const roleTerms = termsOf(terms, 'replicate_role');

  return (
    <div className="grid gap-4 xl:grid-cols-3">
      <TermForm
        onError={onError}
        onCreated={async (name) => {
          await refresh();
          onNotice(`Recorded term ${name}.`);
        }}
      />
      <ExperimentForm
        projectId={projectId}
        samples={samples}
        samplesStatus={samplesStatus}
        sampleError={
          samplesQuery.isError ? errorMessage(samplesQuery.error) : null
        }
        assayTerms={assayTerms}
        targetTerms={targetTerms}
        batchTerms={batchTerms}
        roleTerms={roleTerms}
        onError={onError}
        onCreated={async (revision) => {
          await refresh();
          onCreated(revision);
          onNotice(`Recorded experiment ${revision.display_name}.`);
        }}
      />
      <ReviseForm
        key={head?.experiment_revision_id ?? 'closed'}
        head={head}
        assayTerms={assayTerms}
        onError={onError}
        onConflict={onConflict}
        onRevised={async (revision) => {
          await refresh();
          onRevised(revision);
          onNotice(`Revised classification to ${revision.display_name}.`);
        }}
      />
    </div>
  );
}

function TermForm({
  onCreated,
  onError,
}: {
  onCreated: (displayName: string) => Promise<void>;
  onError: (message: string | null) => void;
}) {
  const [termType, setTermType] = useState<(typeof TERM_TYPES)[number]>('assay');
  const [displayName, setDisplayName] = useState('');
  const mutation = useMutation({
    mutationFn: () =>
      catalogCreateTerm(
        { term_type: termType, display_name: displayName.trim() },
        newIdempotencyHeaders(),
      ),
    onSuccess: async () => {
      const name = displayName.trim();
      setDisplayName('');
      onError(null);
      await onCreated(name);
    },
    onError: (caught) => onError(errorMessage(caught)),
  });

  return (
    <form
      aria-label="Record a term"
      className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
      onSubmit={(event) => {
        event.preventDefault();
        if (displayName.trim() === '') return;
        mutation.mutate();
      }}
    >
      <h3 className="text-sm font-semibold">Record a term</h3>
      <label className="block text-sm font-medium">
        Term type
        <select
          className={FIELD_CLASS}
          value={termType}
          onChange={(event) =>
            setTermType(event.currentTarget.value as (typeof TERM_TYPES)[number])
          }
        >
          {TERM_TYPES.map((type) => (
            <option key={type} value={type}>
              {type}
            </option>
          ))}
        </select>
      </label>
      <label className="block text-sm font-medium">
        Display name
        <input
          className={FIELD_CLASS}
          value={displayName}
          onChange={(event) => setDisplayName(event.currentTarget.value)}
          required
        />
      </label>
      <Button type="submit" variant="primary" disabled={mutation.isPending}>
        Record term
      </Button>
    </form>
  );
}

function ExperimentForm({
  projectId,
  samples,
  samplesStatus,
  sampleError,
  assayTerms,
  targetTerms,
  batchTerms,
  roleTerms,
  onCreated,
  onError,
}: {
  projectId: string;
  samples: readonly CatalogSampleListResponseSamplesItem[];
  samplesStatus: 'loading' | 'ready' | 'error';
  sampleError: string | null;
  assayTerms: readonly CatalogTermSummaryResponse[];
  targetTerms: readonly CatalogTermSummaryResponse[];
  batchTerms: readonly CatalogTermSummaryResponse[];
  roleTerms: readonly CatalogTermSummaryResponse[];
  onCreated: (revision: CatalogExperimentResponse) => Promise<void>;
  onError: (message: string | null) => void;
}) {
  const [stableKey, setStableKey] = useState('');
  const [displayName, setDisplayName] = useState('');
  const [assayId, setAssayId] = useState('');
  const [targetId, setTargetId] = useState('');
  const [batchId, setBatchId] = useState('');
  const [sampleRevisionId, setSampleRevisionId] = useState('');
  const [roleId, setRoleId] = useState('');
  const mutation = useMutation({
    mutationFn: () => {
      const assay = knownTermSlot(assayTerms, assayId);
      const target = knownTermSlot(targetTerms, targetId);
      const batch = knownTermSlot(batchTerms, batchId);
      const role = knownTermSlot(roleTerms, roleId);
      if (
        assay === null ||
        target === null ||
        batch === null ||
        role === null ||
        projectId === ''
      ) {
        throw new ApiError(400, 'catalog_payload_invalid', 'Choose a known term for each classification field.');
      }
      return catalogCreateExperiment(
        projectId,
        {
          stable_key: stableKey.trim(),
          display_name: displayName.trim(),
          assay,
          target,
          batch,
          samples: [
            { sample_revision_id: sampleRevisionId, replicate_role: role },
          ],
        },
        newIdempotencyHeaders(),
      );
    },
    onSuccess: async (revision) => {
      setStableKey('');
      setDisplayName('');
      onError(null);
      await onCreated(revision);
    },
    onError: (caught) => onError(errorMessage(caught)),
  });

  return (
    <form
      aria-label="Record an experiment"
      className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
      onSubmit={(event) => {
        event.preventDefault();
        if (!STABLE_KEY.test(stableKey.trim())) {
          onError('Stable key must start with a letter or digit and stay within the catalog key rules.');
          return;
        }
        mutation.mutate();
      }}
    >
      <h3 className="text-sm font-semibold">Record an experiment</h3>
      {projectId === '' && (
        <p className="text-xs text-[var(--color-text-muted)]">Select a project first.</p>
      )}
      {projectId !== '' && samplesStatus === 'loading' && (
        <p role="status" className="text-xs text-[var(--color-text-muted)]">
          Loading sample revisions…
        </p>
      )}
      {projectId !== '' && samplesStatus === 'error' && (
        <p role="alert" className="text-xs text-[var(--color-error)]">
          {sampleError}
        </p>
      )}
      {projectId !== '' && samplesStatus === 'ready' && samples.length === 0 && (
        <p className="text-xs text-[var(--color-text-muted)]">
          This project has no sample revisions. The catalog page cannot create samples.
        </p>
      )}
      <label className="block text-sm font-medium">
        Stable key
        <input
          className={FIELD_CLASS}
          value={stableKey}
          onChange={(event) => setStableKey(event.currentTarget.value)}
          required
        />
      </label>
      <label className="block text-sm font-medium">
        Display name
        <input
          className={FIELD_CLASS}
          value={displayName}
          onChange={(event) => setDisplayName(event.currentTarget.value)}
          required
        />
      </label>
      <TermSelect label="Assay term" terms={assayTerms} value={assayId} onChange={setAssayId} />
      <TermSelect label="Target term" terms={targetTerms} value={targetId} onChange={setTargetId} />
      <TermSelect label="Batch term" terms={batchTerms} value={batchId} onChange={setBatchId} />
      <label className="block text-sm font-medium">
        Sample revision
        <select
          className={FIELD_CLASS}
          value={sampleRevisionId}
          onChange={(event) => setSampleRevisionId(event.currentTarget.value)}
          required
        >
          <option value="">Select a sample revision</option>
          {samples.map((sample) => (
            <option key={sample.sample_revision_id} value={sample.sample_revision_id}>
              {sample.display_name} · 修订 {sample.revision_number}
            </option>
          ))}
        </select>
      </label>
      <TermSelect label="Replicate role" terms={roleTerms} value={roleId} onChange={setRoleId} />
      <Button
        type="submit"
        variant="primary"
        disabled={
          mutation.isPending ||
          projectId === '' ||
          samplesStatus !== 'ready' ||
          samples.length === 0
        }
      >
        Record experiment
      </Button>
    </form>
  );
}

function TermSelect({
  label,
  terms,
  value,
  onChange,
}: {
  label: string;
  terms: readonly CatalogTermSummaryResponse[];
  value: string;
  onChange: (termId: string) => void;
}) {
  return (
    <label className="block text-sm font-medium">
      {label}
      <select
        className={FIELD_CLASS}
        value={value}
        onChange={(event) => onChange(event.currentTarget.value)}
        required
      >
        <option value="">Select a term</option>
        {terms.map((term) => (
          <option key={term.term_id} value={term.term_id}>
            {term.display_name}
          </option>
        ))}
      </select>
    </label>
  );
}

function ReviseForm({
  head,
  assayTerms,
  onRevised,
  onConflict,
  onError,
}: {
  head: CatalogExperimentResponse | null;
  assayTerms: readonly CatalogTermSummaryResponse[];
  onRevised: (revision: CatalogExperimentResponse) => Promise<void>;
  onConflict: () => Promise<void>;
  onError: (message: string | null) => void;
}) {
  const [displayName, setDisplayName] = useState(head?.display_name ?? '');
  const [assayId, setAssayId] = useState(
    head?.assay.status === 'known' ? head.assay.term_id : '',
  );

  const copied =
    head === null
      ? null
      : {
          target: writeSlot(head.target),
          batch: writeSlot(head.batch),
          samples: samplesForWrite(head.samples),
        };
  const canCopy =
    copied !== null &&
    copied.target !== null &&
    copied.batch !== null &&
    copied.samples !== null;

  const mutation = useMutation({
    mutationFn: () => {
      const target = copied?.target ?? null;
      const batch = copied?.batch ?? null;
      const samples = copied?.samples ?? null;
      if (head === null || target === null || batch === null || samples === null) {
        throw new ApiError(
          400,
          'catalog_payload_invalid',
          'The current head has an ambiguous field. This page does not choose candidates.',
        );
      }
      const assay: WriteSlot | null = knownTermSlot(assayTerms, assayId);
      if (assay === null) {
        throw new ApiError(400, 'catalog_payload_invalid', 'Choose a known assay term.');
      }
      return catalogReviseExperiment(
        head.experiment_id,
        {
          expected_revision_number: head.revision_number,
          display_name: displayName.trim(),
          assay,
          target,
          batch,
          samples,
        },
        newIdempotencyHeaders(),
      );
    },
    onSuccess: async (revision) => {
      onError(null);
      await onRevised(revision);
    },
    onError: async (caught) => {
      if (caught instanceof ApiError && caught.code === 'catalog_revision_conflict') {
        await onConflict();
        return;
      }
      onError(errorMessage(caught));
    },
  });

  return (
    <form
      aria-label="Revise experiment classification"
      className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
      onSubmit={(event) => {
        event.preventDefault();
        mutation.mutate();
      }}
    >
      <h3 className="text-sm font-semibold">Revise the current head</h3>
      <p className="text-xs text-[var(--color-text-muted)]">
        A revision appends a new experiment revision. Older revisions stay as they were.
      </p>
      {head === null && (
        <p className="text-xs text-[var(--color-text-muted)]">
          Open an experiment before revising it.
        </p>
      )}
      {head !== null && !canCopy && (
        <p className="text-xs text-[var(--color-text-muted)]">
          The current head has an ambiguous field. This page does not choose candidates.
        </p>
      )}
      <label className="block text-sm font-medium">
        Display name
        <input
          className={FIELD_CLASS}
          value={displayName}
          onChange={(event) => setDisplayName(event.currentTarget.value)}
          required
          disabled={head === null}
        />
      </label>
      <TermSelect
        label="Assay term"
        terms={assayTerms}
        value={assayId}
        onChange={setAssayId}
      />
      <Button
        type="submit"
        variant="primary"
        disabled={mutation.isPending || head === null || !canCopy}
      >
        Revise classification
      </Button>
    </form>
  );
}
