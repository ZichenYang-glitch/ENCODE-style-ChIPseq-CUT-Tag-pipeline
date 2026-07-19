import {
  useInfiniteQuery,
  useMutation,
  useQuery,
  useQueryClient,
} from '@tanstack/react-query';
import { RefreshCw } from 'lucide-react';
import { useEffect, useRef, useState, type ReactNode } from 'react';
import type {
  ArtifactReferenceResponse,
  RunArtifactDetailResponse,
  RunArtifactsResponse,
} from '../../api/generated/models';
import {
  downloadRunArtifact,
  getRunArtifact,
  listRunArtifacts,
} from '../../api/generated/artifacts/artifacts';
import { Button } from '../../components/Button';
import { ApiError } from '../../api/fetcher';
import { ArtifactInspector } from './ArtifactInspector';
import { ArtifactList } from './ArtifactList';
import {
  flattenArtifactPages,
  isValidArtifactCursor,
  isValidArtifactGeneration,
  isValidArtifactId,
  safeNextArtifactCursor,
  type ArtifactExtractionOutcome,
} from './artifactState';
import {
  safeArtifactDownloadFilename,
  saveArtifactBlob,
} from './artifactDownload';

const ARTIFACT_PAGE_SIZE = 50;

interface ArtifactBrowserProps {
  runId: string;
  runStatus: string;
  outcome: ArtifactExtractionOutcome;
  selectedArtifactId: string | null;
  onSelectArtifact: (artifactId: string) => void;
  onRefreshStatus: () => void | Promise<void>;
}

export function artifactListQueryKey(runId: string, generation: string | null) {
  return ['run-artifacts', runId, generation] as const;
}

export function artifactDetailQueryKey(
  runId: string,
  generation: string | null,
  artifactId: string | null,
) {
  return ['run-artifact', runId, generation, artifactId] as const;
}

function validateListResponse(
  response: RunArtifactsResponse,
  runId: string,
  generation: string,
): RunArtifactsResponse {
  if (
    response.ok !== true ||
    response.run_id !== runId ||
    !isValidArtifactGeneration(response.artifact_generation) ||
    response.artifact_generation !== generation ||
    !Array.isArray(response.artifacts) ||
    response.artifacts.some(
      (artifact) =>
        artifact.run_id !== runId || !isValidArtifactId(artifact.artifact_id),
    ) ||
    (response.next_cursor != null && !isValidArtifactCursor(response.next_cursor))
  ) {
    throw new Error('invalid artifact list response');
  }
  return response;
}

function validateDetailResponse(
  response: RunArtifactDetailResponse,
  runId: string,
  artifactId: string,
  generation: string,
): ArtifactReferenceResponse {
  const artifact = response.artifact;
  if (
    response.ok !== true ||
    response.run_id !== runId ||
    !isValidArtifactGeneration(response.artifact_generation) ||
    response.artifact_generation !== generation ||
    !artifact ||
    artifact.run_id !== runId ||
    artifact.artifact_id !== artifactId
  ) {
    throw new Error('invalid artifact detail response');
  }
  return artifact;
}

function isArtifactGenerationChanged(error: unknown): boolean {
  return (
    error instanceof ApiError &&
    error.code === 'RUN_ARTIFACT_GENERATION_CHANGED'
  );
}

function isArtifactDownloadGenerationChanged(error: unknown): boolean {
  return (
    error instanceof ApiError &&
    error.status === 409 &&
    error.code === 'RUN_ARTIFACT_DOWNLOAD_CONFLICT' &&
    error.issues.some((issue) => issue.path === 'generation')
  );
}

export function ArtifactBrowser({
  runId,
  runStatus,
  outcome,
  selectedArtifactId,
  onSelectArtifact,
  onRefreshStatus,
}: ArtifactBrowserProps) {
  const queryClient = useQueryClient();
  const indexed = runStatus === 'succeeded' && outcome.kind === 'indexed';
  const generation = outcome.kind === 'indexed' ? outcome.generation : null;
  const [staleGeneration, setStaleGeneration] = useState<string | null>(null);
  const staleRunId = useRef<string | null>(null);
  const [statusRefreshFailed, setStatusRefreshFailed] = useState(false);
  const automaticRefreshGeneration = useRef<string | null>(null);
  const generationIsStale =
    generation !== null &&
    staleRunId.current === runId &&
    generation === staleGeneration;
  const validSelection = isValidArtifactId(selectedArtifactId);
  const listQueryKey = artifactListQueryKey(runId, generation);
  const detailQueryKey = artifactDetailQueryKey(
    runId,
    generation,
    selectedArtifactId,
  );

  const listQuery = useInfiniteQuery({
    queryKey: listQueryKey,
    enabled: indexed && !generationIsStale,
    initialPageParam: undefined as string | undefined,
    queryFn: async ({ pageParam }) => {
      if (!generation) throw new Error('artifact generation is unavailable');
      return validateListResponse(
        await listRunArtifacts(runId, {
          after: pageParam,
          generation,
          limit: ARTIFACT_PAGE_SIZE,
        }),
        runId,
        generation,
      );
    },
    getNextPageParam: (lastPage, allPages, _lastPageParam, allPageParams) =>
      safeNextArtifactCursor(lastPage, allPages, allPageParams),
    retry: (failureCount, error) =>
      !(error instanceof ApiError && error.status === 409) && failureCount < 1,
    retryDelay: 100,
  });

  const detailQuery = useQuery({
    queryKey: detailQueryKey,
    enabled: indexed && validSelection && !generationIsStale,
    queryFn: async () => {
      if (!generation) throw new Error('artifact generation is unavailable');
      return validateDetailResponse(
        await getRunArtifact(runId, selectedArtifactId!, { generation }),
        runId,
        selectedArtifactId!,
        generation,
      );
    },
    retry: (failureCount, error) =>
      !(error instanceof ApiError && error.status === 409) && failureCount < 1,
    retryDelay: 100,
  });

  const downloadMutation = useMutation({
    mutationFn: async ({
      artifact,
      viewedGeneration,
    }: {
      artifact: ArtifactReferenceResponse;
      viewedGeneration: string;
    }) => {
      const blob = await downloadRunArtifact(runId, artifact.artifact_id, {
        generation: viewedGeneration,
        revision: artifact.revision,
      });
      saveArtifactBlob(
        blob,
        safeArtifactDownloadFilename(artifact.name, artifact.artifact_id),
      );
      return artifact.artifact_id;
    },
  });
  const resetDownload = downloadMutation.reset;

  useEffect(() => {
    queryClient.removeQueries({
      queryKey: ['run-artifacts', runId],
      predicate: (query) => query.queryKey[2] !== generation,
    });
    queryClient.removeQueries({
      queryKey: ['run-artifact', runId],
      predicate: (query) => query.queryKey[2] !== generation,
    });
    resetDownload();
    const staleBelongsToAnotherRun =
      staleGeneration !== null && staleRunId.current !== runId;
    const newerGenerationArrived =
      staleGeneration !== null &&
      generation !== null &&
      staleGeneration !== generation;
    if (staleBelongsToAnotherRun || newerGenerationArrived) {
      setStaleGeneration(null);
      staleRunId.current = null;
      setStatusRefreshFailed(false);
      automaticRefreshGeneration.current = null;
    }
  }, [generation, queryClient, resetDownload, runId, staleGeneration]);

  const generationChanged =
    isArtifactGenerationChanged(listQuery.error) ||
    isArtifactGenerationChanged(detailQuery.error) ||
    (downloadMutation.variables?.viewedGeneration === generation &&
      isArtifactDownloadGenerationChanged(downloadMutation.error));

  useEffect(() => {
    if (!generationChanged || !generation) return;
    if (automaticRefreshGeneration.current === generation) return;
    automaticRefreshGeneration.current = generation;
    staleRunId.current = runId;
    setStaleGeneration(generation);
    setStatusRefreshFailed(false);
    void Promise.resolve()
      .then(onRefreshStatus)
      .catch(() => setStatusRefreshFailed(true));
  }, [generation, generationChanged, onRefreshStatus, runId]);

  useEffect(() => {
    if (!generationIsStale) return;
    queryClient.removeQueries({
      queryKey: ['run-artifacts', runId, generation],
      exact: true,
    });
    queryClient.removeQueries({
      queryKey: ['run-artifact', runId, generation],
      predicate: (query) => query.queryKey[2] === generation,
    });
    resetDownload();
  }, [generation, generationIsStale, queryClient, resetDownload, runId]);

  async function retryCanonicalStatus(): Promise<void> {
    setStatusRefreshFailed(false);
    try {
      await onRefreshStatus();
    } catch {
      setStatusRefreshFailed(true);
    }
  }

  if (runStatus !== 'succeeded') {
    return (
      <ArtifactStatus title="Artifacts are not available yet">
        Artifacts become available after the run completes successfully.
      </ArtifactStatus>
    );
  }

  if (outcome.kind === 'pending') {
    return (
      <ArtifactStatus title="Indexing artifacts" busy>
        The workflow succeeded. Waiting for the worker to finish indexing persisted artifact metadata.
      </ArtifactStatus>
    );
  }

  if (outcome.kind === 'unconfirmed') {
    return (
      <ArtifactStatus title="Artifact status could not be confirmed">
        The available run history is incomplete or malformed. Refresh the canonical run status before continuing.
        <StatusRefreshButton onClick={onRefreshStatus} />
      </ArtifactStatus>
    );
  }

  if (outcome.kind === 'failed') {
    return (
      <ArtifactStatus title="Artifact indexing failed" tone="error">
        The run remains succeeded, but its artifact metadata could not be indexed.
        {outcome.reasonCode && (
          <span className="mt-1 block font-mono text-xs">Reference: {outcome.reasonCode}</span>
        )}
        <StatusRefreshButton onClick={onRefreshStatus} />
      </ArtifactStatus>
    );
  }

  if (generationChanged || generationIsStale) {
    return (
      <ArtifactStatus title="Artifact generation changed" tone="error">
        Cached artifact pages and details were discarded. Refresh canonical run status before loading the replacement generation.
        {statusRefreshFailed && (
          <span className="mt-1 block">
            The status refresh did not complete. It is safe to retry.
          </span>
        )}
        <StatusRefreshButton onClick={retryCanonicalStatus} />
      </ArtifactStatus>
    );
  }

  const pages = listQuery.data?.pages ?? [];
  let artifacts: ArtifactReferenceResponse[] = [];
  let pageIntegrityFailed = false;
  try {
    artifacts = flattenArtifactPages(pages);
  } catch {
    pageIntegrityFailed = true;
  }
  const hasCachedArtifacts = artifacts.length > 0;

  if (listQuery.isPending) {
    return <ArtifactListSkeleton />;
  }

  if ((listQuery.isError || pageIntegrityFailed) && !hasCachedArtifacts) {
    return (
      <ArtifactStatus title="Artifacts could not be loaded" tone="error">
        The persisted artifact index is temporarily unavailable. Existing run state was not changed.
        <Button className="mt-3" variant="secondary" onClick={() => void listQuery.refetch()}>
          Retry artifact list
        </Button>
      </ArtifactStatus>
    );
  }

  if (outcome.count === 0 && artifacts.length === 0 && !listQuery.isError) {
    return (
      <ArtifactStatus title="No indexed artifacts">
        Artifact indexing completed successfully and produced an empty result set.
      </ArtifactStatus>
    );
  }

  if (artifacts.length === 0) {
    return (
      <ArtifactStatus title="Artifact index is not ready">
        The indexing event and persisted artifact list do not yet agree. Refresh the artifact list before treating this run as empty.
        <Button className="mt-3" variant="secondary" onClick={() => void listQuery.refetch()}>
          Retry artifact list
        </Button>
      </ArtifactStatus>
    );
  }

  return (
    <div className="min-w-0 space-y-3" data-testid="artifact-browser">
      <div className="flex min-w-0 flex-wrap items-baseline justify-between gap-2">
        <div>
          <h3 className="text-sm font-semibold">Indexed artifacts</h3>
          <p className="text-xs text-[var(--color-text-muted)]">
            {artifacts.length} loaded{outcome.count >= artifacts.length ? ` of ${outcome.count}` : ''}
          </p>
        </div>
      </div>

      {listQuery.isError && hasCachedArtifacts && (
        <div
          className="flex flex-wrap items-center justify-between gap-2 border-y border-[var(--color-warning)] bg-[var(--color-warning-bg)] px-2 py-2 text-sm text-[var(--color-warning)]"
          role="status"
        >
          <span>Additional artifact data could not be loaded. Existing rows are preserved.</span>
          <Button variant="secondary" onClick={() => void listQuery.refetch()}>
            Retry artifact list
          </Button>
        </div>
      )}

      <div className="grid min-w-0 gap-4 xl:grid-cols-[minmax(0,3fr)_minmax(18rem,1fr)]">
        <ArtifactList
          artifacts={artifacts}
          selectedArtifactId={validSelection ? selectedArtifactId : null}
          onSelect={onSelectArtifact}
          hasNextPage={Boolean(listQuery.hasNextPage)}
          isFetchingNextPage={listQuery.isFetchingNextPage}
          onLoadMore={() => void listQuery.fetchNextPage()}
        />
        <ArtifactInspector
          artifact={detailQuery.data ?? null}
          selectedArtifactId={selectedArtifactId}
          isLoading={detailQuery.isPending && validSelection}
          isError={detailQuery.isError}
          invalidSelection={selectedArtifactId !== null && !validSelection}
          onRetry={() => void detailQuery.refetch()}
          onDownload={(artifact) => {
            if (generation) {
              downloadMutation.mutate({
                artifact,
                viewedGeneration: generation,
              });
            }
          }}
          isDownloading={
            downloadMutation.isPending &&
            downloadMutation.variables?.artifact.artifact_id ===
            detailQuery.data?.artifact_id
          }
          downloadStatus={
            downloadMutation.variables?.artifact.artifact_id !==
            detailQuery.data?.artifact_id
              ? 'idle'
              : downloadMutation.isError
                ? 'error'
                : downloadMutation.isSuccess
                  ? 'success'
                  : 'idle'
          }
          onBackToList={() => {
            const list = document.getElementById('artifact-list');
            if (list && typeof list.scrollIntoView === 'function') {
              list.scrollIntoView({ block: 'start' });
            }
          }}
        />
      </div>
    </div>
  );
}

function ArtifactStatus({
  title,
  children,
  tone = 'neutral',
  busy = false,
}: {
  title: string;
  children: ReactNode;
  tone?: 'neutral' | 'error';
  busy?: boolean;
}) {
  return (
    <div
      className={`min-h-32 border-y border-[var(--color-border)] px-2 py-5 ${
        tone === 'error' ? 'text-[var(--color-error)]' : 'text-[var(--color-text-muted)]'
      }`}
      role="status"
      aria-live="polite"
      aria-busy={busy || undefined}
    >
      <h3 className="text-sm font-semibold text-[var(--color-text)]">{title}</h3>
      <div className="mt-1 max-w-2xl text-sm">{children}</div>
    </div>
  );
}

function StatusRefreshButton({
  onClick,
}: {
  onClick: () => void | Promise<void>;
}) {
  return (
    <Button
      className="mt-3"
      variant="secondary"
      onClick={() => void onClick()}
    >
      <RefreshCw className="mr-1.5 h-4 w-4" aria-hidden="true" />
      Refresh status
    </Button>
  );
}

function ArtifactListSkeleton() {
  return (
    <div className="min-h-64 space-y-2" role="status" aria-label="Loading artifacts">
      <div className="h-8 animate-pulse rounded bg-gray-100" />
      {Array.from({ length: 5 }, (_, index) => (
        <div key={index} className="h-12 animate-pulse border-b border-[var(--color-border)] bg-gray-50" />
      ))}
    </div>
  );
}
