import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { ArrowLeft, Download, RefreshCw } from 'lucide-react';
import { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import { getArtifactPublication } from '../../api/generated/artifact-publications/artifact-publications';
import { downloadRunArtifact } from '../../api/generated/artifacts/artifacts';
import type { ArtifactPublicationResponse } from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import {
  safeArtifactDownloadFilename,
  saveArtifactBlob,
} from '../run-artifacts/artifactDownload';
import {
  artifactPublicationDetailQueryKey,
  derivedArtifactGenerationStatus,
  validateArtifactPublicationDetailResponse,
} from './artifactPublicationState';

interface ArtifactPublicationDetailProps {
  runId: string;
  artifactId: string;
  generation: string;
}

function isArtifactDownloadGenerationConflict(error: unknown): boolean {
  return (
    error instanceof ApiError &&
    error.status === 409 &&
    error.code === 'RUN_ARTIFACT_DOWNLOAD_CONFLICT' &&
    error.issues.some((issue) => issue.path === 'generation')
  );
}

const DATE_TIME_FORMAT = new Intl.DateTimeFormat(undefined, {
  dateStyle: 'medium',
  timeStyle: 'short',
});

export function ArtifactPublicationDetail({
  runId,
  artifactId,
  generation,
}: ArtifactPublicationDetailProps) {
  const queryClient = useQueryClient();
  const queryKey = artifactPublicationDetailQueryKey(runId, artifactId, generation);
  const [staleGeneration, setStaleGeneration] = useState(false);
  const detailQuery = useQuery({
    queryKey,
    queryFn: async () =>
      validateArtifactPublicationDetailResponse(
        await getArtifactPublication(runId, artifactId, { generation }),
        runId,
        artifactId,
        generation,
      ),
    retry: (failureCount, error) =>
      !(error instanceof ApiError && error.status >= 400 && error.status < 500) &&
      failureCount < 1,
    retryDelay: 100,
  });

  const downloadMutation = useMutation({
    mutationFn: async (publication: ArtifactPublicationResponse) => {
      const blob = await downloadRunArtifact(runId, artifactId, {
        generation: publication.artifact_generation,
        revision: publication.artifact_revision,
      });
      saveArtifactBlob(blob, safeArtifactDownloadFilename('', artifactId));
    },
    onError: async (error) => {
      if (!isArtifactDownloadGenerationConflict(error)) return;
      setStaleGeneration(true);
      await Promise.all([
        queryClient.invalidateQueries({
          queryKey: ['artifact-publications'],
          refetchType: 'all',
        }),
        queryClient.invalidateQueries({ queryKey, exact: true }),
      ]);
    },
  });
  const resetDownload = downloadMutation.reset;

  useEffect(() => {
    setStaleGeneration(false);
    resetDownload();
  }, [artifactId, generation, resetDownload, runId]);

  const publication = detailQuery.data;
  const status = publication
    ? derivedArtifactGenerationStatus(publication)
    : null;

  return (
    <section className="flex min-w-0 flex-1 flex-col gap-3">
      <Panel title="Artifact publication">
        <div className="min-w-0 space-y-3" data-testid="artifact-publication-detail">
          <div className="flex min-w-0 flex-wrap items-start justify-between gap-3">
            <div className="min-w-0">
              <Button asChild className="mb-2 px-2 text-xs" variant="secondary">
                <Link to="/artifacts">
                  <ArrowLeft className="mr-1.5 h-4 w-4" aria-hidden="true" />
                  All artifacts
                </Link>
              </Button>
              <h2 className="break-all text-sm font-semibold">Artifact publication details</h2>
              <p className="mt-1 break-all text-xs text-[var(--color-text-muted)]">
                {artifactId}
              </p>
            </div>
            <Button
              variant="secondary"
              onClick={() => void detailQuery.refetch()}
              disabled={detailQuery.isFetching}
              aria-label="Refresh artifact publication"
            >
              <RefreshCw
                className={`mr-1.5 h-4 w-4 ${detailQuery.isFetching ? 'animate-spin' : ''}`}
                aria-hidden="true"
              />
              Refresh
            </Button>
          </div>

          {detailQuery.isPending ? (
            <DetailSkeleton />
          ) : detailQuery.isError && !publication ? (
            <DetailState title="Artifact publication could not be loaded" tone="error">
              The requested publication is unavailable or its public identity could not be confirmed.
              <Button
                className="mt-3"
                variant="secondary"
                onClick={() => void detailQuery.refetch()}
              >
                Retry publication
              </Button>
            </DetailState>
          ) : publication ? (
            <>
              {detailQuery.isError && (
                <div
                  className="border border-[var(--color-error)] bg-[var(--color-error-bg)] p-2 text-xs"
                  role="alert"
                >
                  Refresh failed. The last confirmed publication metadata is retained.
                </div>
              )}
              <div className="flex min-w-0 flex-wrap items-start justify-between gap-3 border-y border-[var(--color-border)] py-3">
                <div className="min-w-0">
                  <p className="break-all text-base font-semibold">{publication.output_type}</p>
                  <p className="mt-1 break-all text-xs text-[var(--color-text-muted)]">
                    {publication.artifact_id}
                  </p>
                </div>
                <div className="flex shrink-0 flex-wrap items-center gap-2">
                  {staleGeneration ? (
                    <span
                      className="rounded-full border border-[var(--color-error)] bg-[var(--color-error-bg)] px-2 py-1 text-xs font-semibold text-[var(--color-error)]"
                      role="status"
                    >
                      Stale generation
                    </span>
                  ) : status === 'current' ? (
                    <span className="rounded-full border border-emerald-300 bg-emerald-50 px-2 py-1 text-xs font-semibold text-emerald-800">
                      Current
                    </span>
                  ) : (
                    <span className="text-xs font-semibold text-[var(--color-text-muted)]">
                      Superseded · historical metadata only
                    </span>
                  )}
                  {status === 'current' && (
                    <Button
                      variant="secondary"
                      onClick={() => downloadMutation.mutate(publication)}
                      disabled={downloadMutation.isPending || staleGeneration}
                      aria-label={`Download artifact ${publication.artifact_id}`}
                    >
                      <Download className="mr-1.5 h-4 w-4" aria-hidden="true" />
                      {downloadMutation.isPending ? 'Downloading…' : 'Download'}
                    </Button>
                  )}
                </div>
              </div>

              {downloadMutation.isError && !staleGeneration && (
                <p className="text-xs text-[var(--color-error)]" role="status">
                  Download could not be completed. Retry when the current artifact is available.
                </p>
              )}

              <PublicationFacts publication={publication} />
              <AssociatedRunSamples publication={publication} />
            </>
          ) : null}
        </div>
      </Panel>
    </section>
  );
}

function PublicationFacts({ publication }: { publication: ArtifactPublicationResponse }) {
  const status = derivedArtifactGenerationStatus(publication);
  return (
    <dl className="grid min-w-0 grid-cols-[8rem_minmax(0,1fr)] gap-x-3 gap-y-2 text-xs sm:grid-cols-[9rem_minmax(0,1fr)]">
      <dt className="text-[var(--color-text-muted)]">Project</dt>
      <dd className="break-all">{publication.project_id}</dd>
      <dt className="text-[var(--color-text-muted)]">Run</dt>
      <dd className="break-all">{publication.run_id}</dd>
      <dt className="text-[var(--color-text-muted)]">Workflow</dt>
      <dd className="break-all">{publication.workflow_id}</dd>
      <dt className="text-[var(--color-text-muted)]">Artifact type</dt>
      <dd className="break-all">{publication.output_type}</dd>
      <dt className="text-[var(--color-text-muted)]">Published</dt>
      <dd>
        <time dateTime={publication.published_at}>
          {DATE_TIME_FORMAT.format(new Date(publication.published_at))}
        </time>
      </dd>
      <dt className="text-[var(--color-text-muted)]">Generation</dt>
      <dd className="break-all font-mono">{publication.artifact_generation}</dd>
      <dt className="text-[var(--color-text-muted)]">Revision</dt>
      <dd className="break-all font-mono">{publication.artifact_revision}</dd>
      {status === 'superseded' && (
        <>
          <dt className="text-[var(--color-text-muted)]">Current generation</dt>
          <dd className="break-all font-mono">
            {publication.current_artifact_generation}
          </dd>
        </>
      )}
      <dt className="text-[var(--color-text-muted)]">State</dt>
      <dd>{status === 'current' ? 'Current' : 'Superseded'}</dd>
    </dl>
  );
}

function AssociatedRunSamples({ publication }: { publication: ArtifactPublicationResponse }) {
  const binding = publication.run_sample_binding;
  return (
    <section className="border-t border-[var(--color-border)] pt-3" aria-labelledby="associated-run-samples-heading">
      <div className="flex min-w-0 flex-wrap items-baseline justify-between gap-2">
        <h3 id="associated-run-samples-heading" className="text-sm font-semibold">
          Associated run samples
        </h3>
        <span className="text-xs text-[var(--color-text-muted)]">
          {binding.binding_mode} · {binding.provenance}
        </span>
      </div>
      <p className="mt-1 text-xs text-[var(--color-text-muted)]">
        These are input samples associated with the run; they are not per-artifact sample attribution.
      </p>
      {binding.associated_run_samples.length === 0 ? (
        <p className="mt-3 text-sm text-[var(--color-text-muted)]">
          No associated run samples were recorded for this run binding.
        </p>
      ) : (
        <ol className="mt-3 divide-y divide-[var(--color-border)] border-y border-[var(--color-border)]">
          {binding.associated_run_samples.map((sample) => (
            <li key={sample.sample_revision_id} className="grid min-w-0 gap-1 py-2 text-xs sm:grid-cols-[minmax(0,1fr)_minmax(0,1fr)_auto] sm:gap-3">
              <span className="break-all" title={sample.sample_id}>{sample.sample_id}</span>
              <span className="break-all font-mono text-[var(--color-text-muted)]" title={sample.sample_revision_id}>
                {sample.sample_revision_id}
              </span>
              <span className="text-[var(--color-text-muted)]">
                Revision {sample.revision_number} · position {sample.ordinal + 1}
              </span>
            </li>
          ))}
        </ol>
      )}
    </section>
  );
}

function DetailState({
  title,
  tone = 'neutral',
  children,
}: {
  title: string;
  tone?: 'neutral' | 'error';
  children: React.ReactNode;
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

function DetailSkeleton() {
  return (
    <div className="min-h-64 animate-pulse space-y-3" aria-label="Loading artifact publication">
      <div className="h-14 rounded bg-slate-100" />
      <div className="h-32 rounded bg-slate-100" />
      <div className="h-24 rounded bg-slate-100" />
    </div>
  );
}
