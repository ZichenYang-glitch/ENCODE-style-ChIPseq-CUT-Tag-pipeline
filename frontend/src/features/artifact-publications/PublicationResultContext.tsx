import { useQuery } from '@tanstack/react-query';
import { Link } from 'react-router-dom';
import { readPublishedResultContext } from '../../api/generated/result-context/result-context';
import type { PublishedResultContextResponse } from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { IdWithCopy } from '../../components/IdWithCopy';

interface PublicationResultContextProps {
  projectId: string;
  runId: string;
  artifactId: string;
  generation: string;
}

const TRACE_REASON: Record<string, string> = {
  experiment_snapshot_binding_absent: 'This run has no copied experiment revision.',
  experiment_revision_evidence_inconsistent:
    'The copied experiment revision does not match stored evidence.',
};

export function PublicationResultContext({
  projectId,
  runId,
  artifactId,
  generation,
}: PublicationResultContextProps) {
  const query = useQuery({
    queryKey: ['published-result-context', projectId, runId, artifactId, generation],
    queryFn: () =>
      readPublishedResultContext(projectId, runId, artifactId, { generation }),
    retry: (failureCount, error) =>
      !(error instanceof ApiError && error.status >= 400 && error.status < 500) &&
      failureCount < 1,
    retryDelay: 100,
  });

  return (
    <section
      className="border-t border-[var(--color-border)] pt-3"
      aria-labelledby="published-result-context-heading"
    >
      <div className="flex min-w-0 flex-wrap items-baseline justify-between gap-2">
        <h3 id="published-result-context-heading" className="text-sm font-semibold">
          Frozen experiment revision
        </h3>
        <Link
          className="text-xs text-[var(--color-text-muted)] underline"
          to={`/runs/${encodeURIComponent(runId)}`}
        >
          Open run
        </Link>
      </div>
      {query.isPending ? (
        <p className="mt-2 text-xs text-[var(--color-text-muted)]" role="status">
          Loading result context.
        </p>
      ) : null}
      {query.isError ? (
        <p className="mt-2 text-xs text-[var(--color-error)]" role="alert">
          Result context could not be read.
        </p>
      ) : null}
      {query.data ? <ContextFacts context={query.data} /> : null}
    </section>
  );
}

function ContextFacts({ context }: { context: PublishedResultContextResponse }) {
  const trace = context.experiment_trace;
  const revisionId = trace.experiment_revision_id;
  const digest = trace.payload_digest;
  const frozen = trace.trace_status === 'experiment_snapshot' && revisionId && digest;
  const metrics = context.qc_metrics ?? [];
  const failures = context.preflight_failures ?? [];
  return (
    <div className="mt-2 space-y-3 text-xs">
      {frozen ? (
        <dl className="grid min-w-0 grid-cols-[8rem_minmax(0,1fr)] gap-x-3 gap-y-2">
          <dt className="text-[var(--color-text-muted)]">Revision</dt>
          <dd>
            <IdWithCopy value={revisionId} label="experiment revision ID" />
          </dd>
          <dt className="text-[var(--color-text-muted)]">Payload digest</dt>
          <dd>
            <IdWithCopy value={digest} label="experiment payload digest" />
          </dd>
        </dl>
      ) : (
        <p>{TRACE_REASON[trace.reason ?? ''] ?? 'This result could not be traced.'}</p>
      )}
      <p>
        Run status {context.run_status}
        {context.current_stage ? ` · stage ${context.current_stage}` : ''}
      </p>
      {context.reference_evidence === 'present' && context.reference_profile ? (
        <p>
          Reference {context.reference_profile.display_name} ·{' '}
          {context.reference_profile.organism} {context.reference_profile.assembly}
        </p>
      ) : (
        <p>No reference profile was recorded for this run.</p>
      )}
      {context.qc_evidence === 'present' && metrics.length > 0 ? (
        <ul className="space-y-1">
          {metrics.map((metric) => (
            <li key={metric.metric_id}>
              {metric.display_name}: {metric.value} {metric.unit}
            </li>
          ))}
        </ul>
      ) : (
        <p>No QC metrics were recorded for this run.</p>
      )}
      {context.qc_next_cursor ? <p>More QC metrics exist on the run.</p> : null}
      {context.preflight_failure_evidence === 'present' && failures.length > 0 ? (
        <ul className="space-y-1">
          {failures.map((failure) => (
            <li key={failure.event_id}>
              {failure.message}
              {failure.code ? ` (${failure.code})` : ''}
            </li>
          ))}
        </ul>
      ) : (
        <p>No preflight failure was recorded.</p>
      )}
    </div>
  );
}
