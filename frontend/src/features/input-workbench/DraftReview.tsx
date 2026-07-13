import { CheckCircle2, CircleAlert } from 'lucide-react';
import type { WorkbenchSchema } from './schemaContract';
import type { InputDraftController } from './useInputDraft';

interface DraftReviewProps {
  schema: WorkbenchSchema;
  draft: InputDraftController;
}

export function DraftReview({ schema, draft }: DraftReviewProps) {
  const issues = [
    !draft.configValid && 'Config does not satisfy the adapter schema.',
    draft.state.yamlIssue !== null && 'Config YAML still contains an error.',
    draft.state.configFormIssue !== null && draft.state.configFormIssue.message,
    !draft.samplesValid && 'Samples do not satisfy the adapter row schema.',
    !draft.optionsValid && 'Options do not satisfy the adapter schema.',
    draft.state.optionsFormIssue !== null && draft.state.optionsFormIssue.message,
    !draft.review.ok &&
      draft.review.code !== 'DRAFT_SAMPLE_INVALID' &&
      draft.review.message,
  ].filter((value): value is string => Boolean(value));

  return (
    <section className="min-w-0 space-y-3" aria-labelledby="draft-review-title">
      <div>
        <h3 id="draft-review-title" className="text-sm font-semibold">
          Draft review
        </h3>
        <p className="mt-1 text-xs text-[var(--color-text-muted)]">
          This is client-side structural review only. The backend has not
          scientifically validated these inputs.
        </p>
      </div>
      <div
        className={`flex items-start gap-2 rounded border px-3 py-2 text-sm ${
          draft.reviewReady
            ? 'border-emerald-200 bg-emerald-50 text-emerald-800'
            : 'border-amber-200 bg-[var(--color-warning-bg)] text-[var(--color-warning)]'
        }`}
        role="status"
      >
        {draft.reviewReady ? (
          <CheckCircle2 className="mt-0.5 shrink-0" aria-hidden="true" size={16} />
        ) : (
          <CircleAlert className="mt-0.5 shrink-0" aria-hidden="true" size={16} />
        )}
        <div className="min-w-0">
          <p className="font-medium">
            {draft.reviewReady
              ? 'Draft is structurally ready for backend validation.'
              : 'Draft needs attention before backend validation.'}
          </p>
          {(issues.length > 0 || draft.sampleTransportIssue !== null) && (
            <ul className="mt-1 list-disc space-y-0.5 pl-4">
              {issues.map((item) => (
                <li key={item} className="break-words [overflow-wrap:anywhere]">
                  {item}
                </li>
              ))}
              {draft.sampleTransportIssue !== null && (
                <li className="min-w-0 break-words [overflow-wrap:anywhere]">
                  {draft.sampleTransportIssue.message} Row{' '}
                  {draft.sampleTransportIssue.row}. Column{' '}
                  <code className="break-all font-mono">
                    {draft.sampleTransportIssue.column}
                  </code>
                  .
                </li>
              )}
            </ul>
          )}
        </div>
      </div>
      <dl className="grid gap-2 text-xs text-[var(--color-text-muted)] sm:grid-cols-3">
        <div>
          <dt className="font-medium">Schema</dt>
          <dd>{schema.contract.schema_version}</dd>
        </div>
        <div>
          <dt className="font-medium">Sample rows</dt>
          <dd>{draft.state.rows.length}</dd>
        </div>
        <div>
          <dt className="font-medium">Serialized bytes</dt>
          <dd>
            {draft.reviewPreviewAvailable && draft.review.ok
              ? draft.review.byteSize.toLocaleString()
              : 'Unavailable'}
          </dd>
        </div>
      </dl>
      <pre
        data-testid="draft-review-json"
        className="max-h-[34rem] min-h-64 max-w-full overflow-auto whitespace-pre-wrap break-words rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-3 text-xs"
      >
        {draft.reviewPreviewAvailable && draft.review.ok
          ? draft.review.serialized
          : 'A deterministic request preview is unavailable until the structural error is corrected.'}
      </pre>
    </section>
  );
}
