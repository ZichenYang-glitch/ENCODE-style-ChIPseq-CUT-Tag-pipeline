import type { SampleColumn } from '../input-workbench/schemaContract';
import type { ProjectSampleSummaryResponse } from '../../api/generated/models';

export type InputCatalogStatus = 'unavailable' | 'loading' | 'empty' | 'error' | 'ready';

export interface AdmittedInputRevision {
  revisionId: string;
  label: string;
}

export interface WorkflowBindingPreparationProps {
  columns: readonly SampleColumn[];
  rowChoices: Readonly<Record<string, string>>;
  onChangeChoice: (columnKey: string, value: string) => void;
  samples: readonly ProjectSampleSummaryResponse[];
  controlSampleRevisionId: string | null;
  onSelectControl: (revisionId: string | null) => void;
  catalogStatus: InputCatalogStatus;
  catalogError: string | null;
  catalogItems: readonly AdmittedInputRevision[];
  selectedInputRevisionIds: readonly string[];
  unlistedInputCount: number;
  onToggleInput: (revisionId: string, selected: boolean) => void;
  hasPrevious: boolean;
  hasNext: boolean;
  onPrevious?: () => void;
  onNext?: () => void;
  onRetry?: () => void;
}

export function acceptCatalogPage<T>(
  requestedGeneration: number,
  currentGeneration: number,
  page: T,
): T | null {
  return requestedGeneration === currentGeneration ? page : null;
}

export function WorkflowBindingPreparation({
  columns,
  rowChoices,
  onChangeChoice,
  samples,
  controlSampleRevisionId,
  onSelectControl,
  catalogStatus,
  catalogError,
  catalogItems,
  selectedInputRevisionIds,
  unlistedInputCount,
  onToggleInput,
  hasPrevious,
  hasNext,
  onPrevious,
  onNext,
  onRetry,
}: WorkflowBindingPreparationProps) {
  const selectedInputs = new Set(selectedInputRevisionIds);
  const controlListed = samples.some(
    (sample) => sample.sample_revision_id === controlSampleRevisionId,
  );
  return (
    <section className="space-y-3" aria-labelledby="workflow-binding-preparation-title">
      <h2 id="workflow-binding-preparation-title" className="text-sm font-medium">
        Sample row preparation
      </h2>
      <p className="text-xs text-[var(--color-text-muted)]">
        Preparation only. These choices stay on this page and are not included
        in the save request.
      </p>
      <div className="space-y-2">
        {columns.map((column) =>
          column.enumValues ? (
            <label key={column.key} className="block text-sm">
              {column.label}
              <select
                className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
                value={rowChoices[column.key] ?? ''}
                onChange={(event) => onChangeChoice(column.key, event.target.value)}
              >
                <option value="">Select…</option>
                {column.enumValues.filter((option) => option !== '').map((option) => (
                  <option key={option} value={option}>
                    {option}
                  </option>
                ))}
              </select>
            </label>
          ) : (
            <p key={column.key} className="text-sm">
              <span className="font-medium">{column.label}. </span>
              This field is not connected to save yet.
            </p>
          ),
        )}
      </div>
      <label className="block text-sm">
        Control sample revision
        <select
          className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
          value={controlListed ? controlSampleRevisionId ?? '' : ''}
          onChange={(event) => onSelectControl(event.target.value || null)}
        >
          <option value="">No control sample revision</option>
          {samples.map((sample) => (
            <option key={sample.sample_revision_id} value={sample.sample_revision_id}>
              {sample.display_name}
            </option>
          ))}
        </select>
      </label>
      {controlSampleRevisionId && !controlListed && (
        <p role="status">The selected control sample revision is no longer listed.</p>
      )}
      <div className="space-y-2" aria-label="Admitted input revisions">
        <h3 className="text-sm font-medium">Admitted input revisions</h3>
        {catalogStatus === 'unavailable' && (
          <p role="status">The admitted input catalog is not connected yet.</p>
        )}
        {catalogStatus === 'loading' && <p role="status">Loading admitted input revisions.</p>}
        {catalogStatus === 'empty' && (
          <p role="status">This project has no admitted input revisions.</p>
        )}
        {catalogStatus === 'error' && (
          <div className="space-y-2">
            <p role="alert">{catalogError ?? 'The admitted input catalog could not be loaded.'}</p>
            <button type="button" className="min-h-11 text-sm underline" onClick={() => onRetry?.()}>
              Try the input catalog again
            </button>
          </div>
        )}
        {catalogStatus === 'ready' && catalogItems.length === 0 && (
          <p role="status">This project has no admitted input revisions.</p>
        )}
        {catalogStatus === 'ready' && catalogItems.length > 0 && (
          <ul className="space-y-1">
            {catalogItems.map((item) => (
              <li key={item.revisionId}>
                <label className="flex min-h-11 items-center gap-2 text-sm">
                  <input
                    type="checkbox"
                    checked={selectedInputs.has(item.revisionId)}
                    onChange={(event) => onToggleInput(item.revisionId, event.target.checked)}
                  />
                  {item.label}
                </label>
              </li>
            ))}
          </ul>
        )}
        {unlistedInputCount > 0 && (
          <p role="status">
            {unlistedInputCount === 1
              ? 'A selected input revision is no longer listed.'
              : `${unlistedInputCount} selected input revisions are no longer listed.`}
          </p>
        )}
        {(hasPrevious || hasNext) && (
          <div className="flex gap-2">
            <button type="button" disabled={!hasPrevious} onClick={() => onPrevious?.()}>
              Previous input revisions
            </button>
            <button type="button" disabled={!hasNext} onClick={() => onNext?.()}>
              Next input revisions
            </button>
          </div>
        )}
      </div>
    </section>
  );
}
