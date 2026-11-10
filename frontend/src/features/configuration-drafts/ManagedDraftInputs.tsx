import { ReferenceProfileSelector } from '../input-workbench/ReferenceProfileSelector';
import type { ReferenceProfileSummary } from '../../api/runTypes';
import type { ProjectSampleSummaryResponse } from '../../api/generated/models';

interface ManagedDraftInputsProps {
  samples: readonly ProjectSampleSummaryResponse[];
  selectedRevisionIds: readonly string[];
  unlistedCount: number;
  samplesLoading: boolean;
  samplesError: string | null;
  onToggleSample: (revisionId: string, selected: boolean) => void;
  onDropUnlisted: () => void;
  profiles: readonly ReferenceProfileSummary[] | null;
  referenceRevisionId: string | null;
  referenceLoading: boolean;
  referenceError: string | null;
  referenceMissing: boolean;
  onSelectReference: (revisionId: string | null) => void;
  onRefreshReferences: () => void;
  preservedSelectionCount: number;
}

export function ManagedDraftInputs({
  samples,
  selectedRevisionIds,
  unlistedCount,
  samplesLoading,
  samplesError,
  onToggleSample,
  onDropUnlisted,
  profiles,
  referenceRevisionId,
  referenceLoading,
  referenceError,
  referenceMissing,
  onSelectReference,
  onRefreshReferences,
  preservedSelectionCount,
}: ManagedDraftInputsProps) {
  const selected = new Set(selectedRevisionIds);
  return (
    <div className="space-y-3">
      <fieldset className="space-y-2">
        <legend className="text-sm font-medium">Sample revisions</legend>
        <p className="text-xs text-[var(--color-text-muted)]">
          Choose sample revisions from this project. Local file paths are not
          accepted.
        </p>
        {samplesLoading && <p role="status">Loading project samples.</p>}
        {samplesError && <p role="alert">{samplesError}</p>}
        {!samplesLoading && samples.length === 0 && !samplesError && (
          <p role="status">This project has no sample revisions yet.</p>
        )}
        <ul className="space-y-1" aria-label="Project sample revisions">
          {samples.map((sample) => (
            <li key={sample.sample_revision_id}>
              <label className="flex min-h-11 items-center gap-2 text-sm">
                <input
                  type="checkbox"
                  checked={selected.has(sample.sample_revision_id)}
                  onChange={(event) =>
                    onToggleSample(sample.sample_revision_id, event.target.checked)
                  }
                />
                {sample.display_name}
              </label>
            </li>
          ))}
        </ul>
        {unlistedCount > 0 && (
          <label className="flex min-h-11 items-center gap-2 text-sm">
            <input type="checkbox" checked onChange={onDropUnlisted} />
            {unlistedCount === 1
              ? 'A saved sample revision is no longer listed'
              : `${unlistedCount} saved sample revisions are no longer listed`}
          </label>
        )}
      </fieldset>
      <ReferenceProfileSelector
        profiles={profiles}
        selectedRevisionId={referenceMissing ? null : referenceRevisionId}
        loading={referenceLoading}
        refreshing={false}
        errorMessage={referenceError}
        selectionMessage={
          referenceMissing
            ? 'The saved reference revision is not in the current list. It is kept until you choose another.'
            : null
        }
        onSelect={onSelectReference}
        onRefresh={onRefreshReferences}
      />
      {preservedSelectionCount > 0 && (
        <p className="text-xs text-[var(--color-text-muted)]">
          This draft keeps server-managed file selections. Local paths cannot
          be added here.
        </p>
      )}
    </div>
  );
}
