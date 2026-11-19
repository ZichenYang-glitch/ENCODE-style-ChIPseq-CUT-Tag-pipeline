import { useEffect } from 'react';
import { useQuery } from '@tanstack/react-query';
import { loadSelectableExperiments } from './experimentCatalog';
import {
  sampleRevisionSetsEqual,
  type ExperimentValidationHandoff,
} from './handoff';

interface ExperimentRevisionPanelProps {
  handoff: ExperimentValidationHandoff;
  sampleRevisionIds: readonly string[];
  experimentRevisionId: string | null;
  onSampleRevisionIdsChange: (ids: string[]) => void;
  onExperimentRevisionIdChange: (id: string | null) => void;
}

export function ExperimentRevisionPanel({
  handoff,
  sampleRevisionIds,
  experimentRevisionId,
  onSampleRevisionIdsChange,
  onExperimentRevisionIdChange,
}: ExperimentRevisionPanelProps) {
  const experiments = useQuery({
    queryKey: ['catalog-experiments', handoff.projectId],
    queryFn: () => loadSelectableExperiments(handoff.projectId),
    retry: false,
  });
  const choices = experiments.data ?? [];
  const selected =
    choices.find((item) => item.experimentRevisionId === experimentRevisionId) ??
    null;

  useEffect(() => {
    if (experimentRevisionId === null || !experiments.isSuccess) return;
    if (
      !experiments.data.some(
        (item) => item.experimentRevisionId === experimentRevisionId,
      )
    ) {
      onExperimentRevisionIdChange(null);
    }
  }, [
    experimentRevisionId,
    experiments.data,
    experiments.isSuccess,
    onExperimentRevisionIdChange,
  ]);

  return (
    <section
      className="mt-4 space-y-3 border-t border-[var(--color-border)] pt-4"
      aria-labelledby="experiment-revision-title"
      data-testid="experiment-revision-panel"
    >
      <div>
        <h3 id="experiment-revision-title" className="text-sm font-semibold">
          Experiment selection
        </h3>
        <p className="mt-1 text-xs text-[var(--color-text-muted)]">
          Carried from the saved draft. Choosing a revision does not bind it.
          A binding appears only after this page&apos;s validation returns one.
        </p>
      </div>
      <dl className="grid gap-2 text-sm sm:grid-cols-2">
        <div>
          <dt className="text-xs text-[var(--color-text-muted)]">Project</dt>
          <dd data-testid="carried-project-id">
            <code className="break-all">{handoff.projectId}</code>
          </dd>
        </div>
        <div>
          <dt className="text-xs text-[var(--color-text-muted)]">Reference profile</dt>
          <dd data-testid="carried-reference-profile">
            {handoff.referenceProfileRevisionId ?? 'Not copied'}
          </dd>
        </div>
      </dl>
      <ul className="list-disc space-y-1 pl-5 text-xs text-[var(--color-text-muted)]">
        {handoff.notes.map((note) => (
          <li key={note}>{note}</li>
        ))}
      </ul>
      <div className="space-y-2">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <p className="text-sm font-medium">Submitted sample revisions</p>
          <button
            type="button"
            className="min-h-9 rounded-[4px] border border-[var(--color-border-strong)] px-2 text-xs"
            onClick={() => onSampleRevisionIdsChange([...sampleRevisionIds].reverse())}
          >
            Reverse sample revision order
          </button>
        </div>
        <ol className="space-y-2">
          {sampleRevisionIds.map((id, index) => (
            <li key={`${index}:${sampleRevisionIds.length}`}>
              <label className="block text-xs" htmlFor={`submitted-sample-revision-${index}`}>
                Submitted sample revision {index + 1}
                <input
                  id={`submitted-sample-revision-${index}`}
                  className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 font-mono text-xs sm:h-9"
                  value={id}
                  spellCheck={false}
                  onChange={(event) => {
                    const next = [...sampleRevisionIds];
                    next[index] = event.target.value;
                    onSampleRevisionIdsChange(next);
                  }}
                />
              </label>
            </li>
          ))}
        </ol>
      </div>
      <label className="block text-sm font-medium" htmlFor="experiment-revision-choice">
        Experiment revision
        <select
          id="experiment-revision-choice"
          className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
          value={experimentRevisionId ?? ''}
          disabled={experiments.isLoading}
          onChange={(event) =>
            onExperimentRevisionIdChange(event.target.value || null)
          }
        >
          <option value="">Omit experiment</option>
          {choices.map((item) => (
            <option key={item.experimentRevisionId} value={item.experimentRevisionId}>
              {item.displayName}
            </option>
          ))}
        </select>
      </label>
      {experiments.isLoading && <p role="status">Loading experiment revisions.</p>}
      {experiments.isError && (
        <p role="alert">Experiment revisions could not be loaded. None is selected.</p>
      )}
      {experiments.isSuccess && choices.length === 0 && (
        <p role="status">This project has no selectable experiment revisions.</p>
      )}
      <p role="status" data-testid="experiment-selection-preview">
        {selectionPreview(selected, sampleRevisionIds, experimentRevisionId)}
      </p>
    </section>
  );
}

function selectionPreview(
  selected: { displayName: string; sampleRevisionIds: string[] } | null,
  submitted: readonly string[],
  experimentRevisionId: string | null,
): string {
  if (experimentRevisionId === null || selected === null) {
    return 'No experiment revision is selected. Validation will leave the snapshot unbound.';
  }
  if (sampleRevisionSetsEqual(selected.sampleRevisionIds, submitted)) {
    return `${selected.displayName} pins the same sample revisions. Order does not matter. Validation can freeze this revision.`;
  }
  return `${selected.displayName} does not pin this sample revision set. Validation will not freeze a binding.`;
}
