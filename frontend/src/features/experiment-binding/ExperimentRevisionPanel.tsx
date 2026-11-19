import { useEffect } from 'react';
import { useQuery } from '@tanstack/react-query';
import {
  authoringListProjectSamples,
  authoringListProjects,
} from '../../api/generated/authoring/authoring';
import { loadSelectableExperiments } from './experimentCatalog';
import {
  sampleRevisionSetsEqual,
  type ExperimentValidationHandoff,
} from './handoff';

const SAMPLE_REVISION_ID = /^smpr_[0-9a-f]{32}$/;
const PROJECT_ID = /^prj_[0-9a-f]{32}$/;

interface EntryProject {
  projectId: string;
  displayName: string;
}

interface EntrySample {
  sampleRevisionId: string;
  displayName: string;
}

interface ExperimentRevisionPanelProps {
  handoff: ExperimentValidationHandoff | null;
  projectId: string | null;
  sampleRevisionIds: readonly string[];
  experimentRevisionId: string | null;
  onDirectProjectChange?: (projectId: string) => void;
  onSampleRevisionIdsChange: (ids: string[]) => void;
  onExperimentRevisionIdChange: (id: string | null) => void;
}

export function ExperimentRevisionPanel({
  handoff,
  projectId,
  sampleRevisionIds,
  experimentRevisionId,
  onDirectProjectChange,
  onSampleRevisionIdsChange,
  onExperimentRevisionIdChange,
}: ExperimentRevisionPanelProps) {
  const directEntry = handoff === null;
  const projects = useQuery({
    queryKey: ['experiment-entry-projects'],
    queryFn: loadEntryProjects,
    enabled: directEntry,
    retry: false,
  });
  const samples = useQuery({
    queryKey: ['experiment-entry-samples', projectId],
    queryFn: () => loadEntrySamples(projectId ?? ''),
    enabled: directEntry && projectId !== null,
    retry: false,
  });
  const experiments = useQuery({
    queryKey: ['catalog-experiments', projectId],
    queryFn: () => loadSelectableExperiments(projectId ?? ''),
    enabled: projectId !== null,
    retry: false,
  });
  const projectChoices = projects.data ?? [];
  const sampleChoices = samples.data ?? [];
  const choices = experiments.data ?? [];
  const selected =
    choices.find((item) => item.experimentRevisionId === experimentRevisionId) ??
    null;

  useEffect(() => {
    if (!directEntry || projectId !== null || projectChoices.length === 0) return;
    onDirectProjectChange?.(projectChoices[0].projectId);
  }, [directEntry, onDirectProjectChange, projectChoices, projectId]);

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

  function toggleSample(sampleRevisionId: string, checked: boolean) {
    if (checked) {
      if (sampleRevisionIds.includes(sampleRevisionId)) return;
      onSampleRevisionIdsChange([...sampleRevisionIds, sampleRevisionId]);
      return;
    }
    onSampleRevisionIdsChange(
      sampleRevisionIds.filter((id) => id !== sampleRevisionId),
    );
  }

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
          {directEntry
            ? 'Choose a project and sample revisions, then select or omit an experiment. Choosing a revision does not bind it.'
            : "Carried from the saved draft. Choosing a revision does not bind it. A binding appears only after this page's validation returns one."}
        </p>
      </div>
      {handoff !== null && (
        <>
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
        </>
      )}
      {directEntry && (
        <>
          <label className="block text-sm font-medium" htmlFor="experiment-project-choice">
            Project
            <select
              id="experiment-project-choice"
              className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
              value={projectId ?? ''}
              disabled={projects.isLoading}
              onChange={(event) => {
                if (event.target.value !== '') onDirectProjectChange?.(event.target.value);
              }}
            >
              <option value="">Select a project</option>
              {projectChoices.map((project) => (
                <option key={project.projectId} value={project.projectId}>
                  {project.displayName}
                </option>
              ))}
            </select>
          </label>
          {projects.isLoading && <p>Loading projects.</p>}
          {projects.isError && <p>Active projects could not be loaded.</p>}
          <fieldset className="space-y-2">
            <legend className="text-sm font-medium">Project samples</legend>
            {projectId === null && (
              <p className="text-xs text-[var(--color-text-muted)]">
                Select a project to list its sample revisions.
              </p>
            )}
            {samples.isLoading && <p>Loading project samples.</p>}
            {samples.isError && <p>Project samples could not be loaded.</p>}
            {samples.isSuccess && sampleChoices.length === 0 && (
              <p>This project has no sample revisions yet.</p>
            )}
            <ul className="space-y-1" aria-label="Project sample revisions">
              {sampleChoices.map((sample) => (
                <li key={sample.sampleRevisionId}>
                  <label className="flex min-h-11 items-center gap-2 text-sm">
                    <input
                      type="checkbox"
                      checked={sampleRevisionIds.includes(sample.sampleRevisionId)}
                      onChange={(event) =>
                        toggleSample(sample.sampleRevisionId, event.target.checked)
                      }
                    />
                    {sample.displayName}
                  </label>
                </li>
              ))}
            </ul>
          </fieldset>
        </>
      )}
      {(handoff !== null || sampleRevisionIds.length > 0) && (
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
      )}
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

function readEntryProject(project: {
  project_id: string;
  display_name: string;
}): EntryProject | null {
  if (
    !PROJECT_ID.test(project.project_id) ||
    project.project_id === `prj_${'0'.repeat(32)}`
  ) {
    return null;
  }
  const displayName = project.display_name.trim();
  return {
    projectId: project.project_id,
    displayName: displayName.length > 0 ? displayName : project.project_id,
  };
}

function readEntrySample(sample: {
  sample_revision_id: string;
  display_name: string;
}): EntrySample | null {
  if (!SAMPLE_REVISION_ID.test(sample.sample_revision_id)) return null;
  const displayName = sample.display_name.trim();
  return {
    sampleRevisionId: sample.sample_revision_id,
    displayName: displayName.length > 0 ? displayName : sample.sample_revision_id,
  };
}

async function loadEntryProjects(): Promise<EntryProject[]> {
  const seenCursors = new Set<string>();
  const seenProjects = new Set<string>();
  const projects: EntryProject[] = [];
  let cursor: string | null = null;
  for (let page = 0; page < 20; page += 1) {
    const response = await authoringListProjects({
      limit: 50,
      ...(cursor === null ? {} : { cursor }),
    });
    if (!Array.isArray(response.projects)) {
      throw new Error('The project list could not be verified.');
    }
    for (const project of response.projects) {
      const parsed = readEntryProject(project);
      if (parsed === null || seenProjects.has(parsed.projectId)) continue;
      seenProjects.add(parsed.projectId);
      projects.push(parsed);
    }
    const next = response.next_cursor;
    if (typeof next !== 'string' || next.length === 0) return projects;
    if (seenCursors.has(next)) {
      throw new Error('The project list repeated a page.');
    }
    seenCursors.add(next);
    cursor = next;
  }
  throw new Error('The project list did not finish.');
}

async function loadEntrySamples(projectId: string): Promise<EntrySample[]> {
  const seenCursors = new Set<string>();
  const seenSamples = new Set<string>();
  const samples: EntrySample[] = [];
  let cursor: string | null = null;
  for (let page = 0; page < 20; page += 1) {
    const response = await authoringListProjectSamples(projectId, {
      limit: 50,
      ...(cursor === null ? {} : { cursor }),
    });
    if (response.project_id !== projectId || !Array.isArray(response.samples)) {
      throw new Error('The sample list could not be verified.');
    }
    for (const sample of response.samples) {
      const parsed = readEntrySample(sample);
      if (parsed === null || seenSamples.has(parsed.sampleRevisionId)) continue;
      seenSamples.add(parsed.sampleRevisionId);
      samples.push(parsed);
    }
    const next = response.next_cursor;
    if (typeof next !== 'string' || next.length === 0) return samples;
    if (seenCursors.has(next)) {
      throw new Error('The sample list repeated a page.');
    }
    seenCursors.add(next);
    cursor = next;
  }
  throw new Error('The sample list did not finish.');
}
