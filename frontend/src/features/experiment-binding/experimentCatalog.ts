import { catalogListExperiments } from '../../api/generated/catalog/catalog';
import type { CatalogExperimentResponse } from '../../api/generated/models';

const EXPERIMENT_REVISION_ID = /^expr_[0-9a-f]{32}$/;
const SAMPLE_REVISION_ID = /^smpr_[0-9a-f]{32}$/;

export interface SelectableExperiment {
  experimentRevisionId: string;
  displayName: string;
  sampleRevisionIds: string[];
}

export async function loadSelectableExperiments(
  projectId: string,
): Promise<SelectableExperiment[]> {
  const seenCursors = new Set<string>();
  const experiments: SelectableExperiment[] = [];
  let cursor: string | null = null;
  for (let page = 0; page < 20; page += 1) {
    const response = await catalogListExperiments(projectId, {
      limit: 50,
      include_archived: false,
      ...(cursor === null ? {} : { cursor }),
    });
    if (response.project_id !== projectId || !Array.isArray(response.experiments)) {
      throw new Error('The experiment list could not be verified.');
    }
    for (const item of response.experiments) {
      const parsed = readSelectableExperiment(item);
      if (parsed !== null) experiments.push(parsed);
    }
    const next = response.next_cursor;
    if (typeof next !== 'string' || next.length === 0) return experiments;
    if (seenCursors.has(next)) {
      throw new Error('The experiment list repeated a page.');
    }
    seenCursors.add(next);
    cursor = next;
  }
  throw new Error('The experiment list did not finish.');
}

function readSelectableExperiment(
  item: CatalogExperimentResponse,
): SelectableExperiment | null {
  if (!EXPERIMENT_REVISION_ID.test(item.experiment_revision_id)) return null;
  if (!Array.isArray(item.samples)) return null;
  const sampleRevisionIds: string[] = [];
  for (const sample of item.samples) {
    if (!SAMPLE_REVISION_ID.test(sample.sample_revision_id)) return null;
    sampleRevisionIds.push(sample.sample_revision_id);
  }
  const displayName = item.display_name.trim();
  return {
    experimentRevisionId: item.experiment_revision_id,
    displayName: displayName.length > 0 ? displayName : item.experiment_revision_id,
    sampleRevisionIds,
  };
}
