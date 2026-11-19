import type {
  InputFileRevisionSelectionRequest,
  ValidationRequestConfig,
  ValidationRequestOptions,
  ValidationRequestSamplesAnyOfItem,
} from '../../api/generated/models';
import { isJsonObject, isJsonValue, isPlainObject } from '../input-workbench/jsonSafety';

export const EXPERIMENT_VALIDATION_HANDOFF_KEY = 'experimentValidationHandoff';

const LEGACY_PROJECT_ID = 'prj_00000000000000000000000000000000';
const PROJECT_ID = /^prj_[0-9a-f]{32}$/;
const SAMPLE_REVISION_ID = /^smpr_[0-9a-f]{32}$/;
const REFERENCE_PROFILE_REVISION_ID = /^refpr_[0-9a-f]{32}$/;
const INPUT_FILE_REVISION_ID = /^inpfr_[0-9a-f]{32}$/;
const INPUT_USE_KEY = /^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$/;
const EXPERIMENT_REVISION_ID = /^expr_[0-9a-f]{32}$/;
const MAX_SAMPLE_REVISIONS = 1000;
const MAX_INPUT_SELECTIONS = 256;

export interface ExperimentValidationHandoff {
  workflowId: string;
  projectId: string;
  config: ValidationRequestConfig;
  options: ValidationRequestOptions;
  samples: ValidationRequestSamplesAnyOfItem[] | null;
  sampleRevisionIds: string[];
  inputSelections: InputFileRevisionSelectionRequest[];
  referenceProfileRevisionId: string | null;
  notes: readonly string[];
}

export interface DraftCarrySource {
  workflowId: string;
  projectId: unknown;
  config: unknown;
  options: unknown;
  samples: unknown;
  sampleRevisionIds: unknown;
  inputSelections: unknown;
  referenceProfileRevisionId: unknown;
  workflowSamplesPresent: boolean;
}

export interface DraftCarryResult {
  handoff: ExperimentValidationHandoff | null;
  notes: readonly string[];
}

export function sampleRevisionSetsEqual(
  left: readonly string[],
  right: readonly string[],
): boolean {
  if (left.length === 0 || right.length === 0 || left.length !== right.length) {
    return false;
  }
  if (new Set(left).size !== left.length || new Set(right).size !== right.length) {
    return false;
  }
  const rightIds = new Set(right);
  return left.every((id) => rightIds.has(id));
}

export function readDraftCarry(source: DraftCarrySource): DraftCarryResult {
  const notes = [
    'Opening the validation page does not bind an experiment.',
    'A previous draft validation result is not an experiment binding.',
  ];
  if (source.workflowSamplesPresent) {
    notes.push(
      'The workflow sample plan stays on the saved draft. This page submits the sample revision ids.',
    );
  }

  const config = readObject(source.config);
  const options = source.options === undefined ? {} : readObject(source.options);
  const projectId = readProjectId(source.projectId);
  const sampleRevisionIds = readSampleRevisionIds(source.sampleRevisionIds);
  if (config === null) {
    notes.push('The saved configuration is not a JSON object, so it was not copied.');
  }
  if (options === null) {
    notes.push('The saved options are not a JSON object, so they were not copied.');
  }
  if (projectId === null) {
    notes.push('The draft has no active project that this page can submit.');
  }
  if (sampleRevisionIds === null) {
    notes.push(
      'The saved sample revisions are empty, duplicated, or not sample revision ids.',
    );
  }
  if (
    config === null ||
    options === null ||
    projectId === null ||
    sampleRevisionIds === null
  ) {
    return { handoff: null, notes };
  }

  const samples = readSampleSheet(source.samples, notes);
  const inputSelections = readInputSelections(source.inputSelections, notes);
  const referenceProfileRevisionId = readReference(source.referenceProfileRevisionId, notes);

  return {
    handoff: {
      workflowId: source.workflowId,
      projectId,
      config,
      options,
      samples,
      sampleRevisionIds,
      inputSelections,
      referenceProfileRevisionId,
      notes,
    },
    notes,
  };
}

export function readCarriedHandoff(
  state: unknown,
  workflowId: string,
): ExperimentValidationHandoff | null {
  if (!isPlainObject(state)) return null;
  const candidate = state[EXPERIMENT_VALIDATION_HANDOFF_KEY];
  if (!isPlainObject(candidate) || candidate.workflowId !== workflowId) return null;
  const checked = readDraftCarry({
    workflowId,
    projectId: candidate.projectId,
    config: candidate.config,
    options: candidate.options,
    samples: candidate.samples,
    sampleRevisionIds: candidate.sampleRevisionIds,
    inputSelections: candidate.inputSelections,
    referenceProfileRevisionId: candidate.referenceProfileRevisionId,
    workflowSamplesPresent: false,
  });
  if (checked.handoff === null) return null;
  return {
    ...checked.handoff,
    notes: Array.isArray(candidate.notes)
      ? candidate.notes.filter((note): note is string => typeof note === 'string')
      : checked.handoff.notes,
  };
}

export function readExperimentBinding(value: unknown):
  | { status: 'unresolved' }
  | { status: 'experiment_snapshot'; experimentRevisionId: string }
  | null {
  if (!isPlainObject(value)) return null;
  if (value.status === 'unresolved') return { status: 'unresolved' };
  if (
    value.status === 'experiment_snapshot' &&
    typeof value.experiment_revision_id === 'string' &&
    EXPERIMENT_REVISION_ID.test(value.experiment_revision_id)
  ) {
    return {
      status: 'experiment_snapshot',
      experimentRevisionId: value.experiment_revision_id,
    };
  }
  return null;
}

function readObject(value: unknown): ValidationRequestConfig | null {
  if (!isJsonObject(value)) return null;
  return JSON.parse(JSON.stringify(value)) as ValidationRequestConfig;
}

function readProjectId(value: unknown): string | null {
  if (typeof value !== 'string' || !PROJECT_ID.test(value) || value === LEGACY_PROJECT_ID) {
    return null;
  }
  return value;
}

function readSampleRevisionIds(value: unknown): string[] | null {
  if (!Array.isArray(value) || value.length === 0 || value.length > MAX_SAMPLE_REVISIONS) {
    return null;
  }
  if (!value.every((id): id is string => typeof id === 'string' && SAMPLE_REVISION_ID.test(id))) {
    return null;
  }
  if (new Set(value).size !== value.length) return null;
  return [...value];
}

function readSampleSheet(
  value: unknown,
  notes: string[],
): ValidationRequestSamplesAnyOfItem[] | null {
  if (value === null || value === undefined) {
    notes.push(
      'The draft has no sample sheet. Sample revision ids are copied separately.',
    );
    return null;
  }
  if (typeof value === 'string') {
    notes.push('A sample path was not copied. Sample revision ids are copied separately.');
    return null;
  }
  if (!Array.isArray(value) || value.length === 0 || value.length > MAX_SAMPLE_REVISIONS) {
    notes.push('The sample sheet was not copied. Sample revision ids are copied separately.');
    return null;
  }
  const rows: ValidationRequestSamplesAnyOfItem[] = [];
  for (const item of value) {
    if (!isPlainObject(item) || !isJsonValue(item)) {
      notes.push('The sample sheet was not copied. Sample revision ids are copied separately.');
      return null;
    }
    const row: Record<string, string> = {};
    for (const [key, cell] of Object.entries(item)) {
      if (typeof cell !== 'string' || key.length === 0 || /[\u0000\t\n\r]/.test(key + cell)) {
        notes.push('The sample sheet was not copied. Sample revision ids are copied separately.');
        return null;
      }
      row[key] = cell;
    }
    rows.push(row);
  }
  notes.push('The sample sheet was copied into this page.');
  return rows;
}

function readReference(value: unknown, notes: string[]): string | null {
  if (value === null || value === undefined) return null;
  if (typeof value === 'string' && REFERENCE_PROFILE_REVISION_ID.test(value)) return value;
  notes.push('The saved reference profile was not copied.');
  return null;
}

function readInputSelections(
  value: unknown,
  notes: string[],
): InputFileRevisionSelectionRequest[] {
  if (value === undefined || value === null) return [];
  if (!Array.isArray(value) || value.length > MAX_INPUT_SELECTIONS) {
    notes.push('Input revision selections were not copied.');
    return [];
  }
  const selections: InputFileRevisionSelectionRequest[] = [];
  const coordinates = new Set<string>();
  for (const item of value) {
    if (!isPlainObject(item)) {
      notes.push('Input revision selections were not copied.');
      return [];
    }
    const key = item.input_use_key;
    const occurrence = item.occurrence;
    const ids = item.input_file_revision_ids;
    if (
      typeof key !== 'string' ||
      !INPUT_USE_KEY.test(key) ||
      typeof occurrence !== 'number' ||
      !Number.isInteger(occurrence) ||
      occurrence < 0 ||
      occurrence > 1_000_000 ||
      !Array.isArray(ids) ||
      ids.length === 0 ||
      ids.length > 256 ||
      !ids.every((id): id is string => typeof id === 'string' && INPUT_FILE_REVISION_ID.test(id)) ||
      new Set(ids).size !== ids.length
    ) {
      notes.push('Input revision selections were not copied.');
      return [];
    }
    const coordinate = `${key}\u0000${occurrence}`;
    if (coordinates.has(coordinate)) {
      notes.push('Input revision selections were not copied.');
      return [];
    }
    coordinates.add(coordinate);
    selections.push({
      input_use_key: key,
      occurrence,
      input_file_revision_ids: [...ids],
    });
  }
  return selections;
}
