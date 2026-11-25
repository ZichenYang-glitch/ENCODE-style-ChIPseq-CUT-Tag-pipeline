/**
 * Read sealed track versions for one experiment.
 *
 * OpenAPI 0.4.1 does not generate catalogListArchivedTracks or
 * catalogGetArchivedTrack. These helpers keep those names and call the
 * routes in docs/development/experiment-archive-track-history.md through
 * the shared fetcher. When central regeneration adds the generated
 * exports, delegate these two bodies to them and keep the parsers.
 *
 * Expected generated signatures:
 * catalogListArchivedTracks(projectId, experimentId, params?, options?)
 * catalogGetArchivedTrack(projectId, experimentId, archiveVersionId, options?)
 */
import { archivedTrackCollectionPath, archivedTrackResourcePath } from './archivedTrackPaths';
import { ApiError, fetcher } from './fetcher';

const PROJECT_ID = /^prj_[0-9a-f]{32}$/;
const EXPERIMENT_ID = /^exp_[0-9a-f]{32}$/;
const EXPERIMENT_REVISION_ID = /^expr_[0-9a-f]{32}$/;
const ARCHIVE_VERSION_ID = /^[A-Za-z][A-Za-z0-9_.-]{0,127}$/;
const RUN_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$/;
const ARTIFACT_ID = /^[A-Za-z][A-Za-z0-9_.-]{0,127}$/;
const GENERATION = /^artifactgen-[0-9a-f]{64}$/;
const ARTIFACT_REVISION = /^artifactrev-[0-9a-f]{64}$/;
const RESULT_ID = /^rreg_[0-9a-f]{32}$/;
const RESULT_REVISION_ID = /^rrev_[0-9a-f]{32}$/;
const PUBLIC_LABEL = /^[^\u0000-\u001f\u007f-\u009f\\/]{1,255}$/u;
const ASSEMBLY = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$/;

const TRACK_FIELDS = new Set([
  'archive_version_id',
  'project_id',
  'experiment_id',
  'experiment_revision_id',
  'track_key',
  'version_number',
  'display_name',
  'assembly',
  'track_type',
  'source',
  'has_index',
]);

const RUN_FIELDS = new Set([
  'kind',
  'run_id',
  'artifact_id',
  'generation',
  'revision',
  'index_artifact_id',
  'index_generation',
  'index_revision',
]);
const RUN_INDEX_FIELDS = ['index_artifact_id', 'index_generation', 'index_revision'] as const;

const REGISTERED_FIELDS = new Set([
  'kind',
  'result_id',
  'result_revision_id',
  'index_result_id',
  'index_result_revision_id',
]);

export type ArchiveTrackType = 'bigwig' | 'bed' | 'longrange';

export type RunArchiveSource = {
  kind: 'run';
  run_id: string;
  artifact_id: string;
  generation: string;
  revision: string;
  index_artifact_id?: string;
  index_generation?: string;
  index_revision?: string;
};

export type RegisteredArchiveSource = {
  kind: 'registered';
  result_id: string;
  result_revision_id: string;
  index_result_id?: string;
  index_result_revision_id?: string;
};

export type ArchiveSource = RunArchiveSource | RegisteredArchiveSource;

export interface ArchivedTrack {
  archive_version_id: string;
  project_id: string;
  experiment_id: string;
  experiment_revision_id: string;
  track_key: string;
  version_number: number;
  display_name: string;
  assembly: string;
  track_type: ArchiveTrackType;
  source: ArchiveSource;
  has_index: boolean;
}

export interface ArchivedTrackList {
  items: ArchivedTrack[];
  next_offset: number | null;
}

export interface ArchivedTrackScope {
  projectId: string;
  experimentId: string;
  archiveVersionId?: string;
  offset?: number;
}

function rejected(status = 422): never {
  throw new ApiError(status, 'archive_invalid', 'Archived track descriptor was rejected.');
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function exactKeys(value: Record<string, unknown>, allowed: Set<string>): boolean {
  return Object.keys(value).every((key) => allowed.has(key));
}

function text(value: unknown, pattern: RegExp): string | null {
  if (typeof value !== 'string' || !pattern.test(value)) return null;
  return value;
}

function optionalText(
  value: Record<string, unknown>,
  field: string,
  pattern: RegExp,
): string | null | undefined {
  if (value[field] === undefined) return undefined;
  return text(value[field], pattern);
}

function parseSource(value: unknown): ArchiveSource | null {
  if (!isRecord(value)) return null;
  const kind = value.kind;
  if (kind === 'run') {
    if (!exactKeys(value, RUN_FIELDS)) return null;
    const runId = text(value.run_id, RUN_ID);
    const artifactId = text(value.artifact_id, ARTIFACT_ID);
    const generation = text(value.generation, GENERATION);
    const revision = text(value.revision, ARTIFACT_REVISION);
    if (runId === null || artifactId === null || generation === null || revision === null) {
      return null;
    }
    const indexArtifactId = optionalText(value, 'index_artifact_id', ARTIFACT_ID);
    const indexGeneration = optionalText(value, 'index_generation', GENERATION);
    const indexRevision = optionalText(value, 'index_revision', ARTIFACT_REVISION);
    const indexValues = [indexArtifactId, indexGeneration, indexRevision];
    if (indexValues.some((item) => item === null)) return null;
    const indexCount = indexValues.filter((item) => item !== undefined).length;
    if (indexCount !== 0 && indexCount !== RUN_INDEX_FIELDS.length) return null;
    return {
      kind: 'run',
      run_id: runId,
      artifact_id: artifactId,
      generation,
      revision,
      ...(indexCount === 0
        ? {}
        : {
            index_artifact_id: indexArtifactId as string,
            index_generation: indexGeneration as string,
            index_revision: indexRevision as string,
          }),
    };
  }
  if (kind === 'registered') {
    if (!exactKeys(value, REGISTERED_FIELDS)) return null;
    const resultId = text(value.result_id, RESULT_ID);
    const resultRevisionId = text(value.result_revision_id, RESULT_REVISION_ID);
    if (resultId === null || resultRevisionId === null) return null;
    const indexResultId = optionalText(value, 'index_result_id', RESULT_ID);
    const indexResultRevisionId = optionalText(value, 'index_result_revision_id', RESULT_REVISION_ID);
    if (indexResultId === null || indexResultRevisionId === null) return null;
    if ((indexResultId === undefined) !== (indexResultRevisionId === undefined)) return null;
    return {
      kind: 'registered',
      result_id: resultId,
      result_revision_id: resultRevisionId,
      ...(indexResultId === undefined
        ? {}
        : {
            index_result_id: indexResultId,
            index_result_revision_id: indexResultRevisionId as string,
          }),
    };
  }
  return null;
}

export function parseArchivedTrack(value: unknown, scope: ArchivedTrackScope): ArchivedTrack | null {
  if (!isRecord(value) || !exactKeys(value, TRACK_FIELDS)) return null;
  const archiveVersionId = text(value.archive_version_id, ARCHIVE_VERSION_ID);
  const projectId = text(value.project_id, PROJECT_ID);
  const experimentId = text(value.experiment_id, EXPERIMENT_ID);
  const experimentRevisionId = text(value.experiment_revision_id, EXPERIMENT_REVISION_ID);
  const trackKey = text(value.track_key, PUBLIC_LABEL);
  const displayName = text(value.display_name, PUBLIC_LABEL);
  const assembly = text(value.assembly, ASSEMBLY);
  const trackType = value.track_type;
  const source = parseSource(value.source);
  if (
    archiveVersionId === null ||
    projectId === null ||
    experimentId === null ||
    experimentRevisionId === null ||
    trackKey === null ||
    displayName === null ||
    assembly === null ||
    source === null ||
    (trackType !== 'bigwig' && trackType !== 'bed' && trackType !== 'longrange') ||
    typeof value.has_index !== 'boolean' ||
    !Number.isSafeInteger(value.version_number) ||
    (value.version_number as number) < 1
  ) {
    return null;
  }
  if (projectId !== scope.projectId || experimentId !== scope.experimentId) return null;
  if (scope.archiveVersionId !== undefined && archiveVersionId !== scope.archiveVersionId) {
    return null;
  }
  return {
    archive_version_id: archiveVersionId,
    project_id: projectId,
    experiment_id: experimentId,
    experiment_revision_id: experimentRevisionId,
    track_key: trackKey,
    version_number: value.version_number as number,
    display_name: displayName,
    assembly,
    track_type: trackType,
    source,
    has_index: value.has_index,
  };
}

export function parseArchivedTrackList(value: unknown, scope: ArchivedTrackScope): ArchivedTrackList | null {
  if (!isRecord(value) || !exactKeys(value, new Set(['items', 'next_offset']))) return null;
  if (!Array.isArray(value.items) || value.items.length > 100) return null;
  const next = value.next_offset;
  if (next !== null && (!Number.isSafeInteger(next) || (next as number) < 0)) return null;
  if (scope.offset !== undefined && next !== null && (next as number) <= scope.offset) return null;
  const items: ArchivedTrack[] = [];
  const seen = new Set<string>();
  for (const item of value.items) {
    const track = parseArchivedTrack(item, scope);
    if (track === null || seen.has(track.archive_version_id)) return null;
    seen.add(track.archive_version_id);
    items.push(track);
  }
  return { items, next_offset: next as number | null };
}

function requireScope(projectId: string, experimentId: string, archiveVersionId?: string): void {
  if (!PROJECT_ID.test(projectId) || !EXPERIMENT_ID.test(experimentId)) rejected();
  if (archiveVersionId !== undefined && !ARCHIVE_VERSION_ID.test(archiveVersionId)) rejected();
}

function listQuery(params?: { limit?: number; offset?: number }): { suffix: string; offset: number } {
  const search = new URLSearchParams();
  let offset = 0;
  if (params?.limit !== undefined) {
    if (!Number.isSafeInteger(params.limit) || params.limit < 1 || params.limit > 100) rejected(400);
    search.set('limit', String(params.limit));
  }
  if (params?.offset !== undefined) {
    if (!Number.isSafeInteger(params.offset) || params.offset < 0) rejected(400);
    offset = params.offset;
    search.set('offset', String(params.offset));
  }
  const query = search.toString();
  return { suffix: query === '' ? '' : `?${query}`, offset };
}

export const catalogListArchivedTracks = async (
  projectId: string,
  experimentId: string,
  params?: { limit?: number; offset?: number },
  options?: RequestInit,
): Promise<ArchivedTrackList> => {
  requireScope(projectId, experimentId);
  const query = listQuery(params);
  const body = await fetcher<unknown>(
    `${archivedTrackCollectionPath(projectId, experimentId)}${query.suffix}`,
    { ...options, method: 'GET' },
  );
  const parsed = parseArchivedTrackList(body, {
    projectId,
    experimentId,
    offset: query.offset,
  });
  if (parsed === null) rejected();
  return parsed;
};

export const catalogGetArchivedTrack = async (
  projectId: string,
  experimentId: string,
  archiveVersionId: string,
  options?: RequestInit,
): Promise<ArchivedTrack> => {
  requireScope(projectId, experimentId, archiveVersionId);
  const body = await fetcher<unknown>(
    archivedTrackResourcePath(projectId, experimentId, archiveVersionId),
    { ...options, method: 'GET' },
  );
  const parsed = parseArchivedTrack(body, { projectId, experimentId, archiveVersionId });
  if (parsed === null) rejected();
  return parsed;
};
