import type {
  ArtifactPublicationListResponse,
  ArtifactPublicationResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';

export const ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE = 50;

const LEGACY_PROJECT_ID = 'prj_00000000000000000000000000000000';
const PROJECT_ID = /^prj_[0-9a-f]{32}$/;
const RUN_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$/;
const WORKFLOW_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,254}$/;
const LOGICAL_ID = /^[A-Za-z][A-Za-z0-9_.-]{0,127}$/;
const SAMPLE_ID = /^smp_[0-9a-f]{32}$/;
const SAMPLE_REVISION_ID = /^smpr_[0-9a-f]{32}$/;
const ARTIFACT_GENERATION = /^artifactgen-[0-9a-f]{64}$/;
const ARTIFACT_REVISION = /^artifactrev-[0-9a-f]{64}$/;
const ARTIFACT_PUBLICATION_CURSOR = /^artifactpubcur_[A-Za-z0-9_-]+$/;
const TIMEZONE_TIMESTAMP =
  /^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2})(?::(\d{2})(?:\.(\d{1,6}))?)?(?:(Z)|([+-])(\d{2}):(\d{2}))$/;

const LIST_QUERY_PARAMETERS = new Set([
  'project_id',
  'run_id',
  'workflow_id',
  'output_type',
  'associated_run_sample_revision_id',
  'published_from',
  'published_before',
  'current_only',
  'after',
  'limit',
]);

const PUBLICATION_KEYS = new Set([
  'run_id',
  'project_id',
  'workflow_id',
  'artifact_id',
  'output_type',
  'artifact_generation',
  'artifact_revision',
  'published_at',
  'current_artifact_generation',
  'generation_status',
  'run_sample_binding',
]);

export interface ArtifactPublicationFilters {
  projectId: string | null;
  runId: string | null;
  workflowId: string | null;
  outputType: string | null;
  associatedRunSampleRevisionId: string | null;
  publishedFrom: string | null;
  publishedBefore: string | null;
  currentOnly: boolean;
  after: string | null;
  limit: number;
}

export type ArtifactPublicationFilterResult =
  | { ok: true; filters: ArtifactPublicationFilters }
  | { ok: false; filters: null };

export type ArtifactPublicationDetailIdentityResult =
  | { ok: true; runId: string; artifactId: string; generation: string }
  | { ok: false };

export class ArtifactPublicationProtocolError extends Error {
  constructor() {
    super('artifact publication response is invalid');
    this.name = 'ArtifactPublicationProtocolError';
  }
}

export function parseArtifactPublicationFilters(
  searchParams: URLSearchParams,
): ArtifactPublicationFilterResult {
  for (const key of searchParams.keys()) {
    if (!LIST_QUERY_PARAMETERS.has(key)) return { ok: false, filters: null };
  }
  for (const key of LIST_QUERY_PARAMETERS) {
    if (searchParams.getAll(key).length > 1) {
      return { ok: false, filters: null };
    }
  }

  const projectId = optionalToken(searchParams, 'project_id', PROJECT_ID);
  const runId = optionalToken(searchParams, 'run_id', RUN_ID);
  const workflowId = optionalToken(searchParams, 'workflow_id', WORKFLOW_ID);
  const outputType = optionalToken(searchParams, 'output_type', LOGICAL_ID);
  const associatedRunSampleRevisionId = optionalToken(
    searchParams,
    'associated_run_sample_revision_id',
    SAMPLE_REVISION_ID,
  );
  const publishedFrom = optionalTimestamp(searchParams, 'published_from');
  const publishedBefore = optionalTimestamp(searchParams, 'published_before');
  const after = optionalToken(
    searchParams,
    'after',
    ARTIFACT_PUBLICATION_CURSOR,
    1024,
  );
  if (
    projectId === false ||
    runId === false ||
    workflowId === false ||
    outputType === false ||
    associatedRunSampleRevisionId === false ||
    publishedFrom === false ||
    publishedBefore === false ||
    after === false
  ) {
    return { ok: false, filters: null };
  }

  const currentOnlyValue = searchParams.get('current_only');
  if (
    currentOnlyValue !== null &&
    currentOnlyValue !== 'true' &&
    currentOnlyValue !== 'false'
  ) {
    return { ok: false, filters: null };
  }
  const limitValue = searchParams.get('limit');
  const limit = limitValue === null ? ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE : Number(limitValue);
  if (
    !Number.isSafeInteger(limit) ||
    limit < 1 ||
    limit > 100 ||
    (limitValue !== null && String(limit) !== limitValue)
  ) {
    return { ok: false, filters: null };
  }
  if (
    publishedFrom !== null &&
    publishedBefore !== null &&
    publishedBefore <= publishedFrom
  ) {
    return { ok: false, filters: null };
  }

  return {
    ok: true,
    filters: {
      projectId,
      runId,
      workflowId,
      outputType,
      associatedRunSampleRevisionId,
      publishedFrom,
      publishedBefore,
      currentOnly: currentOnlyValue !== 'false',
      after,
      limit,
    },
  };
}

export function parseArtifactPublicationDetailIdentity(
  runId: string | undefined,
  artifactId: string | undefined,
  searchParams: URLSearchParams,
): ArtifactPublicationDetailIdentityResult {
  const keys = [...searchParams.keys()];
  const generations = searchParams.getAll('generation');
  const generation = generations[0];
  if (
    keys.some((key) => key !== 'generation') ||
    generations.length !== 1 ||
    runId === undefined ||
    artifactId === undefined ||
    !RUN_ID.test(runId) ||
    !LOGICAL_ID.test(artifactId) ||
    generation === undefined ||
    !ARTIFACT_GENERATION.test(generation)
  ) {
    return { ok: false };
  }
  return { ok: true, runId, artifactId, generation };
}

export function artifactPublicationQueryKey(filters: ArtifactPublicationFilters) {
  return [
    'artifact-publications',
    filters.projectId,
    filters.runId,
    filters.workflowId,
    filters.outputType,
    filters.associatedRunSampleRevisionId,
    filters.publishedFrom,
    filters.publishedBefore,
    filters.currentOnly,
    filters.after,
    filters.limit,
  ] as const;
}

export function artifactPublicationDetailQueryKey(
  runId: string,
  artifactId: string,
  generation: string,
) {
  return ['artifact-publication', runId, artifactId, generation] as const;
}

export function validateArtifactPublicationListResponse(
  response: ArtifactPublicationListResponse,
  requestedCursor: string | undefined,
): ArtifactPublicationListResponse {
  if (
    !hasExactKeys(response, ['ok', 'publications', 'next_cursor', 'issues']) ||
    response.ok !== true ||
    !Array.isArray(response.publications) ||
    !Array.isArray(response.issues) ||
    response.issues.length !== 0 ||
    !isSafeArtifactPublicationCursor(response.next_cursor) ||
    (requestedCursor !== undefined && response.next_cursor === requestedCursor)
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  const publications = response.publications.map(projectArtifactPublication);
  const identities = publications.map(artifactPublicationIdentity);
  if (new Set(identities).size !== identities.length) {
    throw new ArtifactPublicationProtocolError();
  }
  return {
    ok: true,
    publications,
    next_cursor: response.next_cursor,
    issues: [],
  };
}

export function validateArtifactPublicationDetailResponse(
  response: unknown,
  runId: string,
  artifactId: string,
  generation: string,
): ArtifactPublicationResponse {
  if (
    !hasExactKeys(response, ['ok', 'publication', 'issues']) ||
    response.ok !== true ||
    !Array.isArray(response.issues) ||
    response.issues.length !== 0 ||
    response.publication === null
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  const publication = projectArtifactPublication(response.publication);
  if (
    publication.run_id !== runId ||
    publication.artifact_id !== artifactId ||
    publication.artifact_generation !== generation
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  return publication;
}

export function flattenArtifactPublicationPages(
  pages: ArtifactPublicationListResponse[],
): ArtifactPublicationResponse[] {
  const seen = new Set<string>();
  const publications: ArtifactPublicationResponse[] = [];
  for (const page of pages) {
    const pageIdentities = page.publications.map(artifactPublicationIdentity);
    if (
      new Set(pageIdentities).size !== pageIdentities.length ||
      pageIdentities.some((identity) => seen.has(identity))
    ) break;
    pageIdentities.forEach((identity) => seen.add(identity));
    publications.push(...page.publications);
  }
  return publications;
}

export function artifactPublicationPaginationAnomaly(
  pages: ArtifactPublicationListResponse[],
  pageParams: unknown[],
): boolean {
  const usedCursors = new Set<string>();
  const identities = new Set<string>();
  for (const pageParam of pageParams) {
    if (typeof pageParam !== 'string') continue;
    if (!isSafeArtifactPublicationCursor(pageParam) || usedCursors.has(pageParam)) {
      return true;
    }
    usedCursors.add(pageParam);
  }
  for (const [pageIndex, page] of pages.entries()) {
    for (const publication of page.publications) {
      const identity = artifactPublicationIdentity(publication);
      if (identities.has(identity)) return true;
      identities.add(identity);
    }
    const next = page.next_cursor;
    if (next !== null) {
      const current = pageParams[pageIndex];
      const following = pageParams[pageIndex + 1];
      if (next === current) return true;
      if (following !== undefined) {
        if (next !== following) return true;
      } else if (usedCursors.has(next)) {
        return true;
      }
    }
  }
  return false;
}

export function safeArtifactPublicationNextCursor(
  page: ArtifactPublicationListResponse,
  allPages: ArtifactPublicationListResponse[],
  pageParams: unknown[],
): string | undefined {
  const cursor = page.next_cursor;
  if (cursor === null || !isSafeArtifactPublicationCursor(cursor)) return undefined;
  if (pageParams.includes(cursor)) return undefined;
  if (allPages.slice(0, -1).some((item) => item.next_cursor === cursor)) {
    return undefined;
  }
  return cursor;
}

export function isArtifactPublicationRecoveryError(error: unknown): boolean {
  return (
    error instanceof ArtifactPublicationProtocolError ||
    (error instanceof ApiError &&
      (error.code === 'ARTIFACT_PUBLICATION_CURSOR_INVALID' ||
        error.code === 'ARTIFACT_PUBLICATION_CURSOR_NOT_FOUND' ||
        error.code === 'ARTIFACT_PUBLICATION_DATA_INVALID'))
  );
}

export function derivedArtifactGenerationStatus(
  publication: ArtifactPublicationResponse,
): 'current' | 'superseded' {
  return publication.artifact_generation === publication.current_artifact_generation
    ? 'current'
    : 'superseded';
}

export function artifactPublicationIdentity(
  publication: Pick<
    ArtifactPublicationResponse,
    'run_id' | 'artifact_id' | 'artifact_generation'
  >,
): string {
  return `${publication.run_id}\u0000${publication.artifact_id}\u0000${publication.artifact_generation}`;
}

export function hasActiveArtifactPublicationFilters(
  filters: ArtifactPublicationFilters,
): boolean {
  return Boolean(
    filters.projectId ||
      filters.runId ||
      filters.workflowId ||
      filters.outputType ||
      filters.associatedRunSampleRevisionId ||
      filters.publishedFrom ||
      filters.publishedBefore ||
      !filters.currentOnly ||
      filters.after ||
      filters.limit !== ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE,
  );
}

export function artifactPublicationSearchParams(
  filters: ArtifactPublicationFilters,
): URLSearchParams {
  const params = new URLSearchParams();
  setOptional(params, 'project_id', filters.projectId);
  setOptional(params, 'run_id', filters.runId);
  setOptional(params, 'workflow_id', filters.workflowId);
  setOptional(params, 'output_type', filters.outputType);
  setOptional(
    params,
    'associated_run_sample_revision_id',
    filters.associatedRunSampleRevisionId,
  );
  setOptional(params, 'published_from', filters.publishedFrom);
  setOptional(params, 'published_before', filters.publishedBefore);
  if (!filters.currentOnly) params.set('current_only', 'false');
  setOptional(params, 'after', filters.after);
  if (filters.limit !== ARTIFACT_PUBLICATION_DEFAULT_PAGE_SIZE) {
    params.set('limit', String(filters.limit));
  }
  return params;
}

export function isValidArtifactPublicationRunId(value: string): boolean {
  return RUN_ID.test(value);
}

export function isValidArtifactPublicationArtifactId(value: string): boolean {
  return LOGICAL_ID.test(value);
}

export function isValidArtifactPublicationGeneration(value: string): boolean {
  return ARTIFACT_GENERATION.test(value);
}

function projectArtifactPublication(value: unknown): ArtifactPublicationResponse {
  if (!hasExactKeys(value, PUBLICATION_KEYS)) {
    throw new ArtifactPublicationProtocolError();
  }
  if (
    typeof value.run_id !== 'string' ||
    !RUN_ID.test(value.run_id) ||
    typeof value.project_id !== 'string' ||
    !PROJECT_ID.test(value.project_id) ||
    typeof value.workflow_id !== 'string' ||
    !WORKFLOW_ID.test(value.workflow_id) ||
    typeof value.artifact_id !== 'string' ||
    !LOGICAL_ID.test(value.artifact_id) ||
    typeof value.output_type !== 'string' ||
    !LOGICAL_ID.test(value.output_type) ||
    typeof value.artifact_generation !== 'string' ||
    !ARTIFACT_GENERATION.test(value.artifact_generation) ||
    typeof value.artifact_revision !== 'string' ||
    !ARTIFACT_REVISION.test(value.artifact_revision) ||
    !isTimestamp(value.published_at) ||
    typeof value.current_artifact_generation !== 'string' ||
    !ARTIFACT_GENERATION.test(value.current_artifact_generation) ||
    (value.generation_status !== 'current' &&
      value.generation_status !== 'superseded')
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  const derived =
    value.artifact_generation === value.current_artifact_generation
      ? 'current'
      : 'superseded';
  if (value.generation_status !== derived) {
    throw new ArtifactPublicationProtocolError();
  }
  const runSampleBinding = projectRunSampleBinding(value.run_sample_binding);
  if (
    (value.project_id === LEGACY_PROJECT_ID) !==
    (runSampleBinding.binding_mode === 'legacy_v1')
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  return {
    run_id: value.run_id,
    project_id: value.project_id,
    workflow_id: value.workflow_id,
    artifact_id: value.artifact_id,
    output_type: value.output_type,
    artifact_generation: value.artifact_generation,
    artifact_revision: value.artifact_revision,
    published_at: value.published_at,
    current_artifact_generation: value.current_artifact_generation,
    generation_status: value.generation_status,
    run_sample_binding: runSampleBinding,
  };
}

function projectRunSampleBinding(value: unknown): ArtifactPublicationResponse['run_sample_binding'] {
  if (
    !hasExactKeys(value, [
      'binding_mode',
      'provenance',
      'associated_run_samples',
    ]) ||
    (value.binding_mode !== 'legacy_v1' && value.binding_mode !== 'bound_v1') ||
    (value.provenance !== 'resolved' && value.provenance !== 'unresolved') ||
    !Array.isArray(value.associated_run_samples)
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  const associatedRunSamples = value.associated_run_samples.map((item, index) => {
    const revisionNumber = hasExactKeys(item, [
      'sample_id',
      'sample_revision_id',
      'revision_number',
      'ordinal',
    ]) && typeof item.revision_number === 'number'
      ? item.revision_number
      : null;
    if (
      !hasExactKeys(item, [
        'sample_id',
        'sample_revision_id',
        'revision_number',
        'ordinal',
      ]) ||
      typeof item.sample_id !== 'string' ||
      !SAMPLE_ID.test(item.sample_id) ||
      typeof item.sample_revision_id !== 'string' ||
      !SAMPLE_REVISION_ID.test(item.sample_revision_id) ||
      revisionNumber === null ||
      !Number.isSafeInteger(revisionNumber) ||
      revisionNumber < 1 ||
      item.ordinal !== index
    ) {
      throw new ArtifactPublicationProtocolError();
    }
    return {
      sample_id: item.sample_id,
      sample_revision_id: item.sample_revision_id,
      revision_number: revisionNumber,
      ordinal: item.ordinal,
    };
  });
  if (
    new Set(associatedRunSamples.map((item) => item.sample_revision_id)).size !==
    associatedRunSamples.length
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  const isLegacy = value.binding_mode === 'legacy_v1';
  if (
    (isLegacy &&
      (value.provenance !== 'unresolved' || associatedRunSamples.length !== 0)) ||
    (!isLegacy &&
      (value.provenance !== 'resolved' || associatedRunSamples.length === 0))
  ) {
    throw new ArtifactPublicationProtocolError();
  }
  return {
    binding_mode: value.binding_mode,
    provenance: value.provenance,
    associated_run_samples: associatedRunSamples,
  };
}

function optionalToken(
  searchParams: URLSearchParams,
  key: string,
  pattern: RegExp,
  maxLength?: number,
): string | null | false {
  const value = searchParams.get(key);
  if (value === null) return null;
  if (
    value.trim() !== value ||
    (maxLength !== undefined && value.length > maxLength) ||
    !pattern.test(value)
  ) {
    return false;
  }
  return value;
}

function optionalTimestamp(
  searchParams: URLSearchParams,
  key: string,
): string | null | false {
  const value = searchParams.get(key);
  if (value === null) return null;
  return canonicalTimestamp(value) ?? false;
}

function isTimestamp(value: unknown): value is string {
  return parseTimestamp(value) !== null;
}

interface TimestampParts {
  year: number;
  month: number;
  day: number;
  hour: number;
  minute: number;
  second: number;
  fraction: string;
  offsetMinutes: number;
}

function parseTimestamp(value: unknown): TimestampParts | null {
  if (typeof value !== 'string' || value.length > 64) return null;
  const match = TIMEZONE_TIMESTAMP.exec(value);
  if (!match) return null;
  const year = Number(match[1]);
  const month = Number(match[2]);
  const day = Number(match[3]);
  const hour = Number(match[4]);
  const minute = Number(match[5]);
  const second = Number(match[6] ?? '0');
  const fraction = match[7] ?? '';
  const offsetHour = Number(match[10] ?? '0');
  const offsetMinute = Number(match[11] ?? '0');
  const leapYear = year % 4 === 0 && (year % 100 !== 0 || year % 400 === 0);
  const daysInMonth = [
    31,
    leapYear ? 29 : 28,
    31,
    30,
    31,
    30,
    31,
    31,
    30,
    31,
    30,
    31,
  ];
  if (
    year < 1 ||
    month < 1 ||
    month > 12 ||
    day < 1 ||
    day > daysInMonth[month - 1]! ||
    hour > 23 ||
    minute > 59 ||
    second > 59 ||
    offsetHour > 23 ||
    offsetMinute > 59
  ) {
    return null;
  }
  const offsetDirection = match[9] === '-' ? -1 : match[9] === '+' ? 1 : 0;
  return {
    year,
    month,
    day,
    hour,
    minute,
    second,
    fraction,
    offsetMinutes: offsetDirection * (offsetHour * 60 + offsetMinute),
  };
}

function canonicalTimestamp(value: unknown): string | null {
  const parts = parseTimestamp(value);
  if (!parts) return null;
  const wallClock = new Date(0);
  wallClock.setUTCFullYear(parts.year, parts.month - 1, parts.day);
  wallClock.setUTCHours(parts.hour, parts.minute, parts.second, 0);
  const instant = new Date(
    wallClock.getTime() - parts.offsetMinutes * 60 * 1000,
  );
  const utcYear = instant.getUTCFullYear();
  if (utcYear < 1 || utcYear > 9999) return null;
  const date = [
    String(utcYear).padStart(4, '0'),
    String(instant.getUTCMonth() + 1).padStart(2, '0'),
    String(instant.getUTCDate()).padStart(2, '0'),
  ].join('-');
  const time = [
    String(instant.getUTCHours()).padStart(2, '0'),
    String(instant.getUTCMinutes()).padStart(2, '0'),
    String(instant.getUTCSeconds()).padStart(2, '0'),
  ].join(':');
  return `${date}T${time}.${parts.fraction.padEnd(6, '0')}Z`;
}

function isSafeArtifactPublicationCursor(value: unknown): value is string | null {
  return (
    value === null ||
    (typeof value === 'string' &&
      value.length <= 1024 &&
      ARTIFACT_PUBLICATION_CURSOR.test(value))
  );
}

function hasExactKeys(
  value: unknown,
  expected: Iterable<string>,
): value is Record<string, unknown> {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) {
    return false;
  }
  const expectedKeys = new Set(expected);
  const keys = Object.keys(value);
  return keys.length === expectedKeys.size && keys.every((key) => expectedKeys.has(key));
}

function setOptional(
  params: URLSearchParams,
  key: string,
  value: string | null,
): void {
  if (value !== null) params.set(key, value);
}
