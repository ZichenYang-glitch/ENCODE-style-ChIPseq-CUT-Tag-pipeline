import type { HTMLAttributes } from 'react';
import type { WashuTrackConfig } from 'wuepgg';
import { ApiError } from '../../api/fetcher';
import type { ArchivedTrack } from '../../api/archivedTracks';
import { archivedTrackResourcePath } from '../../api/archivedTrackPaths';
import { assemblyKey, PLAIN_PEAK_MAX_BYTES } from '../../viewer/local/contract';
import { scanPlainPeak } from '../../viewer/local/sliceRead';
import { viewerScaffold } from '../run-artifacts/browserWindow';

const PROJECT_ID = /^prj_[0-9a-f]{32}$/;
const EXPERIMENT_ID = /^exp_[0-9a-f]{32}$/;
const ARCHIVE_VERSION_ID = /^[A-Za-z][A-Za-z0-9_.-]{0,127}$/;
const REGION = /^[A-Za-z0-9_.:-]{1,128}:[0-9]{1,10}-[0-9]{1,10}$/;

export const ARCHIVE_PAGE_LIMIT = 50;

export const ARCHIVE_ERROR_CODES = [
  'archive_version_not_found',
  'archive_version_conflict',
  'archive_offline',
  'archive_integrity',
  'archive_invalid',
  'archive_unauthorized',
] as const;

export type ArchiveErrorCode = (typeof ARCHIVE_ERROR_CODES)[number];

export type ArchiveViewerRequest =
  | { kind: 'absent' }
  | { kind: 'invalid' }
  | {
      kind: 'ready';
      projectId: string;
      experimentId: string;
      archiveVersionId: string;
    };

const ERROR_TEXT: Record<ArchiveErrorCode, string> = {
  archive_version_not_found: '这份封存版本不存在。页面没有改读另一版。',
  archive_version_conflict: '封存版本发生冲突。页面没有改读另一版。',
  archive_offline: '封存副本离线。页面没有改读另一版。',
  archive_integrity: '封存副本未通过校验。页面没有改读另一版。',
  archive_invalid: '封存版本无法读取。页面没有改读另一版。',
  archive_unauthorized: '无权读取这份封存版本。页面没有改读另一版。',
};

export function archiveErrorCopy(code: string | null | undefined): {
  code: ArchiveErrorCode;
  message: string;
} {
  const known = ARCHIVE_ERROR_CODES.find((item) => item === code) ?? 'archive_invalid';
  return { code: known, message: ERROR_TEXT[known] };
}

export function archiveErrorFromUnknown(error: unknown): { code: ArchiveErrorCode; message: string } {
  if (error instanceof ApiError) {
    if (error.status === 401 || error.status === 403) return archiveErrorCopy('archive_unauthorized');
    return archiveErrorCopy(error.code);
  }
  return archiveErrorCopy('archive_invalid');
}

export function archiveViewerRequest(search: string): ArchiveViewerRequest {
  const params = new URLSearchParams(search.startsWith('?') ? search.slice(1) : search);
  const names = ['project', 'experiment', 'archive'] as const;
  if (names.every((name) => !params.has(name))) return { kind: 'absent' };
  if (names.some((name) => params.getAll(name).length !== 1)) return { kind: 'invalid' };
  const projectId = params.get('project') ?? '';
  const experimentId = params.get('experiment') ?? '';
  const archiveVersionId = params.get('archive') ?? '';
  if (
    !PROJECT_ID.test(projectId) ||
    !EXPERIMENT_ID.test(experimentId) ||
    !ARCHIVE_VERSION_ID.test(archiveVersionId)
  ) {
    return { kind: 'invalid' };
  }
  return { kind: 'ready', projectId, experimentId, archiveVersionId };
}

export function archivedViewerHref(
  track: Pick<ArchivedTrack, 'project_id' | 'experiment_id' | 'archive_version_id'>,
): string {
  const params = new URLSearchParams({
    project: track.project_id,
    experiment: track.experiment_id,
    archive: track.archive_version_id,
  });
  return `/viewer.html?${params.toString()}`;
}

export function archivedTrackOpenName(
  track: Pick<ArchivedTrack, 'display_name' | 'version_number'>,
): string {
  return `在 WashU 查看 ${track.display_name} 第 ${track.version_number} 版`;
}

export function archivedVersionText(versionNumber: number): string {
  return `第 ${versionNumber} 版`;
}

export function archivedSourceText(track: ArchivedTrack): string {
  if (track.source.kind === 'run') {
    const base = `运行 ${track.source.run_id} · 产物 ${track.source.artifact_id} · 代次 ${track.source.generation} · 修订 ${track.source.revision}`;
    if (
      track.source.index_artifact_id !== undefined &&
      track.source.index_generation !== undefined &&
      track.source.index_revision !== undefined
    ) {
      return `${base} · 索引产物 ${track.source.index_artifact_id} · 索引代次 ${track.source.index_generation} · 索引修订 ${track.source.index_revision}`;
    }
    return base;
  }
  const base = `登记 ${track.source.result_id} · 修订 ${track.source.result_revision_id}`;
  if (
    track.source.index_result_id !== undefined &&
    track.source.index_result_revision_id !== undefined
  ) {
    return `${base} · 索引登记 ${track.source.index_result_id} · 索引修订 ${track.source.index_result_revision_id}`;
  }
  return base;
}

export function archivedContentPath(track: ArchivedTrack): string {
  return `${archivedTrackResourcePath(track.project_id, track.experiment_id, track.archive_version_id)}/content`;
}

export function archivedIndexPath(track: ArchivedTrack): string | null {
  if (!track.has_index) return null;
  return `${archivedContentPath(track)}.tbi`;
}

export function sameOriginApiUrl(path: string): string | null {
  if (!path.startsWith('/api/')) return null;
  const configured = import.meta.env.VITE_API_BASE_URL;
  if (typeof configured === 'string' && configured.trim() !== '') {
    try {
      const url = new URL(configured, window.location.origin);
      if (url.origin !== window.location.origin) return null;
    } catch {
      return null;
    }
  }
  return new URL(path, window.location.origin).toString();
}

export function archivedHubTracks(track: ArchivedTrack, contentUrl: string, textFile?: File | null): WashuTrackConfig[] {
  const data: WashuTrackConfig = {
    type: track.track_type,
    name: `${track.display_name} ${archivedVersionText(track.version_number)}`,
    url: contentUrl,
    waitToUpdate: true,
  };
  if (textFile) Object.assign(data, { fileObj: textFile, url: '', isText: true, textConfig: {} });
  if (track.track_type === 'longrange') data.options = { clampHeight: true };
  return [
    { type: 'ruler', name: 'Ruler', waitToUpdate: true },
    data,
  ];
}

export function archivedTrackScaffold(track: ArchivedTrack): string | null {
  return viewerScaffold(track.assembly);
}

export function isDrawableRegion(region: string): boolean {
  return REGION.test(region);
}

export function archivedTrackDataProps(track: ArchivedTrack): HTMLAttributes<HTMLElement> {
  const props: Record<string, string> = {
    'data-archive-version-id': track.archive_version_id,
    'data-version-number': String(track.version_number),
    'data-display-name': track.display_name,
    'data-track-key': track.track_key,
    'data-experiment-revision-id': track.experiment_revision_id,
    'data-assembly': track.assembly,
    'data-track-type': track.track_type,
    'data-has-index': track.has_index ? 'true' : 'false',
    'data-project-id': track.project_id,
    'data-experiment-id': track.experiment_id,
    'data-source-kind': track.source.kind,
  };
  if (track.source.kind === 'run') {
    props['data-run-id'] = track.source.run_id;
    props['data-artifact-id'] = track.source.artifact_id;
    props['data-generation'] = track.source.generation;
    props['data-revision'] = track.source.revision;
    props['data-index-artifact-id'] = track.source.index_artifact_id ?? '';
    props['data-index-generation'] = track.source.index_generation ?? '';
    props['data-index-revision'] = track.source.index_revision ?? '';
  } else {
    props['data-result-id'] = track.source.result_id;
    props['data-result-revision-id'] = track.source.result_revision_id;
    props['data-index-result-id'] = track.source.index_result_id ?? '';
    props['data-index-result-revision-id'] = track.source.index_result_revision_id ?? '';
  }
  return props as HTMLAttributes<HTMLElement>;
}

export async function preflightArchivedBytes(
  track: ArchivedTrack,
  contentUrl: string,
  indexUrl: string | null,
  signal: AbortSignal,
): Promise<File | null> {
  if (track.track_type === 'bed' && !track.has_index) {
    // The pinned viewer reads unindexed text through a File. Bound the network
    // request as well as the decoded peak file before handing it to that reader.
    const response = await archiveRange(contentUrl, `bytes=0-${PLAIN_PEAK_MAX_BYTES - 1}`, signal);
    const range = /^bytes 0-(\d+)\/(\d+)$/.exec(response.headers.get('content-range') ?? '');
    const size = range === null ? 0 : Number(range[2]);
    if (size <= 0 || size > PLAIN_PEAK_MAX_BYTES || Number(range?.[1]) !== size - 1) {
      await response.body?.cancel();
      throw new ApiError(409, 'archive_invalid', 'Archived text exceeds the supported range.');
    }
    const bytes = await response.arrayBuffer();
    if (bytes.byteLength !== size) throw new ApiError(409, 'archive_integrity', 'Incomplete archived file.');
    const scanned = await scanPlainPeak(new File([bytes], track.display_name), 'BED3');
    if ('error' in scanned) throw new ApiError(409, 'archive_invalid', scanned.error);
    if (scanned.ok.declaredAssembly && assemblyKey(scanned.ok.declaredAssembly) !== assemblyKey(track.assembly)) {
      throw new ApiError(409, 'archive_version_conflict', 'The track header and archive assembly disagree.');
    }
    return scanned.ok.dataFile;
  }
  const first = await archiveRange(contentUrl, 'bytes=0-0', signal);
  await first.body?.cancel();
  if (indexUrl !== null) {
    const index = await archiveRange(indexUrl, 'bytes=0-1', signal);
    await index.body?.cancel();
  }
  return null;
}

async function archiveRange(url: string, range: string, signal: AbortSignal): Promise<Response> {
  let response: Response;
  try {
    response = await fetch(url, {
      method: 'GET',
      headers: { Range: range },
      credentials: 'same-origin',
      cache: 'no-store',
      signal,
    });
  } catch {
    throw new ApiError(0, 'archive_offline', 'Archived copy could not be reached.');
  }
  if (response.status === 206) return response;
  if (response.ok) {
    await response.body?.cancel();
    throw new ApiError(409, 'archive_invalid', 'The archived response did not honor its byte range.');
  }
  const failure = await failureFromResponse(response);
  throw new ApiError(response.status, failure.code, failure.message);
}

async function failureFromResponse(response: Response): Promise<{ code: ArchiveErrorCode; message: string }> {
  if (response.status === 401 || response.status === 403) return archiveErrorCopy('archive_unauthorized');
  let code: string | null = null;
  try {
    const body = (await response.json()) as { issues?: Array<{ code?: unknown }> };
    const issueCode = body.issues?.[0]?.code;
    if (typeof issueCode === 'string') code = issueCode;
  } catch {
    code = null;
  }
  if (response.status === 404) return archiveErrorCopy(code ?? 'archive_version_not_found');
  if (response.status === 409) return archiveErrorCopy(code ?? 'archive_version_conflict');
  return archiveErrorCopy(code);
}
