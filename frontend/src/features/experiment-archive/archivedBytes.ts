import { ApiError } from '../../api/fetcher';
import type { ArchivedTrack } from '../../api/archivedTracks';
import { assemblyKey, PLAIN_PEAK_MAX_BYTES } from '../../viewer/local/contract';
import { scanPlainPeak } from '../../viewer/local/sliceRead';
import { archiveErrorCopy, type ArchiveErrorCode } from './archiveView';

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
