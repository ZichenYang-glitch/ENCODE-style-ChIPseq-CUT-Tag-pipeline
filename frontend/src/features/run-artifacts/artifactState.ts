import type {
  ArtifactReferenceResponse,
  RunArtifactsResponse,
} from '../../api/generated/models';
import type { RunEventResponse } from '../../api/runTypes';

const ARTIFACT_ID_PATTERN = /^[A-Za-z][A-Za-z0-9_.-]{0,127}$/;
const REASON_CODE_PATTERN = /^[A-Z][A-Z0-9_]{0,127}$/;
const byteFormatter = new Intl.NumberFormat(undefined, {
  maximumFractionDigits: 1,
});

export type ArtifactExtractionOutcome =
  | { kind: 'pending' }
  | { kind: 'unconfirmed' }
  | { kind: 'indexed'; count: number }
  | { kind: 'failed'; reasonCode: string | null };

export function isValidArtifactId(value: string | null): value is string {
  return value !== null && ARTIFACT_ID_PATTERN.test(value);
}

export function artifactExtractionOutcome(
  events: RunEventResponse[],
  truncated = false,
): ArtifactExtractionOutcome {
  const latest = events
    .filter((item) =>
      ['artifacts_indexed', 'artifact_extraction_failed'].includes(item.event_type),
    )
    .sort((left, right) => right.sequence - left.sequence)[0];

  if (!latest) return truncated ? { kind: 'unconfirmed' } : { kind: 'pending' };
  if (latest.event_type === 'artifact_extraction_failed') {
    const reasonCode = latest.context.reason_code;
    return {
      kind: 'failed',
      reasonCode:
        typeof reasonCode === 'string' && REASON_CODE_PATTERN.test(reasonCode)
          ? reasonCode
          : null,
    };
  }

  const count = latest.context.artifact_count;
  if (typeof count !== 'number' || !Number.isSafeInteger(count) || count < 0) {
    return { kind: 'unconfirmed' };
  }
  return { kind: 'indexed', count };
}

export function flattenArtifactPages(
  pages: RunArtifactsResponse[],
): ArtifactReferenceResponse[] {
  const seen = new Set<string>();
  const result: ArtifactReferenceResponse[] = [];
  for (const page of pages) {
    for (const artifact of page.artifacts ?? []) {
      if (seen.has(artifact.artifact_id)) continue;
      seen.add(artifact.artifact_id);
      result.push(artifact);
    }
  }
  return result;
}

export function safeNextArtifactCursor(
  lastPage: RunArtifactsResponse,
  allPages: RunArtifactsResponse[],
  pageParams: unknown[],
): string | undefined {
  const next = lastPage.next_cursor;
  if (typeof next !== 'string' || !isValidArtifactId(next)) return undefined;
  const currentPageParam = pageParams.at(-1);
  if (currentPageParam === next) return undefined;
  if (pageParams.some((value) => value === next)) return undefined;
  if (
    allPages.slice(0, -1).some((page) => page.next_cursor === next)
  ) {
    return undefined;
  }
  return next;
}

export function formatBytes(sizeBytes: number): string {
  if (!Number.isFinite(sizeBytes) || sizeBytes < 0) return 'Unknown size';
  if (sizeBytes < 1024) return `${byteFormatter.format(sizeBytes)} B`;
  const units = ['KiB', 'MiB', 'GiB', 'TiB'];
  let value = sizeBytes / 1024;
  let unit = units[0];
  for (let index = 1; index < units.length && value >= 1024; index += 1) {
    value /= 1024;
    unit = units[index];
  }
  return `${byteFormatter.format(value)} ${unit}`;
}

export function formatProducedTime(value: string, locale?: string): string {
  const timestamp = new Date(value);
  if (Number.isNaN(timestamp.getTime())) return 'Unknown time';
  return new Intl.DateTimeFormat(locale, {
    dateStyle: 'medium',
    timeStyle: 'short',
  }).format(timestamp);
}
