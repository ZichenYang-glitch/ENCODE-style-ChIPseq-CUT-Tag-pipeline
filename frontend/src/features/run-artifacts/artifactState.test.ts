import { describe, expect, it } from 'vitest';
import type {
  ArtifactReferenceResponse,
  RunArtifactsResponse,
} from '../../api/generated/models';
import type { RunEventResponse } from '../../api/runTypes';
import {
  artifactExtractionOutcome,
  flattenArtifactPages,
  formatBytes,
  formatProducedTime,
  isValidArtifactId,
  safeNextArtifactCursor,
} from './artifactState';

function event(
  sequence: number,
  eventType: string,
  context: Record<string, unknown> = {},
): RunEventResponse {
  return {
    event_id: `event-${sequence}`,
    run_id: 'run-1',
    sequence,
    event_type: eventType,
    timestamp: '2026-07-12T12:00:00Z',
    status: 'succeeded',
    stage: 'artifact_extraction',
    message: eventType,
    context,
    issue: null,
  };
}

function artifact(id: string): ArtifactReferenceResponse {
  return {
    artifact_id: id,
    run_id: 'run-1',
    artifact_type: 'file',
    name: `${id}.tsv`,
    uri: `run://runs/run-1/artifacts/${id}`,
    mime_type: 'text/tab-separated-values',
    produced_at: '2026-07-12T12:00:00Z',
    revision: `artifactrev-${'a'.repeat(64)}`,
    relative_path: `results/${id}.tsv`,
    output_type: id,
    size_bytes: 1024,
    metadata: {},
  };
}

function page(
  values: ArtifactReferenceResponse[],
  nextCursor: string | null = null,
): RunArtifactsResponse {
  return {
    ok: true,
    run_id: 'run-1',
    artifacts: values,
    next_cursor: nextCursor,
    issues: [],
  };
}

describe('artifactExtractionOutcome', () => {
  it('uses the latest extraction outcome and validates indexed counts', () => {
    expect(artifactExtractionOutcome([])).toEqual({ kind: 'pending' });
    expect(artifactExtractionOutcome([], true)).toEqual({ kind: 'unconfirmed' });
    expect(
      artifactExtractionOutcome([
        event(3, 'artifact_extraction_failed', { reason_code: 'OLD_FAILURE' }),
        event(4, 'artifacts_indexed', { artifact_count: 2 }),
      ]),
    ).toEqual({ kind: 'indexed', count: 2 });
    expect(
      artifactExtractionOutcome([
        event(5, 'artifacts_indexed', { artifact_count: 2 }),
        event(6, 'artifact_extraction_failed', { reason_code: 'SAFE_CODE' }),
      ]),
    ).toEqual({ kind: 'failed', reasonCode: 'SAFE_CODE' });
  });

  it.each([-1, 1.5, Number.MAX_SAFE_INTEGER + 1, '0', null])(
    'does not treat invalid artifact_count %p as indexed',
    (artifactCount) => {
      expect(
        artifactExtractionOutcome([
          event(1, 'artifacts_indexed', { artifact_count: artifactCount }),
        ]),
      ).toEqual({ kind: 'unconfirmed' });
    },
  );

  it.each(['/private/path', 'SECRET=value', 'ValueError: private detail', 'lowercase'])(
    'does not expose malformed reason code %s',
    (reasonCode) => {
      expect(
        artifactExtractionOutcome([
          event(1, 'artifact_extraction_failed', { reason_code: reasonCode }),
        ]),
      ).toEqual({ kind: 'failed', reasonCode: null });
    },
  );
});

describe('artifact pagination helpers', () => {
  it('deduplicates artifacts across pages without changing first-seen order', () => {
    expect(
      flattenArtifactPages([
        page([artifact('artifact-a'), artifact('artifact-b')], 'artifact-b'),
        page([artifact('artifact-b'), artifact('artifact-c')]),
      ]).map((item) => item.artifact_id),
    ).toEqual(['artifact-a', 'artifact-b', 'artifact-c']);
  });

  it('stops empty, malformed, current, or already-seen cursors', () => {
    const first = page([artifact('artifact-a')], 'artifact-a');
    expect(safeNextArtifactCursor(first, [first], [undefined])).toBe('artifact-a');
    expect(
      safeNextArtifactCursor(first, [first], ['artifact-a']),
    ).toBeUndefined();
    expect(
      safeNextArtifactCursor(page([artifact('artifact-a')], ''), [first], [undefined]),
    ).toBeUndefined();
    expect(
      safeNextArtifactCursor(
        page([artifact('artifact-b')], '../escape'),
        [first],
        [undefined],
      ),
    ).toBeUndefined();
    const repeated = page([artifact('artifact-c')], 'artifact-a');
    expect(
      safeNextArtifactCursor(
        repeated,
        [first, repeated],
        [undefined, 'artifact-a'],
      ),
    ).toBeUndefined();
    expect(
      safeNextArtifactCursor(
        page([artifact('artifact-b')], 'artifact-b'),
        [first],
        [undefined],
      ),
    ).toBe('artifact-b');
  });
});

describe('artifact public formatting', () => {
  it.each([
    ['artifact-a', true],
    ['A.b-c_1', true],
    ['../escape', false],
    ['artifact/a', false],
    ['', false],
    [`a${'x'.repeat(128)}`, false],
  ])('validates artifact ID %s', (value, expected) => {
    expect(isValidArtifactId(value)).toBe(expected);
  });

  it('formats bytes with standard IEC units', () => {
    expect(formatBytes(0)).toBe('0 B');
    expect(formatBytes(1024)).toBe('1 KiB');
    expect(formatBytes(1536)).toBe('1.5 KiB');
  });

  it('formats valid timestamps and safely handles invalid values', () => {
    expect(formatProducedTime('2026-07-12T12:00:00Z', 'en-US')).toContain('2026');
    expect(formatProducedTime('not-a-date', 'en-US')).toBe('Unknown time');
  });
});
