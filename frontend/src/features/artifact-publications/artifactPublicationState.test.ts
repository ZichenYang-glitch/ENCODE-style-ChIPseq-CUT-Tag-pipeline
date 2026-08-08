import { describe, expect, it } from 'vitest';
import type {
  ArtifactPublicationListResponse,
  ArtifactPublicationResponse,
} from '../../api/generated/models';
import {
  artifactPublicationPaginationAnomaly,
  artifactPublicationQueryKey,
  artifactPublicationSearchParams,
  flattenArtifactPublicationPages,
  parseArtifactPublicationDetailIdentity,
  parseArtifactPublicationFilters,
  validateArtifactPublicationDetailResponse,
  validateArtifactPublicationListResponse,
  ArtifactPublicationProtocolError,
} from './artifactPublicationState';

const PROJECT_ID = `prj_${'1'.repeat(32)}`;
const SAMPLE_ID = `smp_${'2'.repeat(32)}`;
const SAMPLE_REVISION_ID = `smpr_${'3'.repeat(32)}`;
const GENERATION_A = `artifactgen-${'a'.repeat(64)}`;
const GENERATION_B = `artifactgen-${'b'.repeat(64)}`;
const REVISION = `artifactrev-${'c'.repeat(64)}`;

function publication(
  overrides: Partial<ArtifactPublicationResponse> = {},
): ArtifactPublicationResponse {
  return {
    run_id: 'run-a',
    project_id: PROJECT_ID,
    workflow_id: 'workflow-a',
    artifact_id: 'artifact-a',
    output_type: 'counts',
    artifact_generation: GENERATION_A,
    artifact_revision: REVISION,
    published_at: '2026-08-07T08:00:00Z',
    current_artifact_generation: GENERATION_A,
    generation_status: 'current',
    run_sample_binding: {
      binding_mode: 'bound_v1',
      provenance: 'resolved',
      associated_run_samples: [
        {
          sample_id: SAMPLE_ID,
          sample_revision_id: SAMPLE_REVISION_ID,
          revision_number: 4,
          ordinal: 0,
        },
      ],
    },
    ...overrides,
  };
}

function page(
  items: ArtifactPublicationResponse[] = [publication()],
  nextCursor: string | null = null,
): ArtifactPublicationListResponse {
  return { ok: true, publications: items, next_cursor: nextCursor, issues: [] };
}

describe('artifact publication URL and protocol state', () => {
  it('normalizes every list filter and includes each in the query key', () => {
    const parsed = parseArtifactPublicationFilters(
      new URLSearchParams({
        project_id: PROJECT_ID,
        run_id: 'run-a',
        workflow_id: 'workflow-a',
        output_type: 'counts',
        associated_run_sample_revision_id: SAMPLE_REVISION_ID,
        published_from: '2026-08-07T16:00:00.123456+08:00',
        published_before: '2026-08-07T16:00:00.123457+08:00',
        current_only: 'false',
        after: 'artifactpubcur_abc_123',
        limit: '100',
      }),
    );
    expect(parsed.ok).toBe(true);
    if (!parsed.ok) return;
    expect(parsed.filters).toEqual({
      projectId: PROJECT_ID,
      runId: 'run-a',
      workflowId: 'workflow-a',
      outputType: 'counts',
      associatedRunSampleRevisionId: SAMPLE_REVISION_ID,
      publishedFrom: '2026-08-07T08:00:00.123456Z',
      publishedBefore: '2026-08-07T08:00:00.123457Z',
      currentOnly: false,
      after: 'artifactpubcur_abc_123',
      limit: 100,
    });
    expect(artifactPublicationQueryKey(parsed.filters).slice(1)).toEqual([
      PROJECT_ID,
      'run-a',
      'workflow-a',
      'counts',
      SAMPLE_REVISION_ID,
      '2026-08-07T08:00:00.123456Z',
      '2026-08-07T08:00:00.123457Z',
      false,
      'artifactpubcur_abc_123',
      100,
    ]);
    const canonicalUrl = artifactPublicationSearchParams(parsed.filters);
    expect(canonicalUrl.get('published_from')).toBe(
      '2026-08-07T08:00:00.123456Z',
    );
    expect(canonicalUrl.get('published_before')).toBe(
      '2026-08-07T08:00:00.123457Z',
    );
    expect(canonicalUrl.get('after')).toBe('artifactpubcur_abc_123');
    expect(parseArtifactPublicationFilters(canonicalUrl)).toEqual(parsed);
  });

  it('pads missing timestamp precision without losing microseconds', () => {
    const parsed = parseArtifactPublicationFilters(
      new URLSearchParams({
        published_from: '2026-08-07T08:00Z',
        published_before: '2026-08-07T08:00:00.1Z',
      }),
    );
    expect(parsed).toEqual({
      ok: true,
      filters: expect.objectContaining({
        publishedFrom: '2026-08-07T08:00:00.000000Z',
        publishedBefore: '2026-08-07T08:00:00.100000Z',
      }),
    });
  });

  it.each([
    'run_id=run-a&run_id=run-b',
    'run_id=%2Fprivate%2Frun',
    'workflow_id=workflow-a%0Asecret',
    'project_id=prj_invalid',
    'associated_run_sample_revision_id=smpr_invalid',
    'published_from=2026-08-08T08%3A00%3A00Z&published_before=2026-08-07T08%3A00%3A00Z',
    'published_from=2026-08-07',
    'published_from=2026-02-31T08%3A00%3A00Z',
    'published_from=2026-08-07T08%3A00%3A00.1234567Z',
    'current_only=maybe',
    'after=cursor-wrong',
    'limit=101',
    'limit=050',
    'unknown=value',
  ])('rejects unsafe, ambiguous, or noncanonical URL state: %s', (query) => {
    expect(parseArtifactPublicationFilters(new URLSearchParams(query))).toEqual({
      ok: false,
      filters: null,
    });
  });

  it('accepts exactly one safe detail identity and rejects extra or repeated state', () => {
    expect(
      parseArtifactPublicationDetailIdentity(
        'run-a',
        'artifact-a',
        new URLSearchParams({ generation: GENERATION_A }),
      ),
    ).toEqual({
      ok: true,
      runId: 'run-a',
      artifactId: 'artifact-a',
      generation: GENERATION_A,
    });
    expect(
      parseArtifactPublicationDetailIdentity(
        'run-a',
        'artifact-a',
        new URLSearchParams(`generation=${GENERATION_A}&generation=${GENERATION_B}`),
      ),
    ).toEqual({ ok: false });
    expect(
      parseArtifactPublicationDetailIdentity(
        '/private',
        'artifact-a',
        new URLSearchParams({ generation: GENERATION_A }),
      ),
    ).toEqual({ ok: false });
  });

  it('strictly validates status evidence, allowlisted fields, and detail identity', () => {
    expect(validateArtifactPublicationListResponse(page(), undefined).publications[0])
      .toMatchObject({ generation_status: 'current' });
    expect(() =>
      validateArtifactPublicationListResponse(
        page([publication({ current_artifact_generation: GENERATION_B })]),
        undefined,
      ),
    ).toThrow(ArtifactPublicationProtocolError);
    expect(() =>
      validateArtifactPublicationListResponse(
        page([{ ...publication(), uri: 'file:///private/result' } as ArtifactPublicationResponse]),
        undefined,
      ),
    ).toThrow(ArtifactPublicationProtocolError);
    expect(() =>
      validateArtifactPublicationDetailResponse(
        { ok: true, publication: publication(), issues: [] },
        'run-other',
        'artifact-a',
        GENERATION_A,
      ),
    ).toThrow(ArtifactPublicationProtocolError);
  });

  it('detects repeated cursors and cross-page compound identities', () => {
    const first = page([publication()], 'artifactpubcur_next');
    const duplicate = page([publication()], null);
    expect(
      artifactPublicationPaginationAnomaly(
        [first, duplicate],
        [undefined, 'artifactpubcur_next'],
      ),
    ).toBe(true);
    expect(flattenArtifactPublicationPages([first, duplicate])).toHaveLength(1);
    expect(() =>
      validateArtifactPublicationListResponse(
        page([publication({ artifact_id: 'artifact-b' })], 'artifactpubcur_next'),
        'artifactpubcur_next',
      ),
    ).toThrow(ArtifactPublicationProtocolError);
  });
});
