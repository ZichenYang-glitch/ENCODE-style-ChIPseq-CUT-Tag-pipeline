import { describe, expect, it } from 'vitest';
import { buildDraftReview } from './draft';

describe('buildDraftReview', () => {
  it('serializes config/options deterministically and sample keys in schema order', () => {
    const result = buildDraftReview(
      { z: 1, nested: { b: true, a: false }, a: 'first' },
      [
        {
          id: 'client-only',
          values: { fastq_1: '/data/S1.fastq.gz', sample: 'S1', layout: 'SE' },
        },
      ],
      { strict_inputs: false },
      ['sample', 'fastq_1', 'layout'],
      2_097_152,
    );

    expect(result.ok).toBe(true);
    if (!result.ok) return;
    expect(result.payload.samples).toEqual([
      { sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' },
    ]);
    expect(result.serialized.indexOf('"a"')).toBeLessThan(
      result.serialized.indexOf('"z"'),
    );
    expect(result.serialized).not.toContain('client-only');
  });

  it('rejects non-JSON values and a UTF-8 payload over the advertised ceiling', () => {
    expect(
      buildDraftReview(
        { invalid: Number.POSITIVE_INFINITY },
        [],
        {},
        [],
        100,
      ),
    ).toMatchObject({ ok: false, code: 'DRAFT_NOT_JSON_SAFE' });
    expect(buildDraftReview({ emoji: '😀' }, [], {}, [], 10)).toMatchObject({
      ok: false,
      code: 'DRAFT_TOO_LARGE',
    });
  });
});
