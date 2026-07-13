import { describe, expect, it } from 'vitest';
import { buildDraftReview } from './draft';

describe('buildDraftReview', () => {
  it('serializes config/options deterministically and sample keys in schema order', () => {
    const result = buildDraftReview(
      {
        z: 1,
        'é': 2,
        a: 'first',
        _: 4,
        A: 5,
        nested: { b: true, a: false },
      },
      [
        {
          id: 'client-only',
          values: { fastq_1: '/data/S1.fastq.gz', sample: 'S1', layout: 'SE' },
        },
      ],
      { strict_inputs: false },
      ['sample', 'fastq_1', 'layout'],
      2_097_152,
      4_096,
    );

    expect(result.ok).toBe(true);
    if (!result.ok) return;
    expect(result.payload.samples).toEqual([
      { sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' },
    ]);
    expect(Object.keys(JSON.parse(result.serialized).config)).toEqual([
      'A',
      '_',
      'a',
      'nested',
      'z',
      'é',
    ]);
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
        4_096,
      ),
    ).toMatchObject({ ok: false, code: 'DRAFT_NOT_JSON_SAFE' });
    expect(
      buildDraftReview({ emoji: '😀' }, [], {}, [], 10, 4_096),
    ).toMatchObject({
      ok: false,
      code: 'DRAFT_TOO_LARGE',
    });
  });

  it('rejects unsafe manually edited sample cells before serialization', () => {
    expect(
      buildDraftReview(
        {},
        [{ id: 'row-1', values: { sample: 'bad\nvalue' } }],
        {},
        ['sample'],
        2_097_152,
        4_096,
      ),
    ).toMatchObject({ ok: false, code: 'DRAFT_SAMPLE_INVALID' });
    expect(
      buildDraftReview(
        {},
        [{ id: 'row-1', values: { sample: 'x'.repeat(4_097) } }],
        {},
        ['sample'],
        2_097_152,
        4_096,
      ),
    ).toMatchObject({ ok: false, code: 'DRAFT_SAMPLE_INVALID' });
  });
});
