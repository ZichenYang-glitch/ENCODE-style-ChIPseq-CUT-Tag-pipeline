import { describe, expect, it } from 'vitest';
import { validateSampleRows } from './sampleValidation';

describe('validateSampleRows', () => {
  it.each([
    ['NUL', 'private-sentinel\0value'],
    ['tab', 'private-sentinel\tvalue'],
    ['carriage return', 'private-sentinel\rvalue'],
    ['line feed', 'private-sentinel\nvalue'],
  ])('rejects a manual cell containing %s without echoing its value', (_label, value) => {
    const result = validateSampleRows(
      [{ sample: value }],
      ['sample'],
      4_096,
    );

    expect(result).toMatchObject({
      code: 'SAMPLE_CELL_INVALID',
      row: 1,
      column: 'sample',
    });
    expect(Object.keys(result ?? {}).sort()).toEqual([
      'code',
      'column',
      'message',
      'row',
    ]);
    expect(result).not.toHaveProperty('value');
    expect(result?.message).not.toContain('private-sentinel');
  });

  it('rejects non-strings and counts Unicode code points for the cell ceiling', () => {
    expect(
      validateSampleRows(
        [{ sample: 42 }],
        ['sample'],
        4_096,
      ),
    ).toMatchObject({ code: 'SAMPLE_CELL_INVALID', row: 1, column: 'sample' });
    expect(
      validateSampleRows(
        [{ sample: '😀'.repeat(4_096) }],
        ['sample'],
        4_096,
      ),
    ).toBeNull();
    expect(
      validateSampleRows(
        [{ sample: '😀'.repeat(4_097) }],
        ['sample'],
        4_096,
      ),
    ).toMatchObject({
      code: 'SAMPLE_CELL_LIMIT_EXCEEDED',
      row: 1,
      column: 'sample',
    });
  });
});
