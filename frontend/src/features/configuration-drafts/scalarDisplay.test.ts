import { describe, expect, it } from 'vitest';
import { formatScalarValue } from './scalarDisplay';

describe('scalar suggestion values', () => {
  it('distinguishes missing, null, empty, zero, and false', () => {
    expect(formatScalarValue({ presence: 'missing' })).toBe('Not set');
    expect(formatScalarValue({ presence: 'null' })).toBe('Null');
    expect(formatScalarValue({ presence: 'present', value: '' })).toBe('Empty');
    expect(formatScalarValue({ presence: 'present', value: 0 })).toBe('0');
    expect(formatScalarValue({ presence: 'present', value: false })).toBe('False');
    expect(formatScalarValue({ presence: 'present', value: true })).toBe('True');
    expect(formatScalarValue({ presence: 'present', value: 'kept' })).toBe('kept');
  });
});
