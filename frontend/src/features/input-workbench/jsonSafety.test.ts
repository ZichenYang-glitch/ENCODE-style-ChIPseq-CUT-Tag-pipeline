import { describe, expect, it } from 'vitest';
import { isJsonValue } from './jsonSafety';

describe('isJsonValue', () => {
  it('accepts the JSON-safe integer boundaries', () => {
    expect(isJsonValue(Number.MAX_SAFE_INTEGER)).toBe(true);
    expect(isJsonValue(Number.MIN_SAFE_INTEGER)).toBe(true);
  });

  it('fails closed for finite integers outside the JSON-safe range', () => {
    expect(isJsonValue(Number.MAX_SAFE_INTEGER + 1)).toBe(false);
    expect(isJsonValue(Number.MIN_SAFE_INTEGER - 1)).toBe(false);
    expect(isJsonValue({ nested: Number.MAX_SAFE_INTEGER + 1 })).toBe(false);
  });
});
