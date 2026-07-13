import { describe, expect, it } from 'vitest';
import { parseConfigYaml, stringifyConfigYaml } from './yamlDraft';

describe('config YAML', () => {
  it('round-trips advanced nested fields without losing JSON values', () => {
    const source = {
      threads: 8,
      advanced: { experimental_mode: 'bounded', flags: [true, false] },
    };
    const parsed = parseConfigYaml(stringifyConfigYaml(source));
    expect(parsed).toEqual({ ok: true, value: source });
  });

  it('preserves the largest JSON-safe integer exactly', () => {
    expect(parseConfigYaml('advanced: 9007199254740991\n')).toEqual({
      ok: true,
      value: { advanced: 9_007_199_254_740_991 },
    });
    expect(parseConfigYaml('advanced: -9007199254740991\n')).toEqual({
      ok: true,
      value: { advanced: -9_007_199_254_740_991 },
    });
  });

  it('stringifies keys in locale-independent UTF-16 code-unit order', () => {
    const text = stringifyConfigYaml({ z: 1, 'é': 2, a: 3, _: 4, A: 5 });
    expect(
      text
        .trim()
        .split('\n')
        .map((line) => line.split(':', 1)[0]),
    ).toEqual(['A', '_', 'a', 'z', 'é']);
  });

  it.each([
    ['invalid syntax', 'threads: ['],
    ['duplicate key', 'threads: 8\nthreads: 4\n'],
    ['multiple documents', 'threads: 8\n---\nthreads: 4\n'],
    ['custom tag', 'value: !unsafe payload\n'],
    ['array root', '- one\n- two\n'],
    ['null root', 'null\n'],
    ['non-finite value', 'value: .inf\n'],
    ['cyclic alias', '&root { self: *root }\n'],
    ['non-string mapping key', '1: value\n'],
    ['nested non-string mapping key', 'advanced:\n  true: value\n'],
    ['colliding mapping keys', '1: first\n"1": second\n'],
    ['unsafe integer', 'advanced: 9007199254740993\n'],
    ['unsafe negative integer', 'advanced: -9007199254740993\n'],
  ])('rejects %s with a controlled issue', (_label, text) => {
    const result = parseConfigYaml(text);
    expect(result.ok).toBe(false);
    if (result.ok) return;
    expect(result.issue.code).toBe('CONFIG_YAML_INVALID');
    expect(result.issue.message).not.toMatch(/Error:|at \w+ \(/);
  });
});
