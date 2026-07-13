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

  it.each([
    ['invalid syntax', 'threads: ['],
    ['duplicate key', 'threads: 8\nthreads: 4\n'],
    ['multiple documents', 'threads: 8\n---\nthreads: 4\n'],
    ['custom tag', 'value: !unsafe payload\n'],
    ['array root', '- one\n- two\n'],
    ['null root', 'null\n'],
    ['non-finite value', 'value: .inf\n'],
    ['cyclic alias', '&root { self: *root }\n'],
  ])('rejects %s with a controlled issue', (_label, text) => {
    const result = parseConfigYaml(text);
    expect(result.ok).toBe(false);
    if (result.ok) return;
    expect(result.issue.code).toBe('CONFIG_YAML_INVALID');
    expect(result.issue.message).not.toMatch(/Error:|at \w+ \(/);
  });
});
