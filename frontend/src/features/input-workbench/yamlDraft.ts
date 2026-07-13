import { parseAllDocuments, stringify } from 'yaml';
import type { ValidationRequestConfig } from '../../api/generated/models';
import { compareJsonKeys, isJsonObject } from './jsonSafety';

const MAX_YAML_ALIASES = 100;
const INVALID_YAML_VALUE = Symbol('invalid-yaml-value');

export interface YamlIssue {
  code: 'CONFIG_YAML_INVALID';
  message: string;
  line: number | null;
  column: number | null;
}

export type ConfigYamlResult =
  | { ok: true; value: ValidationRequestConfig }
  | { ok: false; issue: YamlIssue };

function issue(line: number | null = null, column: number | null = null): ConfigYamlResult {
  return {
    ok: false,
    issue: {
      code: 'CONFIG_YAML_INVALID',
      message:
        line === null
          ? 'Config YAML must contain one JSON-compatible object.'
          : `Config YAML is invalid near line ${line}, column ${column ?? 1}.`,
      line,
      column,
    },
  };
}

function yamlValueToJson(
  value: unknown,
  seen: Set<object> = new Set(),
): unknown | typeof INVALID_YAML_VALUE {
  if (value === null || typeof value === 'string' || typeof value === 'boolean') {
    return value;
  }
  if (typeof value === 'bigint') {
    if (
      value < BigInt(Number.MIN_SAFE_INTEGER) ||
      value > BigInt(Number.MAX_SAFE_INTEGER)
    ) {
      return INVALID_YAML_VALUE;
    }
    return Number(value);
  }
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) return INVALID_YAML_VALUE;
    if (Number.isInteger(value) && !Number.isSafeInteger(value)) {
      return INVALID_YAML_VALUE;
    }
    return value;
  }
  if (typeof value !== 'object' || seen.has(value)) {
    return INVALID_YAML_VALUE;
  }

  seen.add(value);
  if (Array.isArray(value)) {
    const result: unknown[] = [];
    for (const item of value) {
      const converted = yamlValueToJson(item, seen);
      if (converted === INVALID_YAML_VALUE) {
        seen.delete(value);
        return INVALID_YAML_VALUE;
      }
      result.push(converted);
    }
    seen.delete(value);
    return result;
  }
  if (value instanceof Map) {
    const result: Record<string, unknown> = Object.create(null);
    for (const [key, item] of value) {
      if (
        typeof key !== 'string' ||
        Object.prototype.hasOwnProperty.call(result, key)
      ) {
        seen.delete(value);
        return INVALID_YAML_VALUE;
      }
      const converted = yamlValueToJson(item, seen);
      if (converted === INVALID_YAML_VALUE) {
        seen.delete(value);
        return INVALID_YAML_VALUE;
      }
      result[key] = converted;
    }
    seen.delete(value);
    return result;
  }
  seen.delete(value);
  return INVALID_YAML_VALUE;
}

export function parseConfigYaml(text: string): ConfigYamlResult {
  try {
    const documents = parseAllDocuments(text, {
      schema: 'core',
      strict: true,
      uniqueKeys: true,
      intAsBigInt: true,
      prettyErrors: false,
    });
    if (documents.length !== 1) return issue();
    const document = documents[0];
    const problem = document.errors[0] ?? document.warnings[0];
    if (problem) {
      const position = problem.linePos?.[0];
      return issue(position?.line ?? null, position?.col ?? null);
    }
    const yamlValue: unknown = document.toJS({
      mapAsMap: true,
      maxAliasCount: MAX_YAML_ALIASES,
    });
    const value = yamlValueToJson(yamlValue);
    return isJsonObject(value) ? { ok: true, value } : issue();
  } catch {
    return issue();
  }
}

export function stringifyConfigYaml(config: ValidationRequestConfig): string {
  return stringify(config, {
    aliasDuplicateObjects: false,
    lineWidth: 100,
    sortMapEntries: (left, right) =>
      compareJsonKeys(String(left.key), String(right.key)),
  });
}
