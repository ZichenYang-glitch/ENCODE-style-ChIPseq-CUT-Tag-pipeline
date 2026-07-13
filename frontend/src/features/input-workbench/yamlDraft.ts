import { parseAllDocuments, stringify } from 'yaml';
import type { ValidationRequestConfig } from '../../api/generated/models';
import { isJsonObject } from './jsonSafety';

const MAX_YAML_ALIASES = 100;

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

export function parseConfigYaml(text: string): ConfigYamlResult {
  try {
    const documents = parseAllDocuments(text, {
      schema: 'core',
      strict: true,
      uniqueKeys: true,
      prettyErrors: false,
    });
    if (documents.length !== 1) return issue();
    const document = documents[0];
    const problem = document.errors[0] ?? document.warnings[0];
    if (problem) {
      const position = problem.linePos?.[0];
      return issue(position?.line ?? null, position?.col ?? null);
    }
    const value: unknown = document.toJS({
      mapAsMap: false,
      maxAliasCount: MAX_YAML_ALIASES,
    });
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
      String(left.key).localeCompare(String(right.key)),
  });
}
