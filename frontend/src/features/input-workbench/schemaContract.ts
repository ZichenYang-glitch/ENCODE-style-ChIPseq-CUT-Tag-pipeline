import Ajv2020 from 'ajv/dist/2020';
import { getDefaultFormState, type RJSFSchema } from '@rjsf/utils';
import { customizeValidator } from '@rjsf/validator-ajv8';
import type {
  ValidationRequestConfig,
  WorkflowInputLimitsResponse,
  WorkflowSchemaResponse,
} from '../../api/generated/models';
import { isJsonObject, isPlainObject } from './jsonSafety';
import {
  isSampleCellSafe,
  isSampleColumnNameSafe,
} from './sampleValidation';

const SUPPORTED_SCHEMA_VERSION = '1.1.0';
const SUPPORTED_SCHEMA_DIALECT =
  'https://json-schema.org/draft/2020-12/schema';

const AUTHORING_LIMITS_1_1_0: WorkflowInputLimitsResponse = {
  max_request_bytes: 2_097_152,
  max_sample_rows: 1_000,
  max_sample_columns: 64,
  max_sample_column_name_length: 128,
  max_sample_cell_length: 4_096,
};

export const rjsfValidator = customizeValidator({
  AjvClass: Ajv2020,
  ajvOptionsOverrides: {
    allErrors: true,
    strict: false,
  },
});

export interface SampleColumn {
  key: string;
  label: string;
  description: string | null;
  required: boolean;
  enumValues: string[] | null;
}

export interface WorkbenchSchema {
  contract: WorkflowSchemaResponse;
  configSchema: RJSFSchema;
  sampleSchema: RJSFSchema;
  optionSchema: RJSFSchema;
  sampleColumns: SampleColumn[];
  limits: WorkflowInputLimitsResponse;
}

export type WorkbenchSchemaResult =
  | { ok: true; value: WorkbenchSchema }
  | {
      ok: false;
      code: 'AUTHORING_SCHEMA_UNSUPPORTED';
      message: string;
    };

function unsupported(): WorkbenchSchemaResult {
  return {
    ok: false,
    code: 'AUTHORING_SCHEMA_UNSUPPORTED',
    message: 'This workflow authoring contract is not supported by this client.',
  };
}

function hasMode(values: string[], expected: string): boolean {
  return values.includes(expected);
}

function hasSupportedLimits(limits: WorkflowInputLimitsResponse): boolean {
  return (
    limits.max_request_bytes === AUTHORING_LIMITS_1_1_0.max_request_bytes &&
    limits.max_sample_rows === AUTHORING_LIMITS_1_1_0.max_sample_rows &&
    limits.max_sample_columns === AUTHORING_LIMITS_1_1_0.max_sample_columns &&
    limits.max_sample_column_name_length ===
      AUTHORING_LIMITS_1_1_0.max_sample_column_name_length &&
    limits.max_sample_cell_length ===
      AUTHORING_LIMITS_1_1_0.max_sample_cell_length
  );
}

function asObjectSchema(value: unknown): RJSFSchema | null {
  if (!isPlainObject(value) || value.type !== 'object') return null;
  if (value.$schema !== SUPPORTED_SCHEMA_DIALECT) return null;
  return value as RJSFSchema;
}

function sampleSchemaParts(value: unknown): {
  schema: RJSFSchema;
  properties: Record<string, unknown>;
  required: Set<string>;
} | null {
  if (!isPlainObject(value) || value.type !== 'array') return null;
  if (value.$schema !== SUPPORTED_SCHEMA_DIALECT) return null;
  if (!isPlainObject(value.items) || value.items.type !== 'object') return null;
  if (!isPlainObject(value.items.properties)) return null;
  if (value.items.additionalProperties !== false) return null;
  const requiredValues = value.items.required;
  if (
    requiredValues !== undefined &&
    (!Array.isArray(requiredValues) ||
      !requiredValues.every((item) => typeof item === 'string'))
  ) {
    return null;
  }
  const propertyKeys = new Set(Object.keys(value.items.properties));
  if ((requiredValues ?? []).some((item) => !propertyKeys.has(item as string))) {
    return null;
  }
  return {
    schema: value as RJSFSchema,
    properties: value.items.properties,
    required: new Set((requiredValues ?? []) as string[]),
  };
}

function readSampleColumns(
  properties: Record<string, unknown>,
  required: Set<string>,
  limits: WorkflowInputLimitsResponse,
): SampleColumn[] | null {
  if (Object.keys(properties).length > limits.max_sample_columns) return null;
  const columns: SampleColumn[] = [];
  for (const [key, rawDefinition] of Object.entries(properties)) {
    if (!isSampleColumnNameSafe(key, limits.max_sample_column_name_length)) {
      return null;
    }
    if (!isPlainObject(rawDefinition) || rawDefinition.type !== 'string') {
      return null;
    }
    const rawEnum = rawDefinition.enum;
    if (
      rawEnum !== undefined &&
      (!Array.isArray(rawEnum) ||
        !rawEnum.every(
          (item) =>
            typeof item === 'string' &&
            isSampleCellSafe(item, limits.max_sample_cell_length),
        ))
    ) {
      return null;
    }
    columns.push({
      key,
      label:
        typeof rawDefinition.title === 'string' ? rawDefinition.title : key,
      description:
        typeof rawDefinition.description === 'string'
          ? rawDefinition.description
          : null,
      required: required.has(key),
      enumValues: rawEnum === undefined ? null : [...rawEnum],
    });
  }
  return columns.length > 0 ? columns : null;
}

export function readWorkbenchSchema(
  contract: WorkflowSchemaResponse,
): WorkbenchSchemaResult {
  if (
    contract.schema_version !== SUPPORTED_SCHEMA_VERSION ||
    contract.schema_dialect !== SUPPORTED_SCHEMA_DIALECT ||
    !hasMode(contract.authoring_modes.config, 'schema_form') ||
    !hasMode(contract.authoring_modes.config, 'yaml') ||
    !hasMode(contract.authoring_modes.samples, 'tsv_upload') ||
    !hasMode(contract.authoring_modes.samples, 'inline_table') ||
    !hasMode(contract.authoring_modes.options, 'schema_form') ||
    !hasMode(contract.input_modes.config, 'object') ||
    !hasMode(contract.input_modes.samples, 'inline_rows') ||
    !hasMode(contract.input_modes.options, 'object') ||
    !hasSupportedLimits(contract.limits)
  ) {
    return unsupported();
  }

  const configSchema = asObjectSchema(contract.config_schema);
  const optionSchema = asObjectSchema(contract.option_schema);
  const sampleParts = sampleSchemaParts(contract.sample_schema);
  if (!configSchema || !optionSchema || !sampleParts) return unsupported();

  const sampleColumns = readSampleColumns(
    sampleParts.properties,
    sampleParts.required,
    contract.limits,
  );
  if (!sampleColumns) return unsupported();

  try {
    createDefaultObject(configSchema);
    createDefaultObject(optionSchema);
    rjsfValidator.isValid(sampleParts.schema, [], sampleParts.schema);
  } catch {
    return unsupported();
  }

  return {
    ok: true,
    value: {
      contract,
      configSchema,
      sampleSchema: sampleParts.schema,
      optionSchema,
      sampleColumns,
      limits: contract.limits,
    },
  };
}

export function createDefaultObject(schema: RJSFSchema): ValidationRequestConfig {
  const defaults = getDefaultFormState(rjsfValidator, schema);
  if (defaults === undefined) return {};
  if (!isJsonObject(defaults)) {
    throw new Error('Schema defaults must produce one JSON object.');
  }
  return defaults;
}
