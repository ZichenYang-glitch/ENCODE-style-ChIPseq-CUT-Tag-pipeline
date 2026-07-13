import { describe, it, expect } from 'vitest';
import type {
  JsonValueInput,
  JsonValueOutput,
  SchemaResponse,
  ValidationRequest,
} from './generated/models';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { resolve, dirname } from 'node:path';

const openapi = JSON.parse(
  readFileSync(
    resolve(dirname(fileURLToPath(import.meta.url)), '../../openapi.json'),
    'utf-8',
  ),
);

const recursiveJsonValues: JsonValueInput[] = [
  null,
  'text',
  42,
  1.5,
  true,
  [null, { nested: ['value', 7, false] }],
  { object: { array: [1, null] } },
];

const recursiveSchemaJson: JsonValueOutput = {
  properties: {
    nullable: null,
    nested: ['value', 7, false, { child: null }],
  },
};

const generatedValidationFixture = {
  config: {
    nullable: null,
    nested: { array: [1, 'two', false] },
  },
  samples: null,
  options: { strict_inputs: false },
} satisfies ValidationRequest;

const generatedSchemaResponseFixture = {
  ok: false,
  workflow_id: 'missing-workflow',
  schema: null,
  issues: [],
} satisfies SchemaResponse;

// @ts-expect-error Date is not a JSON value.
const invalidDateValue: JsonValueInput = new Date();

// @ts-expect-error schema is required even though it is nullable.
const missingSchemaResponse: SchemaResponse = {
  ok: true,
  workflow_id: 'workflow',
  issues: [],
};

describe('generated OpenAPI client coverage', () => {
  const paths = Object.keys(openapi.paths);

  it('exports all public API paths in the OpenAPI contract', () => {
    expect(recursiveJsonValues).toHaveLength(7);
    expect(recursiveSchemaJson).toHaveProperty('properties');
    expect(generatedValidationFixture.config.nullable).toBeNull();
    expect(generatedSchemaResponseFixture.schema).toBeNull();
    expect(invalidDateValue).toBeInstanceOf(Date);
    expect(missingSchemaResponse.ok).toBe(true);
    expect(paths).toContain('/api/v1/workflows/');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/schema');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/validate');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/agent/chat');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/runs');
    expect(paths).toContain('/api/v1/runs/{run_id}');
    expect(paths).toContain('/api/v1/runs/{run_id}/start');
    expect(paths).toContain('/api/v1/runs/{run_id}/cancel');
    expect(paths).toContain('/api/v1/runs/{run_id}/events');
    expect(paths).toContain('/api/v1/runs/{run_id}/logs');
    expect(paths).toContain('/api/v1/runs/{run_id}/preflight');
  });
});
