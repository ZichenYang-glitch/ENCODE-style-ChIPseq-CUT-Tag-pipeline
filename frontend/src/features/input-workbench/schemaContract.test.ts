import { describe, expect, it } from 'vitest';
import { createAuthoringSchemaFixture } from './test-fixtures';
import {
  createDefaultObject,
  readWorkbenchSchema,
  rjsfValidator,
} from './schemaContract';

describe('readWorkbenchSchema', () => {
  it('accepts the PR135 contract and compiles all three Draft 2020-12 schemas', () => {
    const result = readWorkbenchSchema(createAuthoringSchemaFixture());
    expect(result.ok).toBe(true);
    if (!result.ok) return;

    const config = createDefaultObject(result.value.configSchema);
    const options = createDefaultObject(result.value.optionSchema);
    expect(config).toMatchObject({ outdir: 'results', threads: 8 });
    expect(options).toEqual({ strict_inputs: false });
    expect(
      rjsfValidator.isValid(result.value.configSchema, config, result.value.configSchema),
    ).toBe(true);
    expect(
      rjsfValidator.isValid(result.value.sampleSchema, [
        { sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' },
      ], result.value.sampleSchema),
    ).toBe(true);
    expect(
      rjsfValidator.isValid(result.value.optionSchema, options, result.value.optionSchema),
    ).toBe(true);
  });

  it.each([
    ['schema version', (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
      schema.schema_version = '2.0.0';
    }],
    ['schema dialect', (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
      Object.assign(schema, { schema_dialect: 'http://json-schema.org/draft-07/schema#' });
    }],
    ['required authoring mode', (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
      schema.authoring_modes.config = ['schema_form'];
    }],
    ['sample object root', (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
      schema.sample_schema = { type: 'object' };
    }],
  ])('fails closed for an unsupported %s', (_label, mutate) => {
    const schema = createAuthoringSchemaFixture();
    mutate(schema);
    const result = readWorkbenchSchema(schema);
    expect(result).toMatchObject({ ok: false, code: 'AUTHORING_SCHEMA_UNSUPPORTED' });
  });

  it.each([
    ['max_request_bytes', 2_097_152],
    ['max_sample_rows', 1_000],
    ['max_sample_columns', 64],
    ['max_sample_column_name_length', 128],
    ['max_sample_cell_length', 4_096],
  ] as const)(
    'fails closed when 1.0.0 %s differs from the platform ceiling',
    (field, ceiling) => {
      for (const value of [ceiling - 1, ceiling + 1]) {
        const schema = createAuthoringSchemaFixture();
        schema.limits[field] = value;
        expect(readWorkbenchSchema(schema)).toMatchObject({
          ok: false,
          code: 'AUTHORING_SCHEMA_UNSUPPORTED',
        });
      }
    },
  );

  it('fails closed when the sample schema declares 65 columns', () => {
    const schema = createAuthoringSchemaFixture();
    const items = (schema.sample_schema as {
      items: { properties: Record<string, unknown>; required: string[] };
    }).items;
    items.properties = Object.fromEntries(
      Array.from({ length: 65 }, (_, index) => [
        `column_${index}`,
        { type: 'string' },
      ]),
    );
    items.required = [];

    expect(readWorkbenchSchema(schema)).toMatchObject({
      ok: false,
      code: 'AUTHORING_SCHEMA_UNSUPPORTED',
    });
  });

  it.each([
    ['empty', ''],
    ['129 code points', 'c'.repeat(129)],
    ['NUL', 'bad\0column'],
    ['tab', 'bad\tcolumn'],
    ['carriage return', 'bad\rcolumn'],
    ['line feed', 'bad\ncolumn'],
  ])('fails closed for a %s sample property key', (_label, propertyKey) => {
    const schema = createAuthoringSchemaFixture();
    const items = (schema.sample_schema as {
      items: { properties: Record<string, unknown>; required: string[] };
    }).items;
    items.properties = { [propertyKey]: { type: 'string' } };
    items.required = [];

    expect(readWorkbenchSchema(schema)).toMatchObject({
      ok: false,
      code: 'AUTHORING_SCHEMA_UNSUPPORTED',
    });
  });

  it('requires every required sample key to be declared in properties', () => {
    const schema = createAuthoringSchemaFixture();
    const items = (schema.sample_schema as {
      items: { properties: Record<string, unknown>; required: string[] };
    }).items;
    items.required = [...items.required, 'undeclared'];

    expect(readWorkbenchSchema(schema)).toMatchObject({
      ok: false,
      code: 'AUTHORING_SCHEMA_UNSUPPORTED',
    });
  });

  it.each([
    ['overlong', 'x'.repeat(4_097)],
    ['NUL', 'bad\0value'],
    ['tab', 'bad\tvalue'],
    ['carriage return', 'bad\rvalue'],
    ['line feed', 'bad\nvalue'],
  ])('fails closed for a %s sample enum string', (_label, enumValue) => {
    const schema = createAuthoringSchemaFixture();
    const items = (schema.sample_schema as {
      items: { properties: Record<string, { enum?: string[] }> };
    }).items;
    items.properties.layout.enum = ['SE', enumValue];

    expect(readWorkbenchSchema(schema)).toMatchObject({
      ok: false,
      code: 'AUTHORING_SCHEMA_UNSUPPORTED',
    });
  });
});
