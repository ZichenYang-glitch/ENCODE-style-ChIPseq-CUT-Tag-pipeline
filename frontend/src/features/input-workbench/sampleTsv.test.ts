import { describe, expect, it } from 'vitest';
import { readWorkbenchSchema } from './schemaContract';
import {
  parseSampleTsv,
  type DraftSampleRow,
} from './sampleTsv';
import { createAuthoringSchemaFixture } from './test-fixtures';

function contract() {
  const parsed = readWorkbenchSchema(createAuthoringSchemaFixture());
  if (!parsed.ok) throw new Error('fixture must be supported');
  return parsed.value;
}

describe('parseSampleTsv', () => {
  it('accepts quoted fields, empty optional cells, and CRLF without coercing strings', () => {
    const result = parseSampleTsv(
      'sample\tfastq_1\tlayout\tfastq_2\r\n' +
        'S1\t/data/S1.fastq.gz\tSE\t\r\n',
      contract(),
      () => 'row-1',
    );
    expect(result).toEqual({
      ok: true,
      rows: [
        {
          id: 'row-1',
          values: {
            sample: 'S1',
            fastq_1: '/data/S1.fastq.gz',
            layout: 'SE',
            fastq_2: '',
          },
        },
      ],
    });
  });

  it.each([
    ['empty file', ''],
    ['header only', 'sample\tfastq_1\tlayout\n'],
    ['duplicate header', 'sample\tfastq_1\tlayout\tsample\nS1\t/a\tSE\tS2\n'],
    ['missing required header', 'sample\tfastq_1\nS1\t/a\n'],
    ['unknown header', 'sample\tfastq_1\tlayout\tunknown\nS1\t/a\tSE\tx\n'],
    ['ragged row', 'sample\tfastq_1\tlayout\nS1\t/a\n'],
    ['embedded newline', 'sample\tfastq_1\tlayout\nS1\t"/a\n/b"\tSE\n'],
  ])('rejects %s', (_label, text) => {
    const result = parseSampleTsv(text, contract(), () => 'unused');
    expect(result.ok).toBe(false);
  });

  it('enforces Unicode code-point cell limits without replacing existing rows', () => {
    const schema = createAuthoringSchemaFixture();
    schema.limits.max_sample_cell_length = 2;
    const parsed = readWorkbenchSchema(schema);
    expect(parsed.ok).toBe(true);
    if (!parsed.ok) return;
    const existing: DraftSampleRow[] = [
      { id: 'existing', values: { sample: 'old' } },
    ];
    const result = parseSampleTsv(
      'sample\tfastq_1\tlayout\nS1\t😀😀😀\tSE\n',
      parsed.value,
      () => 'new',
    );
    expect(result.ok).toBe(false);
    expect(existing).toEqual([{ id: 'existing', values: { sample: 'old' } }]);
  });

  it.each([
    [
      'file bytes',
      (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
        schema.limits.max_request_bytes = 4;
      },
      'sample\tfastq_1\tlayout\nS1\t/a\tSE\n',
    ],
    [
      'row count',
      (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
        schema.limits.max_sample_rows = 1;
      },
      'sample\tfastq_1\tlayout\nS1\t/a\tSE\nS2\t/b\tSE\n',
    ],
    [
      'column count',
      (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
        schema.limits.max_sample_columns = 2;
      },
      'sample\tfastq_1\tlayout\nS1\t/a\tSE\n',
    ],
    [
      'header length',
      (schema: ReturnType<typeof createAuthoringSchemaFixture>) => {
        schema.limits.max_sample_column_name_length = 5;
      },
      'sample\tfastq_1\tlayout\nS1\t/a\tSE\n',
    ],
  ])('enforces the advertised %s ceiling', (_label, mutate, text) => {
    const raw = createAuthoringSchemaFixture();
    mutate(raw);
    const parsed = readWorkbenchSchema(raw);
    expect(parsed.ok).toBe(true);
    if (!parsed.ok) return;
    expect(parseSampleTsv(text, parsed.value, () => 'row')).toMatchObject({
      ok: false,
      issue: { code: 'SAMPLE_TSV_LIMIT_EXCEEDED' },
    });
  });
});
