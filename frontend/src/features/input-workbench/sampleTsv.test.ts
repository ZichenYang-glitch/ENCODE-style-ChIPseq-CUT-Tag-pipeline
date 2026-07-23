import { describe, expect, it } from 'vitest';
import { readWorkbenchSchema } from './schemaContract';
import {
  createEmptySampleRow,
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
  it('prefills adapter-declared string constants in a new row', () => {
    const schema = createAuthoringSchemaFixture();
    const items = (schema.sample_schema as {
      items: {
        properties: Record<string, unknown>;
        required: string[];
      };
    }).items;
    items.properties.platform = { type: 'string', const: 'ILLUMINA' };
    items.required.push('platform');
    const parsed = readWorkbenchSchema(schema);
    if (!parsed.ok) throw new Error('fixture must be supported');

    expect(createEmptySampleRow(parsed.value, () => 'row-const')).toEqual({
      id: 'row-const',
      values: {
        sample: '',
        fastq_1: '',
        layout: '',
        fastq_2: '',
        platform: 'ILLUMINA',
      },
    });
  });

  it('rejects a TSV value that conflicts with an adapter-declared constant', () => {
    const schema = createAuthoringSchemaFixture();
    const items = (schema.sample_schema as {
      items: {
        properties: Record<string, unknown>;
        required: string[];
      };
    }).items;
    items.properties.platform = { type: 'string', const: 'ILLUMINA' };
    items.required.push('platform');
    const parsedSchema = readWorkbenchSchema(schema);
    if (!parsedSchema.ok) throw new Error('fixture must be supported');

    expect(
      parseSampleTsv(
        'sample\tfastq_1\tlayout\tfastq_2\tplatform\n' +
          'S1\t/data/S1.fastq.gz\tSE\t\tONT\n',
        parsedSchema.value,
      ),
    ).toEqual({
      ok: false,
      issue: {
        code: 'SAMPLE_TSV_INVALID',
        message: 'A sample cell does not match the adapter-declared constant.',
        row: 2,
        column: 'platform',
      },
    });
  });

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

  it('enforces the real Unicode code-point cell ceiling without replacing existing rows', () => {
    const schema = contract();
    const existing: DraftSampleRow[] = [
      { id: 'existing', values: { sample: 'old' } },
    ];
    const atLimit = parseSampleTsv(
      `sample\tfastq_1\tlayout\nS1\t${'😀'.repeat(4_096)}\tSE\n`,
      schema,
      () => 'accepted',
    );
    expect(atLimit.ok).toBe(true);
    const result = parseSampleTsv(
      `sample\tfastq_1\tlayout\nS1\t${'😀'.repeat(4_097)}\tSE\n`,
      schema,
      () => 'new',
    );
    expect(result.ok).toBe(false);
    expect(existing).toEqual([{ id: 'existing', values: { sample: 'old' } }]);
  });

  it.each([
    [
      'file bytes',
      () => 'x'.repeat(2_097_153),
    ],
    [
      'row count',
      () =>
        `sample\tfastq_1\tlayout\n${Array.from(
          { length: 1_001 },
          (_, index) => `S${index}\t/a\tSE`,
        ).join('\n')}\n`,
    ],
    [
      'column count',
      () =>
        `${Array.from({ length: 65 }, (_, index) => `c${index}`).join('\t')}\n${Array.from(
          { length: 65 },
          () => 'x',
        ).join('\t')}\n`,
    ],
    [
      'header length',
      () => `${'h'.repeat(129)}\nvalue\n`,
    ],
  ])('enforces the advertised %s ceiling', (_label, makeText) => {
    expect(parseSampleTsv(makeText(), contract(), () => 'row')).toMatchObject({
      ok: false,
      issue: { code: 'SAMPLE_TSV_LIMIT_EXCEEDED' },
    });
  });
});
