import Papa from 'papaparse';
import type { ValidationRequestSamplesAnyOfItem } from '../../api/generated/models';
import type { WorkbenchSchema } from './schemaContract';
import {
  isSampleColumnNameSafe,
  validateSampleRows,
} from './sampleValidation';

export interface DraftSampleRow {
  id: string;
  values: ValidationRequestSamplesAnyOfItem;
}

export interface SampleImportIssue {
  code:
    | 'SAMPLE_TSV_EMPTY'
    | 'SAMPLE_TSV_INVALID'
    | 'SAMPLE_TSV_LIMIT_EXCEEDED';
  message: string;
  row: number | null;
  column: string | null;
}

export type SampleTsvResult =
  | { ok: true; rows: DraftSampleRow[] }
  | { ok: false; issue: SampleImportIssue };

function failure(
  code: SampleImportIssue['code'],
  message: string,
  row: number | null = null,
  column: string | null = null,
): SampleTsvResult {
  return { ok: false, issue: { code, message, row, column } };
}

function isBlankRecord(row: string[]): boolean {
  return row.length === 1 && row[0] === '';
}

export function createEmptySampleRow(
  schema: WorkbenchSchema,
  createId: () => string = () => crypto.randomUUID(),
): DraftSampleRow {
  return {
    id: createId(),
    values: Object.fromEntries(schema.sampleColumns.map((column) => [column.key, ''])),
  };
}

export function parseSampleTsv(
  text: string,
  schema: WorkbenchSchema,
  createId: () => string = () => crypto.randomUUID(),
): SampleTsvResult {
  if (new TextEncoder().encode(text).byteLength > schema.limits.max_request_bytes) {
    return failure(
      'SAMPLE_TSV_LIMIT_EXCEEDED',
      'The sample file is larger than the workflow authoring limit.',
    );
  }

  const parsed = Papa.parse<string[]>(text, {
    delimiter: '\t',
    dynamicTyping: false,
    skipEmptyLines: false,
  });
  if (parsed.errors.length > 0) {
    return failure(
      'SAMPLE_TSV_INVALID',
      'The sample file is not valid tab-separated text.',
      parsed.errors[0]?.row === undefined ? null : parsed.errors[0].row + 1,
    );
  }

  const records = parsed.data.map((row) => [...row]);
  while (records.length > 0 && isBlankRecord(records[records.length - 1])) {
    records.pop();
  }
  if (records.length < 2) {
    return failure(
      'SAMPLE_TSV_EMPTY',
      'The sample file must contain a header and at least one data row.',
    );
  }

  const headers = records[0];
  if (
    headers.length === 0 ||
    headers.length > schema.limits.max_sample_columns
  ) {
    return failure(
      'SAMPLE_TSV_LIMIT_EXCEEDED',
      'The sample file has an unsupported number of columns.',
    );
  }
  if (
    headers.some(
      (header) =>
        !isSampleColumnNameSafe(
          header,
          schema.limits.max_sample_column_name_length,
        ),
    )
  ) {
    return failure(
      'SAMPLE_TSV_LIMIT_EXCEEDED',
      'A sample column name is invalid or exceeds the workflow authoring limit.',
    );
  }
  if (new Set(headers).size !== headers.length) {
    return failure('SAMPLE_TSV_INVALID', 'Sample column names must be unique.');
  }

  const declaredKeys = schema.sampleColumns.map((column) => column.key);
  const declaredKeySet = new Set(declaredKeys);
  const unknownHeader = headers.find((header) => !declaredKeySet.has(header));
  if (unknownHeader) {
    return failure(
      'SAMPLE_TSV_INVALID',
      'The sample file contains a column not declared by this workflow.',
      1,
      unknownHeader,
    );
  }
  const missingRequired = schema.sampleColumns.find(
    (column) => column.required && !headers.includes(column.key),
  );
  if (missingRequired) {
    return failure(
      'SAMPLE_TSV_INVALID',
      'The sample file is missing a required workflow column.',
      1,
      missingRequired.key,
    );
  }

  const dataRows = records.slice(1);
  if (dataRows.length > schema.limits.max_sample_rows) {
    return failure(
      'SAMPLE_TSV_LIMIT_EXCEEDED',
      'The sample file has more rows than the workflow authoring limit.',
    );
  }

  const headerIndexes = new Map(headers.map((header, index) => [header, index]));
  const rows: DraftSampleRow[] = [];
  for (const [index, rawRow] of dataRows.entries()) {
    const displayRow = index + 2;
    if (rawRow.length !== headers.length || isBlankRecord(rawRow)) {
      return failure(
        'SAMPLE_TSV_INVALID',
        'Every sample row must contain the same number of columns as the header.',
        displayRow,
      );
    }
    const values: ValidationRequestSamplesAnyOfItem = {};
    for (const column of schema.sampleColumns) {
      const headerIndex = headerIndexes.get(column.key);
      const value = headerIndex === undefined ? '' : rawRow[headerIndex];
      values[column.key] = value;
    }
    rows.push({ id: createId(), values });
  }
  const transportIssue = validateSampleRows(
    rows.map((row) => row.values),
    declaredKeys,
    schema.limits.max_sample_cell_length,
    2,
  );
  if (transportIssue !== null) {
    return failure(
      transportIssue.code === 'SAMPLE_CELL_LIMIT_EXCEEDED'
        ? 'SAMPLE_TSV_LIMIT_EXCEEDED'
        : 'SAMPLE_TSV_INVALID',
      transportIssue.message,
      transportIssue.row,
      transportIssue.column,
    );
  }
  return { ok: true, rows };
}
