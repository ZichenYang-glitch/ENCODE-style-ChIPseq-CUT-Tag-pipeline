const FORBIDDEN_SAMPLE_TRANSPORT_CHARACTERS = /\0|\t|\r|\n/;

export interface SampleTransportIssue {
  code: 'SAMPLE_CELL_INVALID' | 'SAMPLE_CELL_LIMIT_EXCEEDED';
  message: string;
  row: number;
  column: string;
}

export function codePointLength(value: string): number {
  return Array.from(value).length;
}

export function isSampleColumnNameSafe(
  value: string,
  maxLength: number,
): boolean {
  return (
    value.length > 0 &&
    codePointLength(value) <= maxLength &&
    !FORBIDDEN_SAMPLE_TRANSPORT_CHARACTERS.test(value)
  );
}

function validateSampleCell(
  value: unknown,
  maxLength: number,
): Pick<SampleTransportIssue, 'code' | 'message'> | null {
  if (typeof value !== 'string') {
    return {
      code: 'SAMPLE_CELL_INVALID',
      message: 'Every sample cell must remain a string.',
    };
  }
  if (FORBIDDEN_SAMPLE_TRANSPORT_CHARACTERS.test(value)) {
    return {
      code: 'SAMPLE_CELL_INVALID',
      message: 'Sample cells cannot contain tabs, line breaks, or NUL characters.',
    };
  }
  if (codePointLength(value) > maxLength) {
    return {
      code: 'SAMPLE_CELL_LIMIT_EXCEEDED',
      message: 'A sample cell exceeds the workflow authoring limit.',
    };
  }
  return null;
}

export function isSampleCellSafe(value: unknown, maxLength: number): boolean {
  return validateSampleCell(value, maxLength) === null;
}

export function validateSampleRows(
  rows: readonly Readonly<Record<string, unknown>>[],
  columnOrder: readonly string[],
  maxCellLength: number,
  firstRowNumber = 1,
): SampleTransportIssue | null {
  for (const [rowIndex, row] of rows.entries()) {
    for (const column of columnOrder) {
      const value = Object.prototype.hasOwnProperty.call(row, column)
        ? row[column]
        : '';
      const issue = validateSampleCell(value, maxCellLength);
      if (issue !== null) {
        return {
          ...issue,
          row: firstRowNumber + rowIndex,
          column,
        };
      }
    }
  }
  return null;
}
