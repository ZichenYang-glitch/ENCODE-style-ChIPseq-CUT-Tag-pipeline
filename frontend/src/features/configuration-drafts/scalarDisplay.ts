import type { ScalarFieldValue } from './model';

export function formatScalarValue(value: ScalarFieldValue): string {
  if (value.presence === 'missing') return 'Not set';
  if (value.presence === 'null') return 'Null';
  if (typeof value.value === 'string') {
    return value.value === '' ? 'Empty' : value.value;
  }
  if (value.value === false) return 'False';
  if (value.value === true) return 'True';
  if (typeof value.value === 'number') return String(value.value);
  return 'Not set';
}
