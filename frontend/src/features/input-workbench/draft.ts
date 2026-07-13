import type {
  JsonValueInput,
  ValidationRequest,
  ValidationRequestConfig,
  ValidationRequestOptions,
  ValidationRequestSamplesAnyOfItem,
} from '../../api/generated/models';
import { isJsonObject, isJsonValue } from './jsonSafety';
import type { DraftSampleRow } from './sampleTsv';

export type DraftReviewResult =
  | {
      ok: true;
      payload: ValidationRequest & {
        samples: ValidationRequestSamplesAnyOfItem[];
        options: ValidationRequestOptions;
      };
      serialized: string;
      byteSize: number;
    }
  | {
      ok: false;
      code: 'DRAFT_NOT_JSON_SAFE' | 'DRAFT_TOO_LARGE';
      message: string;
    };

function stableJsonValue(value: JsonValueInput): JsonValueInput {
  if (Array.isArray(value)) return value.map(stableJsonValue);
  if (value !== null && typeof value === 'object') {
    return Object.fromEntries(
      Object.keys(value)
        .sort((left, right) => left.localeCompare(right))
        .map((key) => [key, stableJsonValue(value[key])]),
    );
  }
  return value;
}

export function buildDraftReview(
  config: ValidationRequestConfig,
  rows: DraftSampleRow[],
  options: ValidationRequestOptions,
  sampleColumnOrder: string[],
  maxRequestBytes: number,
): DraftReviewResult {
  if (!isJsonObject(config) || !isJsonObject(options)) {
    return {
      ok: false,
      code: 'DRAFT_NOT_JSON_SAFE',
      message: 'The draft contains a value that cannot be sent as JSON.',
    };
  }
  if (
    rows.some((row) =>
      Object.values(row.values).some((value) => typeof value !== 'string'),
    )
  ) {
    return {
      ok: false,
      code: 'DRAFT_NOT_JSON_SAFE',
      message: 'Every sample cell must remain a string.',
    };
  }

  const samples = rows.map((row) =>
    Object.fromEntries(
      sampleColumnOrder.map((key) => [key, row.values[key] ?? '']),
    ),
  );
  const payload = {
    config: stableJsonValue(config) as ValidationRequestConfig,
    samples,
    options: stableJsonValue(options) as ValidationRequestOptions,
  };
  if (!isJsonValue(payload)) {
    return {
      ok: false,
      code: 'DRAFT_NOT_JSON_SAFE',
      message: 'The draft contains a value that cannot be sent as JSON.',
    };
  }

  const serialized = JSON.stringify(payload, null, 2);
  const byteSize = new TextEncoder().encode(serialized).byteLength;
  if (byteSize > maxRequestBytes) {
    return {
      ok: false,
      code: 'DRAFT_TOO_LARGE',
      message: 'The draft is larger than the workflow authoring limit.',
    };
  }
  return { ok: true, payload, serialized, byteSize };
}
