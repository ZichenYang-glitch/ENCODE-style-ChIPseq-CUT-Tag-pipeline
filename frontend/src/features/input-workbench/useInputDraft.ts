import { useCallback, useMemo, useReducer } from 'react';
import type {
  ValidationRequestConfig,
  ValidationRequestOptions,
} from '../../api/generated/models';
import { buildDraftReview } from './draft';
import { isJsonObject } from './jsonSafety';
import {
  createDefaultObject,
  rjsfValidator,
  type WorkbenchSchema,
} from './schemaContract';
import {
  createEmptySampleRow,
  type DraftSampleRow,
  type SampleImportIssue,
} from './sampleTsv';
import {
  parseConfigYaml,
  stringifyConfigYaml,
  type YamlIssue,
} from './yamlDraft';

interface InputDraftState {
  config: ValidationRequestConfig;
  configRevision: number;
  yamlText: string;
  yamlIssue: YamlIssue | null;
  yamlRevision: number;
  options: ValidationRequestOptions;
  rows: DraftSampleRow[];
  sampleImportIssue: SampleImportIssue | null;
  importedFileName: string | null;
}

type InputDraftAction =
  | { type: 'config-form'; value: ValidationRequestConfig }
  | {
      type: 'yaml-edit';
      text: string;
      parsed: ReturnType<typeof parseConfigYaml>;
    }
  | { type: 'yaml-sync'; text: string }
  | { type: 'options-form'; value: ValidationRequestOptions }
  | { type: 'samples-replace'; rows: DraftSampleRow[]; fileName: string }
  | { type: 'sample-import-failed'; issue: SampleImportIssue }
  | { type: 'sample-add'; row: DraftSampleRow }
  | { type: 'sample-remove'; rowId: string }
  | { type: 'sample-cell'; rowId: string; column: string; value: string };

function reducer(state: InputDraftState, action: InputDraftAction): InputDraftState {
  switch (action.type) {
    case 'config-form':
      return {
        ...state,
        config: action.value,
        configRevision: state.configRevision + 1,
      };
    case 'yaml-edit':
      if (!action.parsed.ok) {
        return {
          ...state,
          yamlText: action.text,
          yamlIssue: action.parsed.issue,
        };
      }
      return {
        ...state,
        config: action.parsed.value,
        configRevision: state.configRevision + 1,
        yamlText: action.text,
        yamlIssue: null,
        yamlRevision: state.configRevision + 1,
      };
    case 'yaml-sync':
      return {
        ...state,
        yamlText: action.text,
        yamlIssue: null,
        yamlRevision: state.configRevision,
      };
    case 'options-form':
      return { ...state, options: action.value };
    case 'samples-replace':
      return {
        ...state,
        rows: action.rows,
        importedFileName: action.fileName,
        sampleImportIssue: null,
      };
    case 'sample-import-failed':
      return { ...state, sampleImportIssue: action.issue };
    case 'sample-add':
      return {
        ...state,
        rows: [...state.rows, action.row],
        sampleImportIssue: null,
      };
    case 'sample-remove':
      return {
        ...state,
        rows: state.rows.filter((row) => row.id !== action.rowId),
      };
    case 'sample-cell':
      return {
        ...state,
        rows: state.rows.map((row) =>
          row.id === action.rowId
            ? { ...row, values: { ...row.values, [action.column]: action.value } }
            : row,
        ),
      };
  }
}

function createInitialState(schema: WorkbenchSchema): InputDraftState {
  const config = createDefaultObject(schema.configSchema);
  const options = createDefaultObject(schema.optionSchema);
  return {
    config,
    configRevision: 0,
    yamlText: stringifyConfigYaml(config),
    yamlIssue: null,
    yamlRevision: 0,
    options,
    rows: [],
    sampleImportIssue: null,
    importedFileName: null,
  };
}

export function useInputDraft(schema: WorkbenchSchema) {
  const [state, dispatch] = useReducer(reducer, schema, createInitialState);
  const sampleColumnOrder = useMemo(
    () => schema.sampleColumns.map((column) => column.key),
    [schema.sampleColumns],
  );
  const review = useMemo(
    () =>
      buildDraftReview(
        state.config,
        state.rows,
        state.options,
        sampleColumnOrder,
        schema.limits.max_request_bytes,
      ),
    [
      sampleColumnOrder,
      schema.limits.max_request_bytes,
      state.config,
      state.options,
      state.rows,
    ],
  );

  const configValid = rjsfValidator.isValid(
    schema.configSchema,
    state.config,
    schema.configSchema,
  );
  const optionsValid = rjsfValidator.isValid(
    schema.optionSchema,
    state.options,
    schema.optionSchema,
  );
  const sampleValues = state.rows.map((row) => row.values);
  const samplesValid = rjsfValidator.isValid(
    schema.sampleSchema,
    sampleValues,
    schema.sampleSchema,
  );

  const setConfig = useCallback((value: unknown) => {
    if (isJsonObject(value)) dispatch({ type: 'config-form', value });
  }, []);
  const setOptions = useCallback((value: unknown) => {
    if (isJsonObject(value)) dispatch({ type: 'options-form', value });
  }, []);
  const editYaml = useCallback((text: string) => {
    dispatch({ type: 'yaml-edit', text, parsed: parseConfigYaml(text) });
  }, []);
  const replaceSamples = useCallback(
    (rows: DraftSampleRow[], fileName: string) => {
      dispatch({ type: 'samples-replace', rows, fileName });
    },
    [],
  );
  const failSampleImport = useCallback((issue: SampleImportIssue) => {
    dispatch({ type: 'sample-import-failed', issue });
  }, []);
  const removeSample = useCallback((rowId: string) => {
    dispatch({ type: 'sample-remove', rowId });
  }, []);
  const updateSample = useCallback(
    (rowId: string, column: string, value: string) => {
      dispatch({ type: 'sample-cell', rowId, column, value });
    },
    [],
  );

  return {
    state,
    review,
    configValid,
    optionsValid,
    samplesValid,
    reviewReady:
      state.yamlIssue === null &&
      configValid &&
      optionsValid &&
      samplesValid &&
      review.ok,
    setConfig,
    setOptions,
    editYaml,
    syncYaml() {
      if (
        state.yamlIssue === null &&
        state.yamlRevision !== state.configRevision
      ) {
        dispatch({ type: 'yaml-sync', text: stringifyConfigYaml(state.config) });
      }
    },
    replaceSamples,
    failSampleImport,
    addSample() {
      if (state.rows.length < schema.limits.max_sample_rows) {
        dispatch({ type: 'sample-add', row: createEmptySampleRow(schema) });
      }
    },
    removeSample,
    updateSample,
  };
}

export type InputDraftController = ReturnType<typeof useInputDraft>;
