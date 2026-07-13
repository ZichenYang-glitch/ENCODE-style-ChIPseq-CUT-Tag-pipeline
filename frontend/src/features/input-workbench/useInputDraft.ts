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
import { validateSampleRows } from './sampleValidation';
import {
  parseConfigYaml,
  stringifyConfigYaml,
  type YamlIssue,
} from './yamlDraft';

interface InputDraftState {
  config: ValidationRequestConfig;
  configRevision: number;
  configFormResetRevision: number;
  yamlText: string;
  yamlIssue: YamlIssue | null;
  yamlRevision: number;
  configFormIssue: FormSafetyIssue | null;
  options: ValidationRequestOptions;
  optionsFormResetRevision: number;
  optionsFormIssue: FormSafetyIssue | null;
  rows: DraftSampleRow[];
  sampleImportIssue: SampleImportIssue | null;
  importedFileName: string | null;
}

interface FormSafetyIssue {
  code: 'FORM_VALUE_UNSAFE';
  message: string;
}

function formSafetyIssue(section: 'Config' | 'Options'): FormSafetyIssue {
  return {
    code: 'FORM_VALUE_UNSAFE',
    message: `${section} form data cannot be represented safely as JSON. The previous safe value was kept.`,
  };
}

type InputDraftAction =
  | { type: 'config-form'; value: unknown }
  | { type: 'config-form-fallback-accepted' }
  | {
      type: 'yaml-edit';
      text: string;
      parsed: ReturnType<typeof parseConfigYaml>;
    }
  | { type: 'yaml-sync'; text: string }
  | { type: 'options-form'; value: unknown }
  | { type: 'options-form-fallback-accepted' }
  | { type: 'samples-replace'; rows: DraftSampleRow[]; fileName: string }
  | { type: 'sample-import-failed'; issue: SampleImportIssue }
  | { type: 'sample-add'; row: DraftSampleRow }
  | { type: 'sample-remove'; rowId: string }
  | { type: 'sample-cell'; rowId: string; column: string; value: string };

function reducer(state: InputDraftState, action: InputDraftAction): InputDraftState {
  switch (action.type) {
    case 'config-form':
      if (!isJsonObject(action.value)) {
        return {
          ...state,
          configFormIssue: formSafetyIssue('Config'),
          configFormResetRevision: state.configFormResetRevision + 1,
        };
      }
      return {
        ...state,
        config: action.value,
        configRevision: state.configRevision + 1,
        configFormIssue: null,
      };
    case 'config-form-fallback-accepted':
      return { ...state, configFormIssue: null };
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
        configFormIssue: null,
      };
    case 'yaml-sync':
      return {
        ...state,
        yamlText: action.text,
        yamlIssue: null,
        yamlRevision: state.configRevision,
      };
    case 'options-form':
      if (!isJsonObject(action.value)) {
        return {
          ...state,
          optionsFormIssue: formSafetyIssue('Options'),
          optionsFormResetRevision: state.optionsFormResetRevision + 1,
        };
      }
      return { ...state, options: action.value, optionsFormIssue: null };
    case 'options-form-fallback-accepted':
      return { ...state, optionsFormIssue: null };
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
    configFormResetRevision: 0,
    yamlText: stringifyConfigYaml(config),
    yamlIssue: null,
    yamlRevision: 0,
    configFormIssue: null,
    options,
    optionsFormResetRevision: 0,
    optionsFormIssue: null,
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
        schema.limits.max_sample_cell_length,
      ),
    [
      sampleColumnOrder,
      schema.limits.max_request_bytes,
      schema.limits.max_sample_cell_length,
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
  const sampleValues = useMemo(
    () => state.rows.map((row) => row.values),
    [state.rows],
  );
  const sampleTransportIssue = useMemo(
    () =>
      validateSampleRows(
        sampleValues,
        sampleColumnOrder,
        schema.limits.max_sample_cell_length,
      ),
    [sampleColumnOrder, sampleValues, schema.limits.max_sample_cell_length],
  );
  const samplesValid = rjsfValidator.isValid(
    schema.sampleSchema,
    sampleValues,
    schema.sampleSchema,
  );

  const setConfig = useCallback((value: unknown) => {
    dispatch({ type: 'config-form', value });
  }, []);
  const setOptions = useCallback((value: unknown) => {
    dispatch({ type: 'options-form', value });
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
    sampleTransportIssue,
    reviewPreviewAvailable:
      state.yamlIssue === null &&
      state.configFormIssue === null &&
      state.optionsFormIssue === null &&
      sampleTransportIssue === null &&
      review.ok,
    reviewReady:
      state.yamlIssue === null &&
      state.configFormIssue === null &&
      state.optionsFormIssue === null &&
      sampleTransportIssue === null &&
      configValid &&
      optionsValid &&
      samplesValid &&
      review.ok,
    setConfig,
    acceptConfigFormFallback() {
      dispatch({ type: 'config-form-fallback-accepted' });
    },
    setOptions,
    acceptOptionsFormFallback() {
      dispatch({ type: 'options-form-fallback-accepted' });
    },
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
