import Form from '@rjsf/core';
import type { RJSFSchema, UiSchema } from '@rjsf/utils';
import type { ValidationRequestConfig } from '../../api/generated/models';
import { rjsfValidator } from './schemaContract';
import { isJsonObject, isPlainObject } from './jsonSafety';
import { cascadeSectionSwitches, sectionDisableMode } from './sectionSwitches';

interface SchemaObjectFormProps {
  schema: RJSFSchema;
  value: ValidationRequestConfig;
  resetRevision: number;
  onChange: (value: unknown) => void;
  ariaLabel: string;
}

function fieldUiSchema(
  schema: RJSFSchema,
  value: unknown,
  styleScalarFields = true,
): UiSchema {
  const uiSchema: UiSchema = {
    'ui:submitButtonOptions': { norender: true },
  };
  const mode = sectionDisableMode(schema);
  for (const [key, candidate] of Object.entries(schema.properties ?? {})) {
    if (isPlainObject(candidate)) {
      const child = fieldUiSchema(
        candidate as RJSFSchema,
        isPlainObject(value) ? value[key] : undefined,
        false,
      );
      if (
        mode && key !== 'enabled' && isPlainObject(value) && value.enabled === false
      ) {
        if (mode === 'reset' || candidate.type === 'boolean') {
          child['ui:disabled'] = true;
        }
      }
      uiSchema[key] = child;
    }
    if (
      !styleScalarFields ||
      typeof candidate !== 'object' ||
      candidate === null ||
      Array.isArray(candidate) ||
      'oneOf' in candidate ||
      'anyOf' in candidate ||
      'allOf' in candidate ||
      '$ref' in candidate ||
      'enum' in candidate ||
      'format' in candidate ||
      'contentMediaType' in candidate ||
      'contentEncoding' in candidate
    ) {
      continue;
    }
    const type = candidate.type;
    if (
      type === 'string' ||
      type === 'number' ||
      type === 'integer' ||
      type === 'boolean'
    ) {
      uiSchema[key] = { ...uiSchema[key], 'ui:classNames': 'hw-scalar-field' };
    }
  }
  if (mode) {
    uiSchema.enabled = {
      ...uiSchema.enabled,
      'ui:help': mode === 'reset'
        ? 'Turning this off clears its settings. Turning it on again does not restore them.'
        : 'Turning this off clears its checkboxes. Select them again after turning it on.',
    };
  }
  return uiSchema;
}

export function SchemaObjectForm({
  schema,
  value,
  resetRevision,
  onChange,
  ariaLabel,
}: SchemaObjectFormProps) {
  return (
    <div className="schema-object-form" aria-label={ariaLabel}>
      <Form
        key={resetRevision}
        schema={schema}
        validator={rjsfValidator}
        formData={value}
        experimental_defaultFormStateBehavior={{
          emptyObjectFields: 'skipDefaults',
        }}
        liveOmit={false}
        liveValidate="onBlur"
        omitExtraData={false}
        noHtml5Validate
        showErrorList={false}
        onChange={(event) => onChange(
          isJsonObject(event.formData)
            ? cascadeSectionSwitches(schema, value, event.formData)
            : event.formData,
        )}
        onSubmit={() => undefined}
        uiSchema={fieldUiSchema(schema, value)}
      >
        <span className="hidden" aria-hidden="true" />
      </Form>
    </div>
  );
}
