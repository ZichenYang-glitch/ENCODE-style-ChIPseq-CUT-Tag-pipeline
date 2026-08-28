import Form from '@rjsf/core';
import type { RJSFSchema, UiSchema } from '@rjsf/utils';
import type { ValidationRequestConfig } from '../../api/generated/models';
import { rjsfValidator } from './schemaContract';

interface SchemaObjectFormProps {
  schema: RJSFSchema;
  value: ValidationRequestConfig;
  resetRevision: number;
  onChange: (value: unknown) => void;
  ariaLabel: string;
}

function scalarFieldUiSchema(schema: RJSFSchema): UiSchema {
  const uiSchema: UiSchema = {
    'ui:submitButtonOptions': { norender: true },
  };
  for (const [key, candidate] of Object.entries(schema.properties ?? {})) {
    if (
      typeof candidate !== 'object' ||
      candidate === null ||
      Array.isArray(candidate) ||
      'oneOf' in candidate ||
      'anyOf' in candidate ||
      'allOf' in candidate ||
      '$ref' in candidate
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
      uiSchema[key] = { 'ui:classNames': 'hw-scalar-field' };
    }
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
        onChange={(event) => onChange(event.formData)}
        onSubmit={() => undefined}
        uiSchema={scalarFieldUiSchema(schema)}
      >
        <span className="hidden" aria-hidden="true" />
      </Form>
    </div>
  );
}
