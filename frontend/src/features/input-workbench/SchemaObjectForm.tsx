import Form from '@rjsf/core';
import type { RJSFSchema } from '@rjsf/utils';
import type { ValidationRequestConfig } from '../../api/generated/models';
import { rjsfValidator } from './schemaContract';

interface SchemaObjectFormProps {
  schema: RJSFSchema;
  value: ValidationRequestConfig;
  resetRevision: number;
  onChange: (value: unknown) => void;
  ariaLabel: string;
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
        uiSchema={{
          'ui:submitButtonOptions': { norender: true },
        }}
      >
        <span className="hidden" aria-hidden="true" />
      </Form>
    </div>
  );
}
