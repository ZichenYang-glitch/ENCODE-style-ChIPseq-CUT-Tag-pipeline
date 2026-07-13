import type { InputDraftController } from './useInputDraft';
import type { WorkbenchSchema } from './schemaContract';
import { SchemaObjectForm } from './SchemaObjectForm';

interface OptionsEditorProps {
  schema: WorkbenchSchema;
  draft: InputDraftController;
}

export function OptionsEditor({ schema, draft }: OptionsEditorProps) {
  return (
    <section className="min-w-0 space-y-3" aria-labelledby="options-editor-title">
      <div>
        <h3 id="options-editor-title" className="text-sm font-semibold">
          Platform options
        </h3>
        <p className="mt-1 text-xs text-[var(--color-text-muted)]">
          Fields and defaults are declared by this workflow adapter.
        </p>
      </div>
      <SchemaObjectForm
        schema={schema.optionSchema}
        value={draft.state.options}
        onChange={draft.setOptions}
        ariaLabel="Workflow options form"
      />
    </section>
  );
}
