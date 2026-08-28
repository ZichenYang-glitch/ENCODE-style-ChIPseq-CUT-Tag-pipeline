import { Button } from '../../components/Button';
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
          Adapter options
        </h3>
        <p className="mt-1 text-xs text-[var(--color-text-muted)]">
          Fields and defaults are declared by this workflow adapter.
        </p>
      </div>
      <SchemaObjectForm
        schema={schema.optionSchema}
        value={draft.state.options}
        resetRevision={draft.state.optionsFormResetRevision}
        onChange={draft.setOptions}
        ariaLabel="Workflow options form"
      />
      {draft.state.optionsFormIssue && (
        <div
          className="flex flex-col gap-2 rounded border border-red-200 bg-[var(--color-error-bg)] px-3 py-2 text-sm text-[var(--color-error)] sm:flex-row sm:items-center sm:justify-between"
          role="alert"
        >
          <span>{draft.state.optionsFormIssue.message}</span>
          <Button
            type="button"
            variant="secondary"
            className="shrink-0"
            onClick={draft.acceptOptionsFormFallback}
          >
            Use previous safe options
          </Button>
        </div>
      )}
    </section>
  );
}
