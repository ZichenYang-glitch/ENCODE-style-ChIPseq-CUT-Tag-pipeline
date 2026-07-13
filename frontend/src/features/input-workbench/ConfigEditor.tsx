import CodeMirror from '@uiw/react-codemirror';
import { yaml } from '@codemirror/lang-yaml';
import { EditorView } from '@codemirror/view';
import { Code2, ListTree } from 'lucide-react';
import { Button } from '../../components/Button';
import type { InputDraftController } from './useInputDraft';
import type { WorkbenchSchema } from './schemaContract';
import { SchemaObjectForm } from './SchemaObjectForm';

export type ConfigMode = 'form' | 'yaml';

interface ConfigEditorProps {
  schema: WorkbenchSchema;
  draft: InputDraftController;
  mode: ConfigMode;
  onModeChange: (mode: ConfigMode) => void;
}

export function ConfigEditor({
  schema,
  draft,
  mode,
  onModeChange,
}: ConfigEditorProps) {
  return (
    <section className="min-w-0 space-y-3" aria-labelledby="config-editor-title">
      <div className="flex flex-col gap-2 sm:flex-row sm:items-start sm:justify-between">
        <div>
          <h3 id="config-editor-title" className="text-sm font-semibold">
            Workflow config
          </h3>
          <p className="mt-1 text-xs text-[var(--color-text-muted)]">
            Form coverage is {schema.contract.coverage.config}. Advanced keys
            remain available in YAML.
          </p>
        </div>
        <div
          className="inline-flex w-fit rounded border border-[var(--color-border)] p-0.5"
          aria-label="Config editor mode"
        >
          <Button
            type="button"
            variant={mode === 'form' ? 'primary' : 'secondary'}
            className="gap-1.5 border-0 px-2 py-1"
            aria-label="Form mode"
            aria-pressed={mode === 'form'}
            disabled={draft.state.yamlIssue !== null}
            onClick={() => onModeChange('form')}
          >
            <ListTree aria-hidden="true" size={15} />
            Form
          </Button>
          <Button
            type="button"
            variant={mode === 'yaml' ? 'primary' : 'secondary'}
            className="gap-1.5 border-0 px-2 py-1"
            aria-label="YAML mode"
            aria-pressed={mode === 'yaml'}
            onClick={() => onModeChange('yaml')}
          >
            <Code2 aria-hidden="true" size={15} />
            YAML
          </Button>
        </div>
      </div>

      {mode === 'form' ? (
        <SchemaObjectForm
          schema={schema.configSchema}
          value={draft.state.config}
          onChange={draft.setConfig}
          ariaLabel="Workflow config form"
        />
      ) : (
        <div className="min-w-0 overflow-hidden rounded border border-[var(--color-border)] bg-white">
          <CodeMirror
            value={draft.state.yamlText}
            height="28rem"
            basicSetup={{ lineNumbers: true, foldGutter: true }}
            extensions={[
              yaml(),
              EditorView.lineWrapping,
              EditorView.contentAttributes.of({
                'aria-label': 'Advanced config YAML',
                'aria-describedby': 'yaml-editor-help',
              }),
            ]}
            onChange={draft.editYaml}
          />
        </div>
      )}

      <p id="yaml-editor-help" className="text-xs text-[var(--color-text-muted)]">
        YAML is parsed locally into the same config object. Formatting and
        comments may be normalized after a Form edit.
      </p>
      {draft.state.yamlIssue && (
        <p
          className="rounded border border-red-200 bg-[var(--color-error-bg)] px-3 py-2 text-sm text-[var(--color-error)]"
          role="alert"
        >
          {draft.state.yamlIssue.message}
        </p>
      )}
    </section>
  );
}
