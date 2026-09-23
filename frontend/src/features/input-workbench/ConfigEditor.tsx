import CodeMirror from '@uiw/react-codemirror';
import { useEffect, useRef } from 'react';
import { yaml } from '@codemirror/lang-yaml';
import { EditorView } from '@codemirror/view';
import { Code2, ListTree } from 'lucide-react';
import { Button } from '../../components/Button';
import type { InputDraftController } from './useInputDraft';
import type { WorkbenchSchema } from './schemaContract';
import { SchemaObjectForm } from './SchemaObjectForm';
import { focusConfigIssue, type ConfigFocusRequest, type SafeIssue } from './validationFeedback';

export type ConfigMode = 'form' | 'yaml';

interface ConfigEditorProps {
  schema: WorkbenchSchema;
  draft: InputDraftController;
  mode: ConfigMode;
  onModeChange: (mode: ConfigMode) => void;
  issues?: SafeIssue[];
  focusRequest?: ConfigFocusRequest | null;
}

export function ConfigEditor({
  schema,
  draft,
  mode,
  onModeChange,
  issues,
  focusRequest,
}: ConfigEditorProps) {
  const editor = useRef<HTMLElement>(null);
  useEffect(() => {
    if (!focusRequest || !editor.current) return;
    if (mode === 'yaml') {
      editor.current.querySelector<HTMLElement>('[aria-label="Advanced config YAML"]')?.focus();
    } else {
      focusConfigIssue(editor.current, schema.configSchema, focusRequest.path);
    }
  }, [focusRequest, mode, schema.configSchema]);
  return (
    <section ref={editor} tabIndex={-1} className="min-w-0 max-w-5xl space-y-3" aria-labelledby="config-editor-title">
      <div className="flex flex-col gap-2 sm:flex-row sm:items-start sm:justify-between">
        <div>
          <h3 id="config-editor-title" className="text-sm font-semibold">
            Workflow config
          </h3>
          <p className="mt-1 text-xs text-[var(--color-text-muted)]">
            {schema.contract.coverage.config === 'complete'
              ? 'Form and YAML use the same complete adapter-owned schema.'
              : 'Form coverage is partial. Additional adapter-owned keys remain available in YAML.'}
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
          resetRevision={draft.state.configFormResetRevision}
          onChange={draft.setConfig}
          ariaLabel="Workflow config form"
          issues={issues}
        />
      ) : (
        <div className="min-w-0 overflow-hidden rounded-[4px] border border-[var(--color-border)] bg-[var(--color-surface)]">
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
        {' '}To remove an optional config key, edit YAML. Turning a checkbox off keeps its value as false.
      </p>
      {draft.state.yamlIssue && (
        <p
          className="rounded-[4px] border border-[var(--color-error-border)] bg-[var(--color-error-bg)] px-3 py-2 text-sm text-[var(--color-error)]"
          role="alert"
        >
          {draft.state.yamlIssue.message}
        </p>
      )}
      {draft.state.configFormIssue && (
        <div
          data-testid="config-form-safety-issue"
          className="flex flex-col gap-2 rounded-[4px] border border-[var(--color-error-border)] bg-[var(--color-error-bg)] px-3 py-2 text-sm text-[var(--color-error)] sm:flex-row sm:items-center sm:justify-between"
          role="alert"
        >
          <span>{draft.state.configFormIssue.message}</span>
          <Button
            type="button"
            variant="secondary"
            className="shrink-0"
            onClick={draft.acceptConfigFormFallback}
          >
            Use previous safe config
          </Button>
        </div>
      )}
    </section>
  );
}
