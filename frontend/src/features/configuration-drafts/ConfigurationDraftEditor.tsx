import { useEffect, useRef, useState } from 'react';
import { Button } from '../../components/Button';
import { ConfigEditor, type ConfigMode } from '../input-workbench/ConfigEditor';
import { OptionsEditor } from '../input-workbench/OptionsEditor';
import type { WorkbenchSchema } from '../input-workbench/schemaContract';
import { useInputDraft } from '../input-workbench/useInputDraft';
import { AUTHORING_WRITES_UNAVAILABLE } from './model';

interface ConfigurationDraftEditorProps {
  schema: WorkbenchSchema;
  schemaIdentity: string;
  discardEpoch: number;
  writesEnabled: boolean;
  onDirtyChange: (dirty: boolean) => void;
  onSave?: () => void;
  onDiscard: () => void;
}

export function ConfigurationDraftEditor(props: ConfigurationDraftEditorProps) {
  return (
    <DraftEditorSession
      key={`${props.schemaIdentity}:${props.discardEpoch}`}
      {...props}
    />
  );
}

function DraftEditorSession({
  schema,
  writesEnabled,
  onDirtyChange,
  onSave,
  onDiscard,
}: ConfigurationDraftEditorProps) {
  const draft = useInputDraft(schema);
  const [mode, setMode] = useState<ConfigMode>('form');
  const [baseline, setBaseline] = useState<string | null>(null);
  const interacted = useRef(false);
  const current = stableSnapshot(draft.state.config, draft.state.options);
  const dirty =
    baseline !== null && interacted.current && current !== baseline;
  const saveEnabled = dirty && writesEnabled && onSave !== undefined;

  useEffect(() => {
    if (interacted.current) return;
    setBaseline((existing) => (existing === current ? existing : current));
  }, [current]);

  useEffect(() => {
    onDirtyChange(dirty);
    return () => onDirtyChange(false);
  }, [dirty, onDirtyChange]);

  function markInteraction() {
    interacted.current = true;
  }

  return (
    <div className="input-workbench space-y-4">
      <div
        className="space-y-4"
        onPointerDown={markInteraction}
        onKeyDown={markInteraction}
      >
        <ConfigEditor
          schema={schema}
          draft={draft}
          mode={mode}
          onModeChange={setMode}
        />
        <OptionsEditor schema={schema} draft={draft} />
      </div>
      <p className="text-xs text-[var(--color-text-muted)]">
        Sample selection stays in the existing input workbench. This page does
        not start a run.
      </p>
      {dirty && (
        <p role="status">
          You have unsaved edits. They are not saved, adopted, or validated.
        </p>
      )}
      <div className="flex flex-col gap-2 sm:flex-row">
        <Button
          type="button"
          variant="primary"
          disabled={!saveEnabled}
          aria-describedby={saveEnabled ? undefined : 'draft-save-reason'}
          onClick={() => onSave?.()}
        >
          Save edits
        </Button>
        <Button type="button" variant="secondary" disabled={!dirty} onClick={onDiscard}>
          Discard edits
        </Button>
      </div>
      {!saveEnabled && (
        <p id="draft-save-reason" className="text-xs text-[var(--color-text-muted)]">
          {!writesEnabled
            ? AUTHORING_WRITES_UNAVAILABLE
            : 'Save is available after you change the draft.'}
        </p>
      )}
    </div>
  );
}

function stableSnapshot(config: unknown, options: unknown): string {
  return JSON.stringify({ config, options });
}
