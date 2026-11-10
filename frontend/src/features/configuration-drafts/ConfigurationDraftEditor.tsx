import { useCallback, useEffect, useId, useRef, useState } from 'react';
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

const SCHEMA_RELOAD_REASON =
  'The editor contract changed. Discard these edits and reload before saving.';

export function ConfigurationDraftEditor({
  schema,
  schemaIdentity,
  discardEpoch,
  writesEnabled,
  onDirtyChange,
  onSave,
  onDiscard,
}: ConfigurationDraftEditorProps) {
  const [acceptedIdentity, setAcceptedIdentity] = useState(schemaIdentity);
  const [heldSchema, setHeldSchema] = useState(schema);
  const [reloadEpoch, setReloadEpoch] = useState(0);
  const [keptFor, setKeptFor] = useState<string | null>(null);
  const dirtyRef = useRef(false);
  const identityChanged = schemaIdentity !== acceptedIdentity;

  if (identityChanged && !dirtyRef.current) {
    setAcceptedIdentity(schemaIdentity);
    setHeldSchema(schema);
    setReloadEpoch((value) => value + 1);
  }

  const noteDirty = useCallback(
    (dirty: boolean) => {
      dirtyRef.current = dirty;
      onDirtyChange(dirty);
    },
    [onDirtyChange],
  );

  useEffect(() => {
    if (schemaIdentity !== acceptedIdentity && !dirtyRef.current) {
      setAcceptedIdentity(schemaIdentity);
      setHeldSchema(schema);
      setReloadEpoch((value) => value + 1);
    }
  }, [acceptedIdentity, schema, schemaIdentity, discardEpoch]);

  function reloadFromUpdatedSchema() {
    dirtyRef.current = false;
    setAcceptedIdentity(schemaIdentity);
    setHeldSchema(schema);
    setReloadEpoch((value) => value + 1);
  }

  function discardEdits() {
    if (schemaIdentity !== acceptedIdentity) reloadFromUpdatedSchema();
    onDiscard();
  }

  const activeSchema = identityChanged ? heldSchema : schema;

  return (
    <div className="space-y-4">
      {identityChanged && (
        <div
          role="status"
          className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface-subtle)] p-3"
        >
          <p>
            The editor contract changed. These edits stay here until you
            discard them. They are not checked against the updated contract.
          </p>
          <div className="flex flex-col gap-2 sm:flex-row">
            <Button
              type="button"
              variant="secondary"
              onClick={() => setKeptFor(schemaIdentity)}
            >
              Keep these edits
            </Button>
            <Button type="button" variant="danger" onClick={reloadFromUpdatedSchema}>
              Discard edits and reload
            </Button>
          </div>
          {keptFor === schemaIdentity && (
            <p>These edits remain on the previous contract.</p>
          )}
        </div>
      )}
      <DraftEditorSession
        key={`${acceptedIdentity}:${discardEpoch}:${reloadEpoch}`}
        schema={activeSchema}
        writesEnabled={writesEnabled}
        saveBlockedReason={identityChanged ? SCHEMA_RELOAD_REASON : null}
        onDirtyChange={noteDirty}
        onSave={onSave}
        onDiscard={discardEdits}
      />
    </div>
  );
}

function DraftEditorSession({
  schema,
  writesEnabled,
  saveBlockedReason,
  onDirtyChange,
  onSave,
  onDiscard,
}: {
  schema: WorkbenchSchema;
  writesEnabled: boolean;
  saveBlockedReason: string | null;
  onDirtyChange: (dirty: boolean) => void;
  onSave?: () => void;
  onDiscard: () => void;
}) {
  const reasonId = useId();
  const draft = useInputDraft(schema);
  const [mode, setMode] = useState<ConfigMode>('form');
  const [baseline, setBaseline] = useState<string | null>(null);
  const interacted = useRef(false);
  const current = stableSnapshot(draft.state.config, draft.state.options);
  const invalid =
    draft.state.yamlIssue !== null ||
    draft.state.configFormIssue !== null ||
    draft.state.optionsFormIssue !== null;
  const dirty =
    baseline !== null &&
    interacted.current &&
    (current !== baseline || invalid);
  const saveEnabled =
    dirty &&
    !invalid &&
    saveBlockedReason === null &&
    writesEnabled &&
    onSave !== undefined;

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

  let saveReason = 'Save is available after you change the draft.';
  if (saveBlockedReason) saveReason = saveBlockedReason;
  else if (!writesEnabled) saveReason = AUTHORING_WRITES_UNAVAILABLE;
  else if (dirty && invalid) {
    saveReason = 'Fix the editor errors before saving. These edits are kept.';
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
          aria-describedby={saveEnabled ? undefined : reasonId}
          onClick={() => onSave?.()}
        >
          Save edits
        </Button>
        <Button type="button" variant="secondary" disabled={!dirty} onClick={onDiscard}>
          Discard edits
        </Button>
      </div>
      {!saveEnabled && (
        <p id={reasonId} className="text-xs text-[var(--color-text-muted)]">
          {saveReason}
        </p>
      )}
    </div>
  );
}

function stableSnapshot(config: unknown, options: unknown): string {
  return JSON.stringify({ config, options });
}
