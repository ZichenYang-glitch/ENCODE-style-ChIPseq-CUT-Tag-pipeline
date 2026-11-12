import { useCallback, useEffect, useRef, useState, type ReactNode } from 'react';
import { useBlocker, useLocation } from 'react-router-dom';
import { Button } from '../../components/Button';
import type { WorkbenchSchema } from '../input-workbench/schemaContract';
import type { InputDraftSeed } from '../input-workbench/useInputDraft';
import { AssistantConnectionPanel } from './AssistantConnectionPanel';
import {
  ConfigurationDraftEditor,
  type EditorDocument,
} from './ConfigurationDraftEditor';
import { ConfigurationDraftList } from './ConfigurationDraftList';
import { DraftValidationNotice } from './DraftValidationNotice';
import { ExplanationNotice } from './ExplanationNotice';
import { ProposalReview } from './ProposalReview';
import { StructuredRequestPanel } from './StructuredRequestPanel';
import { UnsavedEditsDialog } from './UnsavedEditsDialog';
import { AUTHORING_WRITES_UNAVAILABLE } from './model';
import type {
  AssistantConnectionView,
  AuthoringCapabilityView,
  DraftListModel,
  DraftValidationView,
  ExplanationView,
  ProposalView,
  StructuredQuestion,
  SuggestionRequestStatus,
} from './model';

export interface ConfigurationDraftWorkbenchProps {
  workflowName: string;
  projectName: string | null;
  schema: WorkbenchSchema | null;
  schemaIdentity: string;
  schemaMessage: string | null;
  list: DraftListModel;
  selectedSid: string | null;
  draftExpired: boolean;
  connection: AssistantConnectionView;
  capability: AuthoringCapabilityView | null;
  proposal: ProposalView | null;
  explanation: ExplanationView | null;
  validation: DraftValidationView;
  requestStatus: SuggestionRequestStatus;
  writesEnabled: boolean;
  apiConnected?: boolean;
  saveUnavailableReason?: string;
  documentSeed?: InputDraftSeed | null;
  externalDirty?: boolean;
  confirmedSnapshot?: string | null;
  serverReloadEpoch?: number;
  remoteConflict?: boolean;
  managedInputs?: ReactNode;
  adoptedReceipts?: ReactNode;
  actionNotice?: string | null;
  onSelectDraft?: (sid: string) => void;
  onListPrevious?: () => void;
  onListNext?: () => void;
  onRetryList?: () => void;
  onSave?: (document: EditorDocument) => void;
  onDiscarded?: () => void;
  onDirtyChange?: (dirty: boolean) => void;
  onReloadRemote?: () => void;
  onKeepRemote?: () => void;
  onConnect?: () => void;
  onClose?: () => void;
  onRequest?: (question: StructuredQuestion) => void;
  onCheckRequest?: () => void;
  onAdopt?: () => void;
  onReject?: () => void;
  onRefreshProposal?: () => void;
  onValidate?: () => void;
  operationPending?: boolean;
  externalPrompt?: boolean;
  onExternalDiscard?: () => void;
  onExternalCancel?: () => void;
}

export function ConfigurationDraftWorkbench({
  workflowName,
  projectName,
  schema,
  schemaIdentity,
  schemaMessage,
  list,
  selectedSid,
  draftExpired,
  connection,
  capability,
  proposal,
  explanation,
  validation,
  requestStatus,
  writesEnabled,
  apiConnected = false,
  saveUnavailableReason,
  documentSeed = null,
  externalDirty = false,
  confirmedSnapshot = null,
  serverReloadEpoch = 0,
  remoteConflict = false,
  managedInputs = null,
  adoptedReceipts = null,
  actionNotice = null,
  onSelectDraft,
  onListPrevious,
  onListNext,
  onRetryList,
  onSave,
  onDiscarded,
  onDirtyChange,
  onReloadRemote,
  onKeepRemote,
  onConnect,
  onClose,
  onRequest,
  onCheckRequest,
  onAdopt,
  onReject,
  onRefreshProposal,
  onValidate,
  operationPending = false,
  externalPrompt = false,
  onExternalDiscard,
  onExternalCancel,
}: ConfigurationDraftWorkbenchProps) {
  const [dirty, setDirty] = useState(false);
  const [discardEpoch, setDiscardEpoch] = useState(0);
  const [localNotice, setLocalNotice] = useState<string | null>(null);
  const [saveAttempted, setSaveAttempted] = useState(false);
  const [prompt, setPrompt] = useState<
    { kind: 'leave' } | { kind: 'select'; sid: string } | { kind: 'external' } | null
  >(null);
  const [editorCanSave, setEditorCanSave] = useState(false);
  const dirtyRef = useRef(false);
  const documentRef = useRef<EditorDocument | null>(null);
  const allowNavigation = useRef(false);
  const location = useLocation();
  dirtyRef.current = dirty;
  const blocker = useBlocker(
    useCallback(
      ({
        currentLocation,
        nextLocation,
      }: {
        currentLocation: { pathname: string };
        nextLocation: { pathname: string };
      }) =>
        dirtyRef.current &&
        !allowNavigation.current &&
        currentLocation.pathname !== nextLocation.pathname,
      [],
    ),
  );

  useEffect(() => {
    allowNavigation.current = false;
  }, [location.pathname, location.search]);

  useEffect(() => {
    if (dirty) setLocalNotice(null);
  }, [dirty]);

  useEffect(() => {
    if (blocker.state === 'blocked') {
      setPrompt((current) => current ?? { kind: 'leave' });
    }
  }, [blocker.state]);

  function selectDraft(sid: string) {
    if (sid === selectedSid) return;
    if (dirtyRef.current) {
      setSaveAttempted(false);
      setPrompt({ kind: 'select', sid });
      return;
    }
    onSelectDraft?.(sid);
  }

  const noteDirty = useCallback((next: boolean) => {
    setDirty(next);
    onDirtyChange?.(next);
  }, [onDirtyChange]);

  const noteDocument = useCallback((document: EditorDocument) => {
    documentRef.current = document;
  }, []);

  useEffect(() => {
    if (externalPrompt) {
      setPrompt((current) => current ?? { kind: 'external' });
      return;
    }
    setPrompt((current) => (current?.kind === 'external' ? null : current));
  }, [externalPrompt]);

  function isEditorDocument(value: unknown): value is EditorDocument {
    if (value == null || typeof value !== 'object') return false;
    if ('nativeEvent' in value || 'currentTarget' in value) return false;
    return 'config' in value && 'options' in value;
  }

  function requestSave(document?: EditorDocument) {
    const payload = isEditorDocument(document) ? document : documentRef.current;
    if (!editorCanSave || !onSave || !payload) return;
    onSave(payload);
    setSaveAttempted(true);
  }

  function discardLocal() {
    setDiscardEpoch((value) => value + 1);
    setSaveAttempted(false);
    setLocalNotice('Unsaved edits discarded.');
    onDiscarded?.();
  }

  function discardAndContinue() {
    const current = prompt;
    setDiscardEpoch((value) => value + 1);
    setSaveAttempted(false);
    setLocalNotice('Unsaved edits discarded.');
    onDiscarded?.();
    setPrompt(null);
    if (current?.kind === 'external') onExternalDiscard?.();
    if (current?.kind === 'select') onSelectDraft?.(current.sid);
    if (blocker.state === 'blocked') {
      allowNavigation.current = true;
      blocker.proceed();
    }
  }

  function cancelPrompt() {
    const current = prompt;
    setPrompt(null);
    setSaveAttempted(false);
    if (current?.kind === 'external') onExternalCancel?.();
    if (blocker.state === 'blocked') blocker.reset();
  }

  return (
    <section
      className="flex min-w-0 flex-1 flex-col gap-3"
      aria-labelledby="configuration-drafts-title"
    >
      <header className="space-y-1">
        <h1 id="configuration-drafts-title" className="text-lg font-semibold">
          Configuration drafts
        </h1>
        <p className="text-sm text-[var(--color-text-muted)]">
          Review your drafts, request a restricted suggestion, and validate by
          hand. This page does not start a run.
        </p>
        <p className="text-sm">
          <span className="text-[var(--color-text-muted)]">Workflow: </span>
          {workflowName}
        </p>
        <p className="text-sm">
          <span className="text-[var(--color-text-muted)]">Project: </span>
          {projectName ?? 'Project list is not connected yet.'}
        </p>
      </header>
      {!writesEnabled && !apiConnected && (
        <p role="status">{AUTHORING_WRITES_UNAVAILABLE}</p>
      )}
      {localNotice && <p role="status">{localNotice}</p>}
      {actionNotice && <p role="status">{actionNotice}</p>}
      {saveAttempted && dirty && (
        <p role="status">
          The save response has not been confirmed. These edits are still
          unsaved.
        </p>
      )}
      <div
        data-testid="configuration-draft-layout"
        className="flex min-w-0 flex-col gap-3 lg:flex-row"
      >
        <div className="flex min-w-0 flex-col gap-3 lg:w-80">
          <ConfigurationDraftList
            model={list}
            selectedSid={selectedSid}
            onSelect={selectDraft}
            onPrevious={onListPrevious}
            onNext={onListNext}
            onRetry={onRetryList}
          />
          <AssistantConnectionPanel
            connection={connection}
            writesEnabled={writesEnabled}
            unavailableMessage={saveUnavailableReason}
            onConnect={onConnect}
            onClose={onClose}
          />
        </div>
        <div className="flex min-w-0 flex-1 flex-col gap-3">
          {draftExpired ? (
            <p role="status">
              This draft has expired. Its configuration is hidden. Adopted
              changes from before it expired stay listed below.
            </p>
          ) : schema ? (
            <>
              <p className="text-sm text-[var(--color-text-muted)]">
                {selectedSid
                  ? 'Edits stay local until a confirmed save.'
                  : 'No server draft is open. Edits in this editor stay on this page.'}
              </p>
              {remoteConflict && (
                <div
                  role="status"
                  className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface-subtle)] p-3"
                >
                  <p>
                    The server has a newer draft. These edits stay here until
                    you choose. Nothing is overwritten.
                  </p>
                  <div className="flex flex-col gap-2 sm:flex-row">
                    <Button type="button" variant="secondary" onClick={() => onKeepRemote?.()}>
                      Keep these edits
                    </Button>
                    <Button type="button" variant="danger" onClick={() => onReloadRemote?.()}>
                      Reload server draft
                    </Button>
                  </div>
                </div>
              )}
              <ConfigurationDraftEditor
                schema={schema}
                schemaIdentity={schemaIdentity}
                discardEpoch={discardEpoch}
                serverReloadEpoch={serverReloadEpoch}
                writesEnabled={writesEnabled}
                saveUnavailableReason={saveUnavailableReason}
                documentSeed={documentSeed}
                externalDirty={externalDirty}
                confirmedSnapshot={confirmedSnapshot}
                managedInputs={managedInputs}
                onDirtyChange={noteDirty}
                onDocumentChange={noteDocument}
                onSaveEnabledChange={setEditorCanSave}
                onSave={onSave ? requestSave : undefined}
                onDiscard={discardLocal}
              />
            </>
          ) : (
            <p role="alert">
              {schemaMessage ?? 'The workflow editor is unavailable.'}
            </p>
          )}
          <StructuredRequestPanel
            capability={capability}
            connectionAvailable={connection.state === 'available'}
            draftExpired={draftExpired}
            writesEnabled={writesEnabled}
            unavailableMessage={saveUnavailableReason}
            dirty={dirty}
            requestStatus={requestStatus}
            operationPending={operationPending}
            onRequest={onRequest}
            onCheckRequest={onCheckRequest}
          />
          <ProposalReview
            proposal={proposal}
            writesEnabled={writesEnabled}
            unavailableMessage={saveUnavailableReason}
            dirty={dirty}
            onAdopt={onAdopt}
            onReject={onReject}
            onRefresh={onRefreshProposal}
            operationPending={operationPending}
          />
          <ExplanationNotice explanation={explanation} />
          <DraftValidationNotice
            validation={validation}
            writesEnabled={writesEnabled}
            unavailableMessage={saveUnavailableReason}
            dirty={dirty}
            draftExpired={draftExpired}
            onValidate={onValidate}
          />
          {adoptedReceipts}
        </div>
      </div>
      <UnsavedEditsDialog
        open={prompt !== null}
        writesEnabled={editorCanSave}
        unavailableMessage={
          writesEnabled
            ? 'Fix the editor errors before saving. These edits are kept.'
            : saveUnavailableReason
        }
        saveAttempted={saveAttempted}
        onSave={() => requestSave()}
        onDiscard={discardAndContinue}
        onCancel={cancelPrompt}
      />
    </section>
  );
}
