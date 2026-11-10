import { useCallback, useEffect, useRef, useState } from 'react';
import { useBlocker, useLocation } from 'react-router-dom';
import type { WorkbenchSchema } from '../input-workbench/schemaContract';
import { AssistantConnectionPanel } from './AssistantConnectionPanel';
import { ConfigurationDraftEditor } from './ConfigurationDraftEditor';
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
  onSelectDraft?: (sid: string) => void;
  onListPrevious?: () => void;
  onListNext?: () => void;
  onRetryList?: () => void;
  onSave?: () => void;
  onConnect?: () => void;
  onRequest?: (question: StructuredQuestion) => void;
  onAdopt?: () => void;
  onReject?: () => void;
  onRefreshProposal?: () => void;
  onValidate?: () => void;
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
  onSelectDraft,
  onListPrevious,
  onListNext,
  onRetryList,
  onSave,
  onConnect,
  onRequest,
  onAdopt,
  onReject,
  onRefreshProposal,
  onValidate,
}: ConfigurationDraftWorkbenchProps) {
  const [dirty, setDirty] = useState(false);
  const [discardEpoch, setDiscardEpoch] = useState(0);
  const [localNotice, setLocalNotice] = useState<string | null>(null);
  const [saveAttempted, setSaveAttempted] = useState(false);
  const [prompt, setPrompt] = useState<
    { kind: 'leave' } | { kind: 'select'; sid: string } | null
  >(null);
  const dirtyRef = useRef(false);
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
  const onDirtyChange = useCallback((next: boolean) => {
    setDirty((current) => (current === next ? current : next));
  }, []);

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

  function requestSave() {
    if (!writesEnabled || !onSave) return;
    onSave();
    setSaveAttempted(true);
  }

  function discardLocal() {
    setDiscardEpoch((value) => value + 1);
    setSaveAttempted(false);
    setLocalNotice('Unsaved edits discarded.');
  }

  function discardAndContinue() {
    const current = prompt;
    setDiscardEpoch((value) => value + 1);
    setSaveAttempted(false);
    setLocalNotice('Unsaved edits discarded.');
    setPrompt(null);
    if (current?.kind === 'select') onSelectDraft?.(current.sid);
    if (blocker.state === 'blocked') {
      allowNavigation.current = true;
      blocker.proceed();
    }
  }

  function cancelPrompt() {
    setPrompt(null);
    setSaveAttempted(false);
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
      {!writesEnabled && (
        <p role="status">{AUTHORING_WRITES_UNAVAILABLE}</p>
      )}
      {localNotice && <p role="status">{localNotice}</p>}
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
            onConnect={onConnect}
          />
        </div>
        <div className="flex min-w-0 flex-1 flex-col gap-3">
          {draftExpired ? (
            <p role="status">
              This draft has expired. Its configuration is hidden. Adoption
              receipts can be reviewed when that interface is connected.
            </p>
          ) : schema ? (
            <>
              <p className="text-sm text-[var(--color-text-muted)]">
                {selectedSid
                  ? 'Edits stay local until a confirmed save.'
                  : 'No server draft is open. Edits in this editor stay on this page.'}
              </p>
              <ConfigurationDraftEditor
                schema={schema}
                schemaIdentity={schemaIdentity}
                discardEpoch={discardEpoch}
                writesEnabled={writesEnabled}
                onDirtyChange={onDirtyChange}
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
            dirty={dirty}
            requestStatus={requestStatus}
            onRequest={onRequest}
          />
          <ProposalReview
            proposal={proposal}
            writesEnabled={writesEnabled}
            dirty={dirty}
            onAdopt={onAdopt}
            onReject={onReject}
            onRefresh={onRefreshProposal}
          />
          <ExplanationNotice explanation={explanation} />
          <DraftValidationNotice
            validation={validation}
            writesEnabled={writesEnabled}
            dirty={dirty}
            draftExpired={draftExpired}
            onValidate={onValidate}
          />
        </div>
      </div>
      <UnsavedEditsDialog
        open={prompt !== null}
        writesEnabled={writesEnabled && onSave !== undefined}
        saveAttempted={saveAttempted}
        onSave={requestSave}
        onDiscard={discardAndContinue}
        onCancel={cancelPrompt}
      />
    </section>
  );
}
