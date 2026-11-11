import { useCallback, useEffect, useRef, useState } from 'react';
import { useQuery, useQueryClient } from '@tanstack/react-query';
import {
  authoringAdoptProposal,
  authoringCloseConnection,
  authoringCreateSession,
  authoringCreateTurn,
  authoringGetSession,
  authoringGetTurn,
  authoringListProjectSamples,
  authoringListProjects,
  authoringListSessions,
  authoringOpenConnection,
  authoringRejectProposal,
  authoringSaveDraft,
  authoringValidateDraft,
} from '../../api/generated/authoring/authoring';
import { ApiError } from '../../api/fetcher';
import type {
  AuthoringAdoptRequest,
  AuthoringDraftDocument,
  AuthoringRejectRequest,
  AuthoringSessionResponse,
  AuthoringTurnCreateRequest,
  AuthoringTurnResponse,
  ProjectSampleSummaryResponse,
  ProjectSummaryResponse,
} from '../../api/generated/models';
import {
  getWorkflowSchema,
  listCompatibleReferenceProfiles,
} from '../../api/generated/workflows/workflows';
import { useClients } from '../../api/client-context';
import { readReferenceProfileSummary } from '../../api/runTypes';
import { useAuth } from '../../app/auth';
import { Button } from '../../components/Button';
import {
  capabilityView,
  classifyAuthoringError,
  connectionView,
  draftSaveDisposition,
  draftWriteRequest,
  type CanonicalDraftContent,
  explanationView,
  isValidationBusinessFailure,
  listLabel,
  proposalView,
  questionRequest,
  sameMembers,
  stableDocumentSnapshot,
  turnCreateBody,
  turnResponseIsCurrent,
  turnScope,
  validationFromIssues,
  validationFromResponse,
  validationView,
  withMissingWorkflowSampleContext,
  workflowSampleRowsAbsent,
  type StoredValidation,
  type TurnScope,
} from '../../features/configuration-drafts/authoringContract';
import { ConfigurationDraftWorkbench } from '../../features/configuration-drafts/ConfigurationDraftWorkbench';
import type { EditorDocument } from '../../features/configuration-drafts/ConfigurationDraftEditor';
import { ManagedDraftInputs } from '../../features/configuration-drafts/ManagedDraftInputs';
import type {
  DraftValidationView,
  ProposalView,
  StructuredQuestion,
  SuggestionRequestStatus,
} from '../../features/configuration-drafts/model';
import { readWorkbenchSchema } from '../../features/input-workbench/schemaContract';
import type { InputDraftSeed } from '../../features/input-workbench/useInputDraft';

const PAGE_LIMIT = 20;

interface ConfigurationDraftsPageProps {
  workflowId: string;
}

type DecisionName = 'adopted' | 'rejected' | 'unconfirmed';

interface PreservedDraft {
  samples: AuthoringDraftDocument['samples'];
  selections: AuthoringDraftDocument['input_selections'];
}

export function ConfigurationDraftsPage({ workflowId }: ConfigurationDraftsPageProps) {
  const { loading: authLoading, principal } = useAuth();
  const userId = principal?.user_id ?? null;
  const queryClient = useQueryClient();
  const { workflowClient } = useClients();
  const scope = useRef(0);
  const preservedRef = useRef<PreservedDraft>({ samples: null, selections: [] });
  const pendingTurn = useRef<AuthoringTurnCreateRequest | null>(null);
  const pendingDecision = useRef<
    | { kind: 'adopt'; proposalId: string; body: AuthoringAdoptRequest }
    | { kind: 'reject'; proposalId: string; body: AuthoringRejectRequest }
    | null
  >(null);
  const busy = useRef(false);
  const operationGeneration = useRef(0);
  const selectedSidRef = useRef<string | null>(null);
  const pendingRemote = useRef<AuthoringSessionResponse | null>(null);

  const [projectId, setProjectId] = useState<string | null>(null);
  const [extraProjects, setExtraProjects] = useState<ProjectSummaryResponse[]>([]);
  const [cursors, setCursors] = useState<Array<string | null>>([null]);
  const [selectedSid, setSelectedSid] = useState<string | null>(null);
  const [seed, setSeed] = useState<InputDraftSeed | null>(null);
  const [seedRevision, setSeedRevision] = useState<number | null>(null);
  const [ignoredRevision, setIgnoredRevision] = useState<number | null>(null);
  const [reloadEpoch, setReloadEpoch] = useState(0);
  const [confirmedSnapshot, setConfirmedSnapshot] = useState<string | null>(null);
  const [sampleIds, setSampleIds] = useState<string[]>([]);
  const [extraSamples, setExtraSamples] = useState<ProjectSampleSummaryResponse[]>([]);
  const [referenceId, setReferenceId] = useState<string | null>(null);
  const [dirty, setDirty] = useState(false);
  const [remoteConflict, setRemoteConflict] = useState(false);
  const [notice, setNotice] = useState<string | null>(null);
  const [connectPending, setConnectPending] = useState(false);
  const [bridgeBlocked, setBridgeBlocked] = useState(false);
  const [requestStatus, setRequestStatus] = useState<SuggestionRequestStatus>('idle');
  const [turn, setTurn] = useState<AuthoringTurnResponse | null>(null);
  const [decision, setDecision] = useState<DecisionName | null>(null);
  const [storedValidation, setStoredValidation] = useState<StoredValidation | null>(null);
  const [validationOverride, setValidationOverride] = useState<DraftValidationView | null>(null);
  const [creating, setCreating] = useState(false);
  const [operationPending, setOperationPending] = useState(false);
  const [pendingProjectId, setPendingProjectId] = useState<string | null>(null);
  const [pendingCreate, setPendingCreate] = useState(false);
  selectedSidRef.current = selectedSid;

  const cursor = cursors[cursors.length - 1] ?? null;
  const stillCurrent = useCallback((captured: number, sid: string) => {
    return captured === scope.current && sid === selectedSidRef.current;
  }, []);

  useEffect(() => {
    abandonOperation();
    pendingRemote.current = null;
    scope.current += 1;
    setProjectId(null);
    setExtraProjects([]);
    setCursors([null]);
    setSelectedSid(null);
    setSeed(null);
    setSeedRevision(null);
    setIgnoredRevision(null);
    setSampleIds([]);
    setExtraSamples([]);
    setReferenceId(null);
    setRemoteConflict(false);
    setNotice(null);
    setBridgeBlocked(false);
    setConnectPending(false);
    setRequestStatus('idle');
    setTurn(null);
    setDecision(null);
    setStoredValidation(null);
    setValidationOverride(null);
    pendingTurn.current = null;
    pendingDecision.current = null;
  }, [userId]);

  const detailQuery = useQuery({
    queryKey: ['workflow', workflowId],
    queryFn: () => workflowClient.getWorkflow(workflowId),
    retry: false,
  });
  const schemaQuery = useQuery({
    queryKey: ['workflow-authoring-schema', workflowId],
    queryFn: () => getWorkflowSchema(workflowId),
    retry: false,
  });
  const projectsQuery = useQuery({
    queryKey: ['authoring-projects', userId],
    queryFn: () => authoringListProjects({ limit: 50 }),
    enabled: userId !== null,
    retry: false,
  });
  const sessionsQuery = useQuery({
    queryKey: ['authoring-sessions', userId, projectId, workflowId, cursor],
    queryFn: () =>
      authoringListSessions({
        project_id: projectId,
        workflow_id: workflowId,
        ...(typeof cursor === 'string' ? { cursor } : {}),
        limit: PAGE_LIMIT,
      }),
    enabled: userId !== null && projectId !== null,
    retry: false,
  });
  const sessionQuery = useQuery({
    queryKey: ['authoring-session', userId, selectedSid],
    queryFn: () => authoringGetSession(selectedSid ?? ''),
    enabled: userId !== null && selectedSid !== null,
    retry: false,
  });
  const samplesQuery = useQuery({
    queryKey: ['authoring-samples', userId, projectId],
    queryFn: () => authoringListProjectSamples(projectId ?? '', { limit: 50 }),
    enabled: userId !== null && projectId !== null,
    retry: false,
  });
  const referenceQuery = useQuery({
    queryKey: ['workflow-reference-profiles', workflowId],
    queryFn: () => listCompatibleReferenceProfiles(workflowId),
    retry: false,
  });

  const firstProjectId = projectsQuery.data?.projects[0]?.project_id ?? null;
  useEffect(() => {
    if (projectId || !firstProjectId) return;
    setProjectId(firstProjectId);
  }, [firstProjectId, projectId]);

  const projects = [...(projectsQuery.data?.projects ?? []), ...extraProjects];
  const session = sessionQuery.data ?? null;
  useEffect(() => {
    if (!session || session.sid !== selectedSid || !session.draft || session.draft_state !== 'active') {
      return;
    }
    if (session.draft_revision === seedRevision) return;
    if (dirty) {
      if (session.draft_revision !== ignoredRevision) setRemoteConflict(true);
      return;
    }
    preservedRef.current = {
      samples: session.draft.samples ?? null,
      selections: session.draft.input_selections ?? [],
    };
    setSampleIds([...(session.draft.sample_revision_ids ?? [])]);
    setReferenceId(session.draft.reference_profile_revision_id ?? null);
    setSeed({ config: session.draft.config, options: session.draft.options });
    setSeedRevision(session.draft_revision);
    setIgnoredRevision(null);
    setRemoteConflict(false);
    setConfirmedSnapshot(null);
    setReloadEpoch((value) => value + 1);
  }, [dirty, ignoredRevision, seedRevision, selectedSid, session]);

  const workflowName =
    detailQuery.data?.ok === true
      ? detailQuery.data.workflow?.metadata.name ?? 'Workflow name unavailable'
      : detailQuery.isError || detailQuery.data?.ok === false
        ? 'Workflow name unavailable'
        : 'Loading workflow';
  const parsed = schemaQuery.data?.schema ? readWorkbenchSchema(schemaQuery.data.schema) : null;
  const schema = parsed?.ok ? parsed.value : null;
  const schemaMessage = schemaQuery.isLoading
    ? null
    : schemaQuery.isError
      ? 'The workflow editor could not be loaded.'
      : schemaQuery.data && !schemaQuery.data.schema
        ? 'This workflow does not publish an editor contract.'
        : parsed && !parsed.ok
          ? parsed.message
          : null;
  const projectName = projects.find((project) => project.project_id === projectId)?.display_name ?? null;
  const activeDraft = session?.draft_state === 'active' ? session.draft : null;
  const inputsDirty = Boolean(
    activeDraft &&
      (!sameMembers(sampleIds, activeDraft.sample_revision_ids ?? []) ||
        referenceId !== (activeDraft.reference_profile_revision_id ?? null)),
  );
  const missing =
    sessionQuery.isError && classifyAuthoringError(sessionQuery.error).kind === 'not-found';
  const expired = session?.draft_state === 'expired' || session?.draft === null;
  const draftOpen = Boolean(session && session.draft_state === 'active' && session.draft && !missing);
  const writesEnabled = Boolean(userId) && draftOpen;
  const capability = capabilityView(session?.authoring_capability);
  const samples = [...(samplesQuery.data?.samples ?? []), ...extraSamples];
  const profiles = referenceQuery.data
    ? referenceQuery.data.profiles.flatMap((profile) => {
        const summary = readReferenceProfileSummary(profile);
        return summary ? [summary] : [];
      })
    : null;
  const listedSampleIds = new Set(samples.map((sample) => sample.sample_revision_id));
  const unlistedCount = sampleIds.filter((id) => !listedSampleIds.has(id)).length;
  const referenceMissing = Boolean(
    referenceId && profiles && !profiles.some((profile) => profile.revision_id === referenceId),
  );

  function sessionKey(sid: string) {
    return ['authoring-session', userId, sid] as const;
  }

  function claimOperation(): number | null {
    if (busy.current) return null;
    busy.current = true;
    setOperationPending(true);
    return operationGeneration.current;
  }

  function releaseOperation(generation: number) {
    if (generation !== operationGeneration.current) return;
    busy.current = false;
    setOperationPending(false);
  }

  function abandonOperation() {
    operationGeneration.current += 1;
    busy.current = false;
    setOperationPending(false);
    setConnectPending(false);
    setCreating(false);
  }

  function confirmSavedDraft(
    saved: AuthoringSessionResponse,
    sent: string,
    sentSamples: readonly string[],
    sentReference: string | null,
  ) {
    if (!saved.draft) return;
    pendingRemote.current = null;
    setConfirmedSnapshot(sent);
    setSeedRevision(saved.draft_revision);
    setSeed({ config: saved.draft.config, options: saved.draft.options });
    setSampleIds((current) =>
      sameMembers(current, sentSamples)
        ? [...(saved.draft?.sample_revision_ids ?? [])]
        : current,
    );
    setReferenceId((current) =>
      current === sentReference
        ? (saved.draft?.reference_profile_revision_id ?? null)
        : current,
    );
    preservedRef.current = {
      samples: saved.draft.samples ?? null,
      selections: saved.draft.input_selections ?? [],
    };
    setRemoteConflict(false);
    setNotice(null);
    rememberSession(saved);
  }

  function rememberSession(next: AuthoringSessionResponse) {
    preservedRef.current = {
      samples: next.draft?.samples ?? null,
      selections: next.draft?.input_selections ?? [],
    };
    queryClient.setQueryData(sessionKey(next.sid), next);
  }

  function applyTurn(result: AuthoringTurnResponse, requested: TurnScope) {
    const current = queryClient.getQueryData<AuthoringSessionResponse>(sessionKey(requested.sid));
    if (!current || !turnResponseIsCurrent(requested, turnScope(current, requested.userId))) {
      setNotice('This response is not for the current draft.');
      setRequestStatus('idle');
      return;
    }
    setTurn(result);
    setDecision(null);
    if (result.status === 'generating') setRequestStatus('unconfirmed');
    else {
      setRequestStatus(result.status === 'failed' ? 'failed' : 'idle');
      pendingTurn.current = null;
    }
  }

  async function createDraft() {
    if (!userId || !projectId) return;
    const generation = claimOperation();
    if (generation === null) return;
    const captured = scope.current;
    const principalId = userId;
    setCreating(true);
    try {
      const created = await authoringCreateSession({
        project_id: projectId,
        workflow_id: workflowId,
      });
      if (captured !== scope.current || principalId !== principal?.user_id) return;
      pendingRemote.current = null;
      rememberSession(created);
      setSelectedSid(created.sid);
      setSeedRevision(null);
      setTurn(null);
      setDecision(null);
      setStoredValidation(null);
      setValidationOverride(null);
      setNotice(null);
      await queryClient.invalidateQueries({
        queryKey: ['authoring-sessions', principalId, projectId, workflowId],
      });
    } catch (error) {
      if (captured === scope.current) setNotice(classifyAuthoringError(error).message);
    } finally {
      if (generation === operationGeneration.current) setCreating(false);
      releaseOperation(generation);
    }
  }

  async function saveDraft(document: EditorDocument) {
    if (!draftOpen || !session || seedRevision === null || !selectedSid || busy.current) return;
    const body = draftWriteRequest({
      expectedRevision: seedRevision,
      config: document.config,
      options: document.options,
      sampleRevisionIds: sampleIds,
      referenceProfileRevisionId: referenceId,
      preservedSamples: preservedRef.current.samples,
      preservedInputSelections: preservedRef.current.selections,
    });
    if (!body) {
      setNotice('These edits cannot be sent. They are still on this page.');
      return;
    }
    const captured = scope.current;
    const sid = selectedSid;
    const sentDocument = stableDocumentSnapshot(document.config, document.options);
    const sentSamples = [...sampleIds];
    const sentReference = referenceId;
    const previousRevision = seedRevision;
    const generation = claimOperation();
    if (generation === null) return;
    try {
      const saved = await authoringSaveDraft(sid, body);
      if (!stillCurrent(captured, sid)) return;
      const disposition = draftSaveDisposition({
        sent: body,
        observed: saved.draft,
        previousRevision,
        observedRevision: saved.draft_revision,
        acknowledged: true,
      });
      if (disposition === 'confirmed') {
        confirmSavedDraft(saved, sentDocument, sentSamples, sentReference);
        return;
      }
      if (disposition === 'conflict') noteUnmatchedSave(saved);
      else setNotice('The save response was not confirmed. These edits are still unsaved.');
    } catch (error) {
      if (!stillCurrent(captured, sid)) return;
      const failure = classifyAuthoringError(error);
      if (failure.kind === 'network') {
        await recoverSave(
          captured,
          sid,
          previousRevision,
          body,
          sentDocument,
          sentSamples,
          sentReference,
        );
        return;
      }
      if (failure.kind === 'stale') {
        setRemoteConflict(true);
        setNotice(failure.message);
        await queryClient.invalidateQueries({ queryKey: sessionKey(sid) });
        return;
      }
      setNotice(failure.message);
    } finally {
      releaseOperation(generation);
    }
  }

  function noteUnmatchedSave(observed: AuthoringSessionResponse) {
    pendingRemote.current = observed;
    setRemoteConflict(true);
    setNotice('The server has a newer draft. These edits are still here.');
  }

  async function recoverSave(
    captured: number,
    sid: string,
    previousRevision: number,
    sent: CanonicalDraftContent,
    sentDocument: string,
    sentSamples: readonly string[],
    sentReference: string | null,
  ) {
    try {
      const recovered = await authoringGetSession(sid);
      if (!stillCurrent(captured, sid)) return;
      const disposition = draftSaveDisposition({
        sent,
        observed: recovered.draft,
        previousRevision,
        observedRevision: recovered.draft_revision,
        acknowledged: false,
      });
      if (disposition === 'confirmed') {
        confirmSavedDraft(recovered, sentDocument, sentSamples, sentReference);
        return;
      }
      if (disposition === 'conflict') {
        noteUnmatchedSave(recovered);
        return;
      }
      setNotice('The save response was not confirmed. These edits are still unsaved.');
    } catch {
      if (stillCurrent(captured, sid)) {
        setNotice('The save response was not confirmed. These edits are still unsaved.');
      }
    }
  }

  async function connect() {
    if (!draftOpen || !selectedSid || connectPending) return;
    const generation = claimOperation();
    if (generation === null) return;
    const captured = scope.current;
    const sid = selectedSid;
    setConnectPending(true);
    try {
      const opened = await authoringOpenConnection(sid);
      if (!stillCurrent(captured, sid)) return;
      setBridgeBlocked(false);
      setNotice(null);
      rememberSession(opened);
    } catch (error) {
      if (!stillCurrent(captured, sid)) return;
      const failure = classifyAuthoringError(error);
      if (failure.kind === 'bridge-unavailable') setBridgeBlocked(true);
      setNotice(failure.message);
      if (failure.kind === 'network' || failure.kind === 'bridge-unavailable') {
        await queryClient.invalidateQueries({ queryKey: sessionKey(sid) });
      }
    } finally {
      if (generation === operationGeneration.current) setConnectPending(false);
      releaseOperation(generation);
    }
  }

  async function closeConnection() {
    if (!selectedSid) return;
    const generation = claimOperation();
    if (generation === null) return;
    const captured = scope.current;
    const sid = selectedSid;
    try {
      const closed = await authoringCloseConnection(sid);
      if (!stillCurrent(captured, sid)) return;
      setBridgeBlocked(false);
      setNotice(null);
      rememberSession(closed);
    } catch (error) {
      if (stillCurrent(captured, sid)) setNotice(classifyAuthoringError(error).message);
    } finally {
      releaseOperation(generation);
    }
  }

  async function requestSuggestion(question: StructuredQuestion) {
    if (requestStatus === 'unconfirmed') return;
    const generation = claimOperation();
    if (generation === null) return;
    if (!draftOpen || !session || !capability || seedRevision === null || !selectedSid) {
      releaseOperation(generation);
      return;
    }
    const built = questionRequest(question, capability);
    if (!built.ok) {
      setNotice(built.message);
      releaseOperation(generation);
      return;
    }
    const notes = [
      built.value.omittedNonIntegerPreferences
        ? 'Preferences that are not whole numbers were left out of this request.'
        : null,
    ].filter((item): item is string => item !== null);
    setNotice(notes.length > 0 ? notes.join(' ') : null);
    const body = turnCreateBody(newClientTurnId(), seedRevision, built.value.question);
    pendingTurn.current = body;
    try {
      await sendTurn(body, turnScope(session, userId ?? ''));
    } finally {
      releaseOperation(generation);
    }
  }

  async function sendTurn(body: AuthoringTurnCreateRequest, requested: TurnScope) {
    const captured = scope.current;
    try {
      const result = await authoringCreateTurn(requested.sid, body);
      if (captured !== scope.current) return;
      applyTurn(result, requested);
    } catch (error) {
      if (captured !== scope.current) return;
      const failure = classifyAuthoringError(error);
      if (failure.kind === 'network') {
        setRequestStatus('unconfirmed');
        setNotice(failure.message);
        return;
      }
      pendingTurn.current = null;
      setRequestStatus('failed');
      setNotice(failure.message);
      if (failure.kind === 'capability-changed') {
        setTurn(null);
        await queryClient.invalidateQueries({ queryKey: sessionKey(requested.sid) });
      }
    }
  }

  async function checkRequest() {
    if (!session || !selectedSid) return;
    const generation = claimOperation();
    if (generation === null) return;
    try {
      if (turn?.turn_id && (turn.status === 'generating' || requestStatus === 'unconfirmed') && !pendingTurn.current) {
        const captured = scope.current;
        const requested = turnScope(session, userId ?? '');
        try {
          const fresh = await authoringGetTurn(selectedSid, turn.turn_id);
          if (captured !== scope.current) return;
          applyTurn(fresh, requested);
        } catch (error) {
          if (captured === scope.current) setNotice(classifyAuthoringError(error).message);
        }
        return;
      }
      if (!pendingTurn.current) return;
      await sendTurn(pendingTurn.current, turnScope(session, userId ?? ''));
    } finally {
      releaseOperation(generation);
    }
  }

  async function adopt() {
    const generation = claimOperation();
    if (generation === null) return;
    if (!turn?.proposal_id || typeof turn.proposal_digest !== 'string' || seedRevision === null || !selectedSid) {
      releaseOperation(generation);
      return;
    }
    const body: AuthoringAdoptRequest = {
      expected_revision: seedRevision,
      proposal_digest: turn.proposal_digest,
    };
    pendingDecision.current = { kind: 'adopt', proposalId: turn.proposal_id, body };
    try {
      await sendDecision();
    } finally {
      releaseOperation(generation);
    }
  }

  async function reject() {
    const generation = claimOperation();
    if (generation === null) return;
    if (!turn?.proposal_id || typeof turn.proposal_digest !== 'string' || !selectedSid) {
      releaseOperation(generation);
      return;
    }
    const body: AuthoringRejectRequest = { proposal_digest: turn.proposal_digest };
    pendingDecision.current = { kind: 'reject', proposalId: turn.proposal_id, body };
    try {
      await sendDecision();
    } finally {
      releaseOperation(generation);
    }
  }

  async function sendDecision() {
    const pending = pendingDecision.current;
    if (!pending || !selectedSid || !session) return;
    const captured = scope.current;
    const sid = selectedSid;
    const requested = turnScope(session, userId ?? '');
    try {
      if (pending.kind === 'adopt') {
        await authoringAdoptProposal(sid, pending.proposalId, pending.body);
        if (!stillCurrent(captured, sid)) return;
        setDecision('adopted');
        setNotice(null);
        try {
          if (turn?.turn_id) {
            const fresh = await authoringGetTurn(sid, turn.turn_id);
            if (stillCurrent(captured, sid)) applyTurn(fresh, requested);
          }
          const reloaded = await authoringGetSession(sid);
          if (!stillCurrent(captured, sid)) return;
          if (dirty) {
            setRemoteConflict(true);
            setNotice(
              'Adoption was confirmed. These edits are still here, so the server draft was not loaded over them.',
            );
          } else {
            setSeedRevision(null);
            rememberSession(reloaded);
          }
        } catch {
          if (stillCurrent(captured, sid)) {
            setNotice('Adoption was confirmed. The updated draft has not been loaded.');
          }
        }
      } else {
        const fresh = await authoringRejectProposal(sid, pending.proposalId, pending.body);
        if (!stillCurrent(captured, sid)) return;
        applyTurn(fresh, requested);
        setNotice(null);
      }
      pendingDecision.current = null;
    } catch (error) {
      if (!stillCurrent(captured, sid)) return;
      const failure = classifyAuthoringError(error);
      if (failure.kind === 'network') {
        setDecision('unconfirmed');
        setNotice(failure.message);
        return;
      }
      if (failure.kind === 'expired' || failure.reasonCode === 'PROPOSAL_INVALIDATED') {
        setDecision(null);
        setTurn((current) =>
          current
            ? {
                ...current,
                adoptable: false,
                proposal_state: failure.reasonCode === 'PROPOSAL_EXPIRED' ? current.proposal_state : 'invalidated',
                reason_code: failure.reasonCode,
              }
            : current,
        );
      }
      if (failure.kind === 'stale') setRemoteConflict(true);
      pendingDecision.current = null;
      setNotice(failure.message);
    }
  }

  async function refreshProposal() {
    const generation = claimOperation();
    if (generation === null) return;
    try {
      if (pendingDecision.current) {
        await sendDecision();
        return;
      }
      if (!turn?.turn_id || !session || !selectedSid) return;
      const captured = scope.current;
      const requested = turnScope(session, userId ?? '');
      try {
        const fresh = await authoringGetTurn(selectedSid, turn.turn_id);
        if (captured !== scope.current) return;
        applyTurn(fresh, requested);
      } catch (error) {
        if (captured === scope.current) setNotice(classifyAuthoringError(error).message);
      }
    } finally {
      releaseOperation(generation);
    }
  }

  function validationWithSampleContext(
    stored: StoredValidation,
    issues: readonly { code: string; message: string }[],
    samples: unknown,
  ): StoredValidation {
    const summary = withMissingWorkflowSampleContext(stored.view.summary ?? '', issues, samples);
    if (summary === stored.view.summary) return stored;
    return { ...stored, view: { ...stored.view, summary } };
  }

  async function validateDraft() {
    if (!draftOpen || !session || seedRevision === null || !selectedSid) return;
    const generation = claimOperation();
    if (generation === null) return;
    const captured = scope.current;
    const sid = selectedSid;
    const revision = seedRevision;
    const schemaDigest = session.schema_digest;
    const samples = session.draft?.samples;
    try {
      const result = await authoringValidateDraft(sid, { expected_revision: revision });
      if (!stillCurrent(captured, sid)) return;
      setStoredValidation(
        validationWithSampleContext(
          validationFromResponse(result, schemaDigest),
          result.issues ?? [],
          samples,
        ),
      );
      setValidationOverride(null);
      setNotice(null);
    } catch (error) {
      if (!stillCurrent(captured, sid)) return;
      const failure = classifyAuthoringError(error);
      if (failure.kind === 'network') {
        setValidationOverride({ status: 'unconfirmed' });
        setNotice(failure.message);
        return;
      }
      if (isValidationBusinessFailure(failure) && error instanceof ApiError) {
        setStoredValidation(
          validationWithSampleContext(
            validationFromIssues(revision, schemaDigest, error.issues),
            error.issues,
            samples,
          ),
        );
        setValidationOverride(null);
        setNotice(null);
        return;
      }
      setValidationOverride({ status: 'request-failed', summary: failure.message });
      setNotice(failure.message);
    } finally {
      releaseOperation(generation);
    }
  }

  function reloadRemote() {
    const pending = pendingRemote.current;
    pendingRemote.current = null;
    if (pending && pending.sid === selectedSidRef.current) rememberSession(pending);
    const current = pending && pending.sid === selectedSidRef.current ? pending : session;
    if (!current?.draft) return;
    preservedRef.current = {
      samples: current.draft.samples ?? null,
      selections: current.draft.input_selections ?? [],
    };
    setSampleIds([...(current.draft.sample_revision_ids ?? [])]);
    setReferenceId(current.draft.reference_profile_revision_id ?? null);
    setSeed({ config: current.draft.config, options: current.draft.options });
    setSeedRevision(current.draft_revision);
    setIgnoredRevision(null);
    setRemoteConflict(false);
    setConfirmedSnapshot(null);
    setReloadEpoch((value) => value + 1);
  }

  if (authLoading || schemaQuery.isLoading) {
    return (
      <section
        className="min-h-64 min-w-0 flex-1 animate-pulse rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
        aria-label="Loading configuration drafts"
      >
        <div className="h-5 w-48 rounded-[4px] bg-[var(--color-skeleton)]" />
        <div className="mt-4 h-10 rounded-[4px] bg-[var(--color-skeleton)]" />
        <div className="mt-4 h-64 rounded-[4px] bg-[var(--color-skeleton)]" />
      </section>
    );
  }

  if (!userId) {
    return (
      <section className="flex min-w-0 flex-1 flex-col gap-3" aria-labelledby="configuration-drafts-title">
        <h1 id="configuration-drafts-title" className="text-lg font-semibold">
          Configuration drafts
        </h1>
        <p role="status">Sign in to view your configuration drafts.</p>
      </section>
    );
  }

  const draftLoading = Boolean(selectedSid) && !sessionQuery.isError && (
    sessionQuery.isLoading ||
    (session?.draft_state === 'active' &&
      seedRevision !== session.draft_revision &&
      !dirty &&
      !remoteConflict)
  );
  const serverConnection = session?.connection ?? { state: 'closed' as const };
  let connection = connectionView(serverConnection, {
    draftRestored: Boolean(session?.requires_new_connection),
    connectPending,
  });
  if (bridgeBlocked && connection.state !== 'available' && !connectPending) {
    connection = { ...connection, state: 'bridge-unavailable' };
  }
  const proposal = withDecision(turn ? proposalView(turn, capability) : null, decision);
  const boundValidation = validationView({
    stored: storedValidation,
    currentRevision: seedRevision,
    schemaDigest: session?.schema_digest ?? null,
    dirty,
  });
  const validation =
    boundValidation.status === 'stale' ? boundValidation : validationOverride ?? boundValidation;
  const listStatus = projectsQuery.isLoading || (projectId !== null && sessionsQuery.isLoading)
    ? 'loading'
    : projectsQuery.isError || sessionsQuery.isError
      ? 'error'
      : !projectId || (sessionsQuery.data?.sessions.length ?? 0) === 0
        ? 'empty'
        : 'ready';

  return (
    <div className="flex min-w-0 flex-1 flex-col gap-3">
      <div className="flex min-w-0 flex-col gap-2 sm:flex-row sm:items-end">
        <label className="min-w-0 flex-1 text-sm font-medium" htmlFor="configuration-draft-project">
          Project
          <select
            id="configuration-draft-project"
            className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
            value={projectId ?? ''}
            onChange={(event) => requestProject(event.target.value)}
          >
            {projects.length === 0 ? (
              <option value="">No projects yet</option>
            ) : (
              projects.map((project) => (
                <option key={project.project_id} value={project.project_id}>
                  {project.display_name}
                </option>
              ))
            )}
          </select>
        </label>
        <Button
          type="button"
          variant="primary"
          disabled={!projectId || creating}
          onClick={() => requestCreate()}
        >
          New draft
        </Button>
        {projectsQuery.data?.next_cursor && (
          <Button
            type="button"
            variant="secondary"
            onClick={() => void loadMoreProjects()}
          >
            More projects
          </Button>
        )}
      </div>
      {projectsQuery.isError && (
        <p role="alert">{classifyAuthoringError(projectsQuery.error).message}</p>
      )}
      {!projectsQuery.isLoading && !projectsQuery.isError && projects.length === 0 && (
        <p role="status">No projects are available for this login.</p>
      )}
      <ConfigurationDraftWorkbench
        workflowName={workflowName}
        projectName={projectName}
        schema={missing || draftLoading ? null : schema}
        schemaIdentity={session?.schema_digest ?? `local:${workflowId}`}
        schemaMessage={
          missing
            ? 'This draft is not available.'
            : draftLoading
              ? 'Loading this draft.'
              : sessionQuery.isError
                ? classifyAuthoringError(sessionQuery.error).message
                : schemaMessage
        }
        list={{
          status: listStatus,
          drafts: (sessionsQuery.data?.sessions ?? []).map((item) =>
            listLabel(item, projectName ?? 'Project', workflowName),
          ),
          pageLimit: PAGE_LIMIT,
          hasPrevious: cursors.length > 1,
          hasNext: Boolean(sessionsQuery.data?.next_cursor),
          errorMessage: projectsQuery.isError
            ? classifyAuthoringError(projectsQuery.error).message
            : sessionsQuery.isError
              ? classifyAuthoringError(sessionsQuery.error).message
              : undefined,
        }}
        selectedSid={selectedSid}
        draftExpired={Boolean(session && expired && !missing)}
        connection={connection}
        capability={capability}
        proposal={proposal}
        explanation={turn ? explanationView(turn, capability) : null}
        validation={draftOpen || storedValidation || validationOverride ? validation : { status: 'idle' }}
        requestStatus={requestStatus}
        writesEnabled={writesEnabled}
        apiConnected
        saveUnavailableReason={
          session && expired
            ? 'This draft has expired, so it cannot be saved.'
            : 'Open one of your drafts before this action.'
        }
        documentSeed={seed}
        externalDirty={inputsDirty}
        confirmedSnapshot={confirmedSnapshot}
        serverReloadEpoch={reloadEpoch}
        remoteConflict={remoteConflict}
        managedInputs={
          draftOpen ? (
            <ManagedDraftInputs
              samples={samples}
              selectedRevisionIds={sampleIds}
              unlistedCount={unlistedCount}
              samplesLoading={samplesQuery.isLoading}
              samplesError={
                samplesQuery.isError ? classifyAuthoringError(samplesQuery.error).message : null
              }
              onToggleSample={(revisionId, selected) =>
                setSampleIds((current) =>
                  selected
                    ? current.includes(revisionId)
                      ? current
                      : [...current, revisionId]
                    : current.filter((id) => id !== revisionId),
                )
              }
              onDropUnlisted={() => setSampleIds((current) => current.filter((id) => listedSampleIds.has(id)))}
              profiles={profiles}
              referenceRevisionId={referenceId}
              referenceLoading={referenceQuery.isLoading}
              referenceError={
                referenceQuery.isError
                  ? 'The reference list could not be loaded.'
                  : null
              }
              referenceMissing={referenceMissing}
              onSelectReference={setReferenceId}
              onRefreshReferences={() => void referenceQuery.refetch()}
              preservedSelectionCount={preservedRef.current.selections.length}
              sampleRowsAbsent={workflowSampleRowsAbsent(activeDraft?.samples)}
            />
          ) : null
        }
        actionNotice={notice}
        onSelectDraft={selectDraft}
        onListPrevious={() => setCursors((current) => current.slice(0, -1))}
        onListNext={() => {
          const next = sessionsQuery.data?.next_cursor;
          if (typeof next === 'string') setCursors((current) => [...current, next]);
        }}
        onRetryList={() => void sessionsQuery.refetch()}
        onSave={(document) => void saveDraft(document)}
        onDirtyChange={setDirty}
        onDiscarded={() => {
          if (!activeDraft) return;
          setSampleIds([...(activeDraft.sample_revision_ids ?? [])]);
          setReferenceId(activeDraft.reference_profile_revision_id ?? null);
        }}
        onReloadRemote={reloadRemote}
        onKeepRemote={() => {
          if (!session) return;
          setIgnoredRevision(session.draft_revision);
          setRemoteConflict(false);
        }}
        onConnect={draftOpen ? () => void connect() : undefined}
        onClose={connection.state === 'available' ? () => void closeConnection() : undefined}
        onRequest={(question) => void requestSuggestion(question)}
        onCheckRequest={() => void checkRequest()}
        onAdopt={() => void adopt()}
        onReject={() => void reject()}
        onRefreshProposal={() => void refreshProposal()}
        onValidate={() => void validateDraft()}
        operationPending={operationPending}
        externalPrompt={pendingProjectId !== null || pendingCreate}
        onExternalDiscard={confirmExternalNavigation}
        onExternalCancel={() => {
          setPendingProjectId(null);
          setPendingCreate(false);
        }}
      />
    </div>
  );

  function requestProject(nextProjectId: string) {
    if (!nextProjectId || nextProjectId === projectId) return;
    if (dirty) {
      setPendingProjectId(nextProjectId);
      return;
    }
    applyProject(nextProjectId);
  }

  function requestCreate() {
    if (dirty) {
      setPendingCreate(true);
      return;
    }
    void createDraft();
  }

  function confirmExternalNavigation() {
    const nextProject = pendingProjectId;
    const create = pendingCreate;
    setPendingProjectId(null);
    setPendingCreate(false);
    if (nextProject) applyProject(nextProject);
    if (create) void createDraft();
  }

  function applyProject(nextProjectId: string) {
    abandonOperation();
    pendingRemote.current = null;
    scope.current += 1;
    setProjectId(nextProjectId);
    setCursors([null]);
    setSelectedSid(null);
    setExtraSamples([]);
    setTurn(null);
    setStoredValidation(null);
    setValidationOverride(null);
    pendingTurn.current = null;
    pendingDecision.current = null;
  }

  function selectDraft(sid: string) {
    abandonOperation();
    pendingRemote.current = null;
    scope.current += 1;
    setSelectedSid(sid);
    setSeed(null);
    setSeedRevision(null);
    setIgnoredRevision(null);
    setRemoteConflict(false);
    setTurn(null);
    setDecision(null);
    setStoredValidation(null);
    setValidationOverride(null);
    setRequestStatus('idle');
    setNotice(null);
    setBridgeBlocked(false);
    pendingTurn.current = null;
    pendingDecision.current = null;
  }

  async function loadMoreProjects() {
    const next = projectsQuery.data?.next_cursor;
    if (typeof next !== 'string') return;
    const page = await authoringListProjects({ cursor: next, limit: 50 });
    setExtraProjects((current) => {
      const seen = new Set(current.map((project) => project.project_id));
      return [...current, ...page.projects.filter((project) => !seen.has(project.project_id))];
    });
  }
}

function withDecision(proposal: ProposalView | null, decision: DecisionName | null): ProposalView | null {
  if (!proposal || !decision) return proposal;
  if (decision === 'unconfirmed') return { ...proposal, state: 'unconfirmed', adoptable: false };
  if (decision === 'adopted') return { ...proposal, state: 'adopted', adoptable: false };
  return { ...proposal, state: 'rejected', adoptable: false };
}

function newClientTurnId(): string {
  return crypto.randomUUID().replace(/-/g, '');
}
