import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { useParams } from 'react-router-dom';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import { ApiError } from '../../api/fetcher';
import type {
  AuthoringSessionResponse,
  AuthoringTurnCreateRequest,
  AuthoringTurnResponse,
  AuthoringValidationResponse,
} from '../../api/generated/models';
import { ConfigurationDraftsPage } from './configuration-drafts';
import { createDraftEditorContract } from '../../features/configuration-drafts/testFixtures';
import { WORKFLOW_ID } from '../../features/input-workbench/test-fixtures';
import { renderWithRouter } from '../../test/test-utils';

const PROJECT = `prj_${'a'.repeat(32)}`;
const SID = `aas_${'b'.repeat(32)}`;
const OTHER_SID = `aas_${'c'.repeat(32)}`;
const DIGEST = 'ab'.repeat(32);
const PROPOSAL = `prop_${'d'.repeat(32)}`;
const TURN = `atr_${'e'.repeat(32)}`;
const SECRET_PATH = '/tmp/should-not-render/samples.tsv';

const authoring = vi.hoisted(() => ({
  listProjects: vi.fn(),
  listSamples: vi.fn(),
  listSessions: vi.fn(),
  createSession: vi.fn(),
  getSession: vi.fn(),
  saveDraft: vi.fn(),
  openConnection: vi.fn(),
  closeConnection: vi.fn(),
  createTurn: vi.fn(),
  getTurn: vi.fn(),
  adopt: vi.fn(),
  reject: vi.fn(),
  validate: vi.fn(),
}));

vi.mock('../../api/generated/authoring/authoring', () => ({
  authoringListProjects: authoring.listProjects,
  authoringListProjectSamples: authoring.listSamples,
  authoringListSessions: authoring.listSessions,
  authoringCreateSession: authoring.createSession,
  authoringGetSession: authoring.getSession,
  authoringSaveDraft: authoring.saveDraft,
  authoringOpenConnection: authoring.openConnection,
  authoringCloseConnection: authoring.closeConnection,
  authoringCreateTurn: authoring.createTurn,
  authoringGetTurn: authoring.getTurn,
  authoringAdoptProposal: authoring.adopt,
  authoringRejectProposal: authoring.reject,
  authoringValidateDraft: authoring.validate,
}));

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: vi.fn(async () => ({
    ok: true,
    workflow_id: WORKFLOW_ID,
    schema: createDraftEditorContract(),
    issues: [],
  })),
  listCompatibleReferenceProfiles: vi.fn(async () => ({
    ok: true,
    workflow_id: WORKFLOW_ID,
    profiles: [],
    issues: [],
  })),
}));

vi.mock('../../api/authClient', async () => {
  const actual = await vi.importActual<typeof import('../../api/authClient')>(
    '../../api/authClient',
  );
  return {
    ...actual,
    fetchSessionState: vi.fn(async () => ({
      ok: true,
      setup_required: false,
      authenticated: true,
      principal: { user_id: 'user-1', username: 'ada', role: 'admin' },
    })),
  };
});

vi.mock('@uiw/react-codemirror', () => ({
  default: ({
    value,
    onChange,
  }: {
    value: string;
    onChange: (text: string) => void;
  }) => (
    <textarea
      aria-label="Advanced config YAML"
      value={value}
      onChange={(event) => onChange(event.target.value)}
    />
  ),
}));

describe('configuration draft API wiring', () => {
  beforeEach(() => {
    vi.clearAllMocks();
    authoring.listProjects.mockResolvedValue({
      ok: true,
      projects: [{ project_id: PROJECT, display_name: 'Alpha project' }],
    });
    authoring.listSamples.mockResolvedValue({ ok: true, samples: [] });
    authoring.listSessions.mockResolvedValue({
      ok: true,
      sessions: [summary(SID), summary(OTHER_SID, 2)],
    });
    authoring.getSession.mockImplementation(async (sid: string) =>
      sid === OTHER_SID ? draftSession(OTHER_SID, { limit: 8 }) : draftSession(SID, { limit: 4 }),
    );
    authoring.getTurn.mockRejectedValue(new ApiError(404, 'NOT_FOUND', 'missing'));
  });

  it('restores a draft after sign-in without opening the assistant', async () => {
    renderPage();
    await openDraft(SID);
    expect(await screen.findByLabelText('Local limit')).toHaveValue(4);
    expect(screen.getByText(/draft is still here after signing in/i)).toBeInTheDocument();
    expect(screen.queryByText(SECRET_PATH)).not.toBeInTheDocument();
    expect(document.querySelector('input[type="file"]')).toBeNull();
    expect(authoring.openConnection).not.toHaveBeenCalled();
    expect(authoring.createTurn).not.toHaveBeenCalled();
    expect(authoring.createSession).not.toHaveBeenCalled();
  });

  it('keeps unsaved edits when save conflicts or the response is lost', async () => {
    const user = userEvent.setup();
    renderPage();
    await openDraft(SID);
    await editLimit(user, '5');
    authoring.saveDraft.mockRejectedValueOnce(new ApiError(409, 'STALE_DRAFT', 'raw server text'));
    await user.click(screen.getByRole('button', { name: 'Save edits' }));
    expect(await screen.findByText(/draft changed on the server/i)).toBeInTheDocument();
    expect(screen.getByText(/nothing is overwritten/i)).toBeInTheDocument();
    expect(screen.getByLabelText('Local limit')).toHaveValue(5);
    expect(authoring.saveDraft).toHaveBeenCalledTimes(1);
    expect(authoring.saveDraft.mock.calls[0]?.[1]).toMatchObject({ expected_revision: 1 });

    authoring.saveDraft.mockRejectedValueOnce(new TypeError('socket lost'));
    authoring.getSession.mockResolvedValueOnce(draftSession(SID, { limit: 5 }));
    await user.click(screen.getByRole('button', { name: 'Save edits' }));
    expect(
      await screen.findByText(
        'The save response was not confirmed. These edits are still unsaved.',
      ),
    ).toBeInTheDocument();
    expect(screen.getByLabelText('Local limit')).toHaveValue(5);
    expect(authoring.saveDraft).toHaveBeenCalledTimes(2);
  });

  it('keeps manual save available when the assistant bridge is unavailable', async () => {
    const user = userEvent.setup();
    renderPage();
    await openDraft(SID);
    authoring.openConnection.mockRejectedValueOnce(
      new ApiError(503, 'BRIDGE_UNAVAILABLE', 'raw bridge text'),
    );
    await user.click(screen.getByRole('button', { name: 'Connect assistant' }));
    expect(await screen.findByText('Assistant unavailable')).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Validate draft' })).toBeEnabled();
    await editLimit(user, '6');
    expect(screen.getByRole('button', { name: 'Save edits' })).toBeEnabled();
  });

  it('adopts and rejects with the server digest and does not adopt an explanation', async () => {
    const user = userEvent.setup();
    renderPage();
    await openDraft(SID);
    await connectReady();
    await user.click(screen.getByRole('radio', { name: 'Explain a validation issue' }));
    await user.click(screen.getByRole('radio', { name: 'Select a sample revision.' }));
    authoring.createTurn.mockResolvedValueOnce(explanationTurn());
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(await screen.findByText(/at least one sample revision/i)).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Adopt suggestion' })).not.toBeInTheDocument();
    const explainBody = authoring.createTurn.mock.calls[0]?.[1] as AuthoringTurnCreateRequest;
    expect(explainBody.question).toEqual({
      intent: 'explain_issues',
      field_tokens: [],
      preferences: [],
    });
    expect(JSON.stringify(explainBody)).not.toMatch(/SAMPLE_SELECTION_REQUIRED/);

    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    await user.click(screen.getByRole('checkbox', { name: 'Example limit' }));
    await user.type(screen.getByRole('textbox', { name: /Preferred value/ }), '0');
    authoring.createTurn.mockResolvedValueOnce(proposalTurn());
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(await screen.findByRole('button', { name: 'Adopt suggestion' })).toBeEnabled();
    authoring.adopt.mockResolvedValueOnce({
      ok: true,
      replayed: false,
      adoption_receipt: {
        proposal_id: PROPOSAL,
        base_revision: 1,
        result_revision: 2,
        changed: true,
        adopted_at: '2026-11-10T00:00:00Z',
        draft_expires_at: '2026-11-17T00:00:00Z',
      },
      draft_state: 'active',
      draft_revision: 2,
      expires_at: '2026-11-17T00:00:00Z',
      purge_after: '2026-11-24T00:00:00Z',
    });
    authoring.getTurn.mockResolvedValueOnce({ ...proposalTurn(), proposal_state: 'adopted', adoptable: false });
    authoring.getSession.mockResolvedValueOnce({
      ...draftSession(SID, { limit: 4 }, { draft_revision: 2 }),
      connection: { state: 'ready' },
      requires_new_connection: false,
    });
    await user.click(screen.getByRole('button', { name: 'Adopt suggestion' }));
    await waitFor(() => expect(authoring.adopt).toHaveBeenCalledTimes(1));
    expect(authoring.adopt.mock.calls[0]?.[2]).toEqual({
      expected_revision: 1,
      proposal_digest: DIGEST,
    });
    expect(await screen.findByText('Suggestion adopted.')).toBeInTheDocument();

    authoring.createTurn.mockResolvedValueOnce(proposalTurn());
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(await screen.findByRole('button', { name: 'Reject suggestion' })).toBeEnabled();
    authoring.reject.mockResolvedValueOnce({
      ...proposalTurn(),
      proposal_state: 'rejected',
      adoptable: false,
    });
    await user.click(screen.getByRole('button', { name: 'Reject suggestion' }));
    expect(await screen.findByText(/suggestion rejected/i)).toBeInTheDocument();
    expect(authoring.reject.mock.calls[0]?.[2]).toEqual({ proposal_digest: DIGEST });
  });

  it('replays the same suggestion request after a lost response and ignores a late one', async () => {
    const user = userEvent.setup();
    const bodies: AuthoringTurnCreateRequest[] = [];
    authoring.createTurn.mockImplementation(async (_sid: string, body: AuthoringTurnCreateRequest) => {
      bodies.push(structuredClone(body));
      if (bodies.length === 1) throw new TypeError('lost');
      return proposalTurn();
    });
    renderPage();
    await openDraft(SID);
    await connectReady();
    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    await user.click(screen.getByRole('checkbox', { name: 'Example limit' }));
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(await screen.findByText(/nothing was marked complete/i)).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Check this request' }));
    await screen.findByRole('button', { name: 'Adopt suggestion' });
    expect(bodies).toHaveLength(2);
    expect(bodies[1]).toEqual(bodies[0]);
  });

  it('ignores a suggestion response that arrives after switching drafts', async () => {
    const user = userEvent.setup();
    let resolveLate: (value: AuthoringTurnResponse) => void = () => undefined;
    authoring.createTurn.mockImplementation(
      () =>
        new Promise<AuthoringTurnResponse>((resolve) => {
          resolveLate = resolve;
        }),
    );
    renderPage();
    await openDraft(SID);
    await connectReady();
    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    await user.click(screen.getByRole('checkbox', { name: 'Example limit' }));
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    await user.click(screen.getByRole('button', { name: /revision 2/ }));
    resolveLate(proposalTurn(77));
    expect(await screen.findByLabelText('Local limit')).toHaveValue(8);
    expect(screen.queryByText('77')).not.toBeInTheDocument();
  });

  it('does not present an old validation result as current after an edit', async () => {
    const user = userEvent.setup();
    authoring.validate.mockResolvedValueOnce({
      draft_revision: 1,
      ok: true,
      snapshot: null,
    } satisfies AuthoringValidationResponse);
    renderPage();
    await openDraft(SID);
    await user.click(screen.getByRole('button', { name: 'Validate draft' }));
    expect(await screen.findByText(/no execution snapshot/i)).toBeInTheDocument();
    await editLimit(user, '3');
    expect(await screen.findByText(/not for the current draft/i)).toBeInTheDocument();
    expect(screen.queryByText(/no execution snapshot/i)).not.toBeInTheDocument();
    expect(screen.queryByText(/validation finished/i)).not.toBeInTheDocument();
  });

  it('distinguishes an expired draft from a missing one', async () => {
    const user = userEvent.setup();
    authoring.listSessions.mockResolvedValue({
      ok: true,
      sessions: [
        { ...summary(SID), draft_state: 'expired' },
        summary(OTHER_SID, 2),
      ],
    });
    authoring.getSession.mockImplementation(async (sid: string) => {
      if (sid === SID) {
        return {
          ...draftSession(SID, { limit: 4 }),
          draft_state: 'expired',
          draft: null,
        };
      }
      throw new ApiError(404, 'NOT_FOUND', 'raw missing text');
    });
    renderPage();
    await user.click(await screen.findByRole('button', { name: /Expired/ }));
    expect(await screen.findByText(/configuration is hidden/i)).toBeInTheDocument();
    expect(screen.queryByLabelText('Local limit')).not.toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: /revision 2/ }));
    expect(await screen.findByText(/not available/i)).toBeInTheDocument();
    expect(screen.queryByText(/raw missing text/i)).not.toBeInTheDocument();
  });
});

function renderPage() {
  return renderWithRouter(
    [
      {
        path: '/workflows/:workflowId/configuration-drafts',
        element: <DraftRoute />,
      },
    ],
    { initialEntries: [`/workflows/${WORKFLOW_ID}/configuration-drafts`] },
  );
}

function DraftRoute() {
  const params = useParams();
  return <ConfigurationDraftsPage workflowId={params.workflowId ?? ''} />;
}

async function openDraft(sid: string) {
  const user = userEvent.setup();
  const revision = sid === OTHER_SID ? 2 : 1;
  await user.click(await screen.findByRole('button', { name: new RegExp(`revision ${revision}`) }));
  await screen.findByLabelText('Local limit');
}

async function editLimit(user: ReturnType<typeof userEvent.setup>, value: string) {
  const limit = await screen.findByLabelText('Local limit');
  await user.clear(limit);
  await user.type(limit, value);
  expect(
    await screen.findByText(
      'You have unsaved edits. They are not saved, adopted, or validated.',
    ),
  ).toBeInTheDocument();
}

async function connectReady() {
  const user = userEvent.setup();
  authoring.openConnection.mockImplementation(async (sid: string) => ({
    ...(await authoring.getSession(sid)),
    connection: { state: 'ready' },
    requires_new_connection: false,
  }));
  await user.click(screen.getByRole('button', { name: 'Connect assistant' }));
  expect(await screen.findByText('Assistant available')).toBeInTheDocument();
}

function summary(sid: string, revision = 1) {
  return {
    sid,
    project_id: PROJECT,
    workflow_id: WORKFLOW_ID,
    draft_state: 'active' as const,
    draft_revision: revision,
    expires_at: '2026-11-17T00:00:00Z',
    purge_after: '2026-11-24T00:00:00Z',
    connection: { state: 'closed' as const },
    requires_new_connection: true,
  };
}

function draftSession(
  sid: string,
  config: { limit: number },
  overrides: Partial<AuthoringSessionResponse> = {},
): AuthoringSessionResponse {
  const other = sid === OTHER_SID;
  return {
    ok: true,
    sid,
    project_id: PROJECT,
    workflow_id: WORKFLOW_ID,
    draft_state: 'active',
    draft_revision: other ? 2 : 1,
    expires_at: '2026-11-17T00:00:00Z',
    purge_after: '2026-11-24T00:00:00Z',
    schema_digest: 'a'.repeat(64),
    projection_digest: other ? 'c'.repeat(64) : 'b'.repeat(64),
    draft: {
      config,
      options: {},
      samples: SECRET_PATH,
      sample_revision_ids: [],
      input_selections: [],
      reference_profile_revision_id: null,
    },
    authoring_capability: {
      available: true,
      projection_version: 'hw-a1/1',
      allowed_intents: ['suggest', 'explain_issues'],
      allowed_issue_codes: [
        { code: 'SAMPLE_SELECTION_REQUIRED', label: 'Select a sample revision.' },
      ],
      fields: [
        {
          token: other ? 'resource.other' : 'resource.limit',
          label: other ? 'Other limit' : 'Example limit',
          type: 'integer',
          constraints: { type: 'integer', minimum: 0, maximum: 8 },
        },
      ],
    },
    connection: { state: 'closed' },
    requires_new_connection: true,
    ...overrides,
  };
}

function proposalTurn(after = 2): AuthoringTurnResponse {
  return {
    ok: true,
    turn_id: TURN,
    client_turn_id: 'client',
    status: 'completed',
    result_kind: 'proposal',
    proposal_id: PROPOSAL,
    proposal_digest: DIGEST,
    adoptable: true,
    proposal_state: 'pending',
    proposal: {
      changes: [
        {
          field_token: 'resource.limit',
          value: after,
          before: { present: true, token: 'resource.limit', value: 4 },
          after: { present: true, token: 'resource.limit', value: after },
        },
      ],
      reason_codes: [],
    },
  };
}

function explanationTurn(): AuthoringTurnResponse {
  return {
    ok: true,
    turn_id: TURN,
    client_turn_id: 'client',
    status: 'completed',
    result_kind: 'explanation',
    adoptable: false,
    explanation: { reason_codes: ['SAMPLE_SELECTION_REQUIRED'] },
  };
}
