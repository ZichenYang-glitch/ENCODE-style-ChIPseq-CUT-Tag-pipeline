import { screen, waitFor } from '@testing-library/react';
import userEvent, { type UserEvent } from '@testing-library/user-event';
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { ApiError } from '../../api/fetcher';
import { appRoutes } from '../../app/router';
import { sampleRevisionSetsEqual } from '../../features/experiment-binding/handoff';
import { WORKFLOW_ID, createAuthoringSchemaFixture } from '../../features/input-workbench/test-fixtures';
import { renderWithRouter } from '../../test/test-utils';

const PROJECT = `prj_${'1'.repeat(32)}`;
const SID = `aas_${'1'.repeat(32)}`;
const REV_A = `smpr_${'a'.repeat(32)}`;
const REV_B = `smpr_${'b'.repeat(32)}`;
const REV_OTHER = `smpr_${'c'.repeat(32)}`;
const EXPR = `expr_${'e'.repeat(32)}`;
const REFERENCE = `refpr_${'1'.repeat(32)}`;
const SNAPSHOT_ID = 'vsnap_0123456789abcdef0123456789abcdef';
const SAVED_SAMPLE_IDS = [REV_A, REV_B];

const PROFILE = {
  profile_id: `refp_${'1'.repeat(32)}`,
  revision_id: REFERENCE,
  revision_number: 1,
  display_name: 'Synthetic reference',
  organism: 'synthetic',
  assembly: 'test',
  identity_sha256: 'ab'.repeat(32),
};

interface ValidateBody {
  config?: { threads?: number };
  samples?: Array<{ sample?: string }>;
  project_id?: string;
  sample_revision_ids?: string[];
  experiment_revision_id?: string;
  reference_profile_revision_id?: string;
}

const generatedMocks = vi.hoisted(() => ({
  createRun: vi.fn(),
  getWorkflowSchema: vi.fn(),
  validateWorkflow: vi.fn(),
  listCompatibleReferenceProfiles: vi.fn(),
  openConnection: vi.fn(),
  createSession: vi.fn(),
  saveDraft: vi.fn(),
  validateDraft: vi.fn(),
  getSession: vi.fn(),
  listProjects: vi.fn(),
  listSessions: vi.fn(),
  listSamples: vi.fn(),
  listInputRevisions: vi.fn(async (projectId: string) => ({
    ok: true as const,
    project_id: projectId,
    revisions: [],
    next_cursor: null,
  })),
  listReceipts: vi.fn(async (sid: string) => ({
    ok: true as const,
    sid,
    draft_state: 'active' as const,
    receipts: [],
    next_cursor: null,
  })),
  catalogListExperiments: vi.fn(),
}));

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: generatedMocks.getWorkflowSchema,
  validateWorkflow: generatedMocks.validateWorkflow,
  listCompatibleReferenceProfiles: generatedMocks.listCompatibleReferenceProfiles,
  getWorkflow: vi.fn(),
  listWorkflows: vi.fn(),
}));

vi.mock('../../api/generated/runs/runs', () => ({
  createRun: generatedMocks.createRun,
  cancelRun: vi.fn(),
  getRun: vi.fn(),
  listRunEvents: vi.fn(),
  listRunLogs: vi.fn(),
  startRun: vi.fn(),
}));

vi.mock('../../api/generated/catalog/catalog', () => ({
  catalogListExperiments: generatedMocks.catalogListExperiments,
}));

vi.mock('../../api/generated/authoring/authoring', () => ({
  authoringListProjects: generatedMocks.listProjects,
  authoringListProjectSamples: generatedMocks.listSamples,
  authoringListProjectInputRevisions: generatedMocks.listInputRevisions,
  authoringListSessions: generatedMocks.listSessions,
  authoringCreateSession: generatedMocks.createSession,
  authoringGetSession: generatedMocks.getSession,
  authoringSaveDraft: generatedMocks.saveDraft,
  authoringOpenConnection: generatedMocks.openConnection,
  authoringCloseConnection: vi.fn(),
  authoringCreateTurn: vi.fn(),
  authoringGetTurn: vi.fn(),
  authoringListReceipts: generatedMocks.listReceipts,
  authoringAdoptProposal: vi.fn(),
  authoringRejectProposal: vi.fn(),
  authoringValidateDraft: generatedMocks.validateDraft,
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

let frozenBindings = 0;

function snapshot(binding: { status: 'unresolved' } | {
  status: 'experiment_snapshot';
  experiment_revision_id: string;
  payload_digest: string;
}) {
  return {
    ok: true,
    workflow_id: WORKFLOW_ID,
    value: null,
    snapshot: {
      snapshot_id: SNAPSHOT_ID,
      workflow_id: WORKFLOW_ID,
      schema_version: '1.1.0',
      adapter_version: '0.3.0',
      payload_digest: 'a'.repeat(64),
      validated_at: '2026-07-14T00:00:00.000Z',
      expires_at: '2026-07-14T00:30:00.000Z',
      project_id: PROJECT,
      binding_mode: 'bound_v1',
      provenance: 'resolved',
      sample_revision_ids: [REV_A, REV_B],
      binding_digest: 'b'.repeat(64),
      input_binding: {
        mode: 'compatibility_unresolved_v1',
        adapter_contract_version: null,
        digest: 'c'.repeat(64),
        fully_managed: false,
        input_uses: [],
      },
      reference_profile: PROFILE,
      experiment_binding: binding,
    },
    issues: [],
  };
}

function savedSession() {
  return {
    ok: true,
    sid: SID,
    project_id: PROJECT,
    workflow_id: WORKFLOW_ID,
    draft_state: 'active',
    draft_revision: 3,
    expires_at: '2026-11-17T00:00:00Z',
    purge_after: '2026-11-24T00:00:00Z',
    schema_digest: 'a'.repeat(64),
    projection_digest: 'b'.repeat(64),
    draft: {
      config: { threads: 4, outdir: 'results' },
      options: { strict_inputs: false },
      samples: [{ sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' }],
      sample_revision_ids: SAVED_SAMPLE_IDS,
      input_selections: [],
      reference_profile_revision_id: REFERENCE,
    },
    authoring_capability: {
      available: true,
      projection_version: 'hw-a1/1',
      allowed_intents: ['suggest', 'explain_issues'],
      allowed_issue_codes: [],
      fields: [
        {
          token: 'resource.threads',
          label: 'Threads',
          type: 'integer',
          constraints: { type: 'integer', minimum: 1, maximum: 16 },
        },
      ],
    },
    connection: { state: 'closed' },
    requires_new_connection: true,
  };
}

beforeAll(async () => {
  await import('../workflows/configuration-drafts');
  await import('../workflows/new-run');
});

describe('experiment binding entry', () => {
  beforeEach(() => {
    frozenBindings = 0;
    generatedMocks.getWorkflowSchema.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      schema: createAuthoringSchemaFixture(),
      issues: [],
    });
    generatedMocks.listCompatibleReferenceProfiles.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      profiles: [PROFILE],
      issues: [],
    });
    generatedMocks.listProjects.mockResolvedValue({
      ok: true,
      projects: [{ project_id: PROJECT, display_name: 'Synthetic project' }],
    });
    generatedMocks.listSessions.mockResolvedValue({
      ok: true,
      sessions: [{
        sid: SID,
        project_id: PROJECT,
        workflow_id: WORKFLOW_ID,
        draft_state: 'active',
        draft_revision: 3,
        expires_at: '2026-11-17T00:00:00Z',
        purge_after: '2026-11-24T00:00:00Z',
        connection: { state: 'closed' },
        requires_new_connection: true,
      }],
      next_cursor: null,
    });
    generatedMocks.listSamples.mockResolvedValue({
      ok: true,
      project_id: PROJECT,
      samples: [],
      next_cursor: null,
    });
    generatedMocks.getSession.mockResolvedValue(savedSession());
    generatedMocks.catalogListExperiments.mockResolvedValue({
      ok: true,
      project_id: PROJECT,
      revision_view: 'current',
      experiments: [{
        experiment_id: `exp_${'1'.repeat(32)}`,
        experiment_revision_id: EXPR,
        revision_number: 1,
        display_name: 'Synthetic experiment',
        samples: [
          { sample_revision_id: REV_A },
          { sample_revision_id: REV_B },
        ],
      }],
      next_cursor: null,
    });
    generatedMocks.validateWorkflow.mockClear();
    generatedMocks.validateWorkflow.mockImplementation(async (_workflowId: string, body: ValidateBody) => {
      if (body.experiment_revision_id) {
        const ids = body.sample_revision_ids ?? [];
        if (
          body.experiment_revision_id !== EXPR ||
          !sampleRevisionSetsEqual(ids, SAVED_SAMPLE_IDS)
        ) {
          throw new ApiError(
            409,
            'catalog_sample_set_mismatch',
            'The sample revisions do not match the experiment.',
            [{
              code: 'catalog_sample_set_mismatch',
              message: 'The sample revisions do not match the experiment.',
              severity: 'error',
            }],
          );
        }
        frozenBindings += 1;
        return snapshot({
          status: 'experiment_snapshot',
          experiment_revision_id: body.experiment_revision_id,
          payload_digest: 'd'.repeat(64),
        });
      }
      return snapshot({ status: 'unresolved' });
    });
    generatedMocks.createRun.mockResolvedValue({
      ok: false,
      workflow_id: WORKFLOW_ID,
      run: null,
      issues: [],
    });
    generatedMocks.saveDraft.mockReset();
    generatedMocks.validateDraft.mockReset();
    generatedMocks.createRun.mockClear();
  });

  async function openCarriedValidation(user: UserEvent) {
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/configuration-drafts`],
    });
    await user.click(await screen.findByRole('button', { name: /Synthetic project/ }));
    expect(generatedMocks.validateDraft).not.toHaveBeenCalled();
    expect(generatedMocks.saveDraft).not.toHaveBeenCalled();
    expect(generatedMocks.validateWorkflow).not.toHaveBeenCalled();
    await user.click(await screen.findByRole('link', { name: 'Open experiment validation' }));
    expect(await screen.findByTestId('carried-project-id')).toHaveTextContent(PROJECT);
    expect(screen.getByTestId('carried-reference-profile')).toHaveTextContent(REFERENCE);
    await waitFor(() => {
      expect(screen.getByTestId('draft-review-json')).toHaveTextContent('"threads": 4');
      expect(screen.getByTestId('draft-review-json')).toHaveTextContent('S1');
    });
    expect(screen.getByLabelText('Submitted sample revision 1')).toHaveValue(REV_A);
    expect(screen.getByLabelText('Submitted sample revision 2')).toHaveValue(REV_B);
    expect(screen.getByTestId('experiment-binding-status')).toHaveTextContent(
      'No experiment binding. Choosing a revision or opening this page does not create one.',
    );
    await waitFor(() => {
      expect(screen.getByRole('button', { name: 'Validate current inputs' })).toBeEnabled();
    });
  }

  it('carries the saved draft, freezes a reversed set, and refuses a different revision', async () => {
    const user = userEvent.setup();
    await openCarriedValidation(user);

    await user.click(screen.getByRole('button', { name: 'Reverse sample revision order' }));
    expect(screen.getByLabelText('Submitted sample revision 1')).toHaveValue(REV_B);
    expect(screen.getByLabelText('Submitted sample revision 2')).toHaveValue(REV_A);

    const choice = screen.getByRole('combobox', { name: 'Experiment revision' });
    await waitFor(() => expect(choice).toBeEnabled());
    await user.selectOptions(choice, 'Synthetic experiment');
    expect(screen.getByTestId('experiment-selection-preview')).toHaveTextContent(
      'Synthetic experiment pins the same sample revisions. Order does not matter. Validation can freeze this revision.',
    );
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    expect(await screen.findByText(`Bound to experiment revision ${EXPR} for this validation snapshot.`)).toBeInTheDocument();
    const boundBody = generatedMocks.validateWorkflow.mock.calls.at(-1)?.[1] as ValidateBody;
    expect(boundBody.project_id).toBe(PROJECT);
    expect(boundBody.sample_revision_ids).toEqual([REV_B, REV_A]);
    expect(boundBody.experiment_revision_id).toBe(EXPR);
    expect(boundBody.config?.threads).toBe(4);
    expect(boundBody.reference_profile_revision_id).toBe(REFERENCE);
    expect(boundBody.samples?.[0]?.sample).toBe('S1');
    expect(frozenBindings).toBe(1);

    await user.selectOptions(choice, '');
    expect(screen.getByTestId('experiment-binding-status')).toHaveTextContent(
      'No experiment binding. Choosing a revision or opening this page does not create one.',
    );
    await waitFor(() => {
      expect(screen.getByTestId('create-validated-run-button')).toBeDisabled();
      expect(screen.getByTestId('validation-feedback')).not.toHaveTextContent(
        'This exact draft can create one run.',
      );
    });
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    expect(await screen.findByText('Not bound. This validation did not select an experiment revision.')).toBeInTheDocument();
    const unboundBody = generatedMocks.validateWorkflow.mock.calls.at(-1)?.[1] as ValidateBody;
    expect(unboundBody.experiment_revision_id).toBeUndefined();
    expect(unboundBody.sample_revision_ids).toEqual([REV_B, REV_A]);
    expect(frozenBindings).toBe(1);

    const first = screen.getByLabelText('Submitted sample revision 1');
    await user.clear(first);
    await user.type(first, REV_OTHER);
    await user.selectOptions(choice, 'Synthetic experiment');
    expect(screen.getByTestId('experiment-selection-preview')).toHaveTextContent(
      'Synthetic experiment does not pin this sample revision set. Validation will not freeze a binding.',
    );
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    expect(await screen.findByText('No experiment binding was saved.')).toBeInTheDocument();
    expect(await screen.findByText('catalog_sample_set_mismatch')).toBeInTheDocument();
    const rejectedBody = generatedMocks.validateWorkflow.mock.calls.at(-1)?.[1] as ValidateBody;
    expect(rejectedBody.sample_revision_ids).toEqual([REV_OTHER, REV_A]);
    expect(rejectedBody.experiment_revision_id).toBe(EXPR);
    expect(frozenBindings).toBe(1);
    expect(generatedMocks.createRun).not.toHaveBeenCalled();
    expect(generatedMocks.validateDraft).not.toHaveBeenCalled();
    expect(generatedMocks.saveDraft).not.toHaveBeenCalled();
    expect(SAVED_SAMPLE_IDS).toEqual([REV_A, REV_B]);
  }, 30_000);

  it('creates a run from the frozen snapshot id only', async () => {
    const user = userEvent.setup();
    await openCarriedValidation(user);
    await user.click(screen.getByRole('button', { name: 'Reverse sample revision order' }));
    const choice = screen.getByRole('combobox', { name: 'Experiment revision' });
    await waitFor(() => expect(choice).toBeEnabled());
    await user.selectOptions(choice, 'Synthetic experiment');
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    await screen.findByText(`Bound to experiment revision ${EXPR} for this validation snapshot.`);
    await user.click(screen.getByRole('button', { name: 'Create run from validated inputs' }));
    await waitFor(() => expect(generatedMocks.createRun).toHaveBeenCalledTimes(1));
    expect(generatedMocks.createRun).toHaveBeenCalledWith(WORKFLOW_ID, {
      snapshot_id: SNAPSHOT_ID,
    });
    expect(Object.keys(generatedMocks.createRun.mock.calls[0][1] as object)).toEqual(['snapshot_id']);
    expect(generatedMocks.validateDraft).not.toHaveBeenCalled();
  }, 30_000);

  it('keeps the carried sample sheet and experiment selector across tabs', async () => {
    const user = userEvent.setup();
    await openCarriedValidation(user);
    await user.click(screen.getByRole('button', { name: 'Reverse sample revision order' }));
    await user.click(screen.getByRole('tab', { name: 'Config' }));
    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('carried-project-id')).toHaveTextContent(PROJECT);
    expect(screen.getByLabelText('Submitted sample revision 1')).toHaveValue(REV_B);
    expect(screen.getByLabelText('Submitted sample revision 2')).toHaveValue(REV_A);
    expect(document.getElementById('experiment-revision-choice')).toBeInTheDocument();
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent('S1');
  });

  it('leaves a lost validation response unconfirmed', async () => {
    generatedMocks.validateWorkflow.mockRejectedValue(new TypeError('network down'));
    const user = userEvent.setup();
    await openCarriedValidation(user);
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    expect(await screen.findByText('VALIDATION_UNAVAILABLE')).toBeInTheDocument();
    expect(screen.getByTestId('experiment-binding-status')).toHaveTextContent(
      'The validation result is unconfirmed.',
    );
    expect(screen.getByTestId('experiment-binding-status')).not.toHaveTextContent(
      'No experiment binding was saved.',
    );
  });

  it('lets the ordinary entry choose a project, samples, and an experiment', async () => {
    generatedMocks.listSamples.mockResolvedValue({
      ok: true,
      project_id: PROJECT,
      samples: [
        {
          sample_id: `smp_${'a'.repeat(32)}`,
          display_name: 'Sample A',
          sample_revision_id: REV_A,
        },
        {
          sample_id: `smp_${'b'.repeat(32)}`,
          display_name: 'Sample B',
          sample_revision_id: REV_B,
        },
      ],
      next_cursor: null,
    });
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    expect(await screen.findByRole('heading', { name: 'Input workbench' })).toBeInTheDocument();
    expect(document.getElementById('experiment-revision-choice')).toBeInTheDocument();
    const project = await screen.findByRole('combobox', { name: 'Project' });
    await waitFor(() => expect(project).toHaveValue(PROJECT));
    await user.click(screen.getByRole('checkbox', { name: 'Sample A' }));
    await user.click(screen.getByRole('checkbox', { name: 'Sample B' }));
    expect(screen.getByLabelText('Submitted sample revision 1')).toHaveValue(REV_A);
    expect(screen.getByLabelText('Submitted sample revision 2')).toHaveValue(REV_B);
    const choice = screen.getByRole('combobox', { name: 'Experiment revision' });
    await waitFor(() => expect(choice).toBeEnabled());
    await user.selectOptions(choice, 'Synthetic experiment');
    expect(screen.getByTestId('experiment-selection-preview')).toHaveTextContent(
      'Order does not matter',
    );
    expect(screen.getByTestId('experiment-binding-status')).toHaveTextContent(
      'No experiment binding. Choosing a revision or opening this page does not create one.',
    );
  });
});
