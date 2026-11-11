import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { createRun } from '../../api/generated/runs/runs';
import {
  getWorkflowSchema,
  validateWorkflow,
} from '../../api/generated/workflows/workflows';
import { appRoutes } from '../../app/router';
import { WORKFLOW_ID } from '../../features/input-workbench/test-fixtures';
import { createDraftEditorContract } from '../../features/configuration-drafts/testFixtures';
import { renderWithRouter } from '../../test/test-utils';

const generatedMocks = vi.hoisted(() => ({
  createRun: vi.fn(),
  getWorkflowSchema: vi.fn(),
  validateWorkflow: vi.fn(),
  listCompatibleReferenceProfiles: vi.fn(),
  openConnection: vi.fn(),
  createSession: vi.fn(),
  saveDraft: vi.fn(),
  validateDraft: vi.fn(),
  listProjects: vi.fn(),
  listSessions: vi.fn(),
  listSamples: vi.fn(),
}));

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: generatedMocks.getWorkflowSchema,
  validateWorkflow: generatedMocks.validateWorkflow,
  listCompatibleReferenceProfiles: generatedMocks.listCompatibleReferenceProfiles,
}));

vi.mock('../../api/generated/runs/runs', () => ({
  createRun: generatedMocks.createRun,
}));

vi.mock('../../api/generated/authoring/authoring', () => ({
  authoringListProjects: generatedMocks.listProjects,
  authoringListProjectSamples: generatedMocks.listSamples,
  authoringListSessions: generatedMocks.listSessions,
  authoringCreateSession: generatedMocks.createSession,
  authoringGetSession: vi.fn(),
  authoringSaveDraft: generatedMocks.saveDraft,
  authoringOpenConnection: generatedMocks.openConnection,
  authoringCloseConnection: vi.fn(),
  authoringCreateTurn: vi.fn(),
  authoringGetTurn: vi.fn(),
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

beforeAll(async () => {
  await import('../workflows/configuration-drafts');
});

describe('configuration drafts route', () => {
  beforeEach(() => {
    generatedMocks.getWorkflowSchema.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      schema: createDraftEditorContract(),
      issues: [],
    });
    generatedMocks.listCompatibleReferenceProfiles.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      profiles: [],
      issues: [],
    });
    generatedMocks.listProjects.mockResolvedValue({ ok: true, projects: [] });
    generatedMocks.listSessions.mockResolvedValue({ ok: true, sessions: [] });
    generatedMocks.listSamples.mockResolvedValue({ ok: true, samples: [] });
    generatedMocks.validateWorkflow.mockReset();
    generatedMocks.createRun.mockReset();
    generatedMocks.openConnection.mockReset();
    generatedMocks.createSession.mockReset();
    generatedMocks.saveDraft.mockReset();
    generatedMocks.validateDraft.mockReset();
  });

  it('opens the wired draft page without connecting or starting a run', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });

    const entry = await screen.findAllByRole('link', {
      name: 'Configuration drafts',
    });
    await user.click(entry[0]);
    expect(router.state.location.pathname).toBe(
      `/workflows/${WORKFLOW_ID}/configuration-drafts`,
    );
    expect(
      await screen.findByRole('heading', { level: 1, name: 'Configuration drafts' }),
    ).toBeInTheDocument();
    expect(
      screen.getByText('ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase'),
    ).toBeInTheDocument();
    expect(screen.getByRole('combobox', { name: 'Project' })).toBeInTheDocument();
    expect(screen.getByText(/no projects are available/i)).toBeInTheDocument();
    expect(screen.queryByText(/draft list is not connected yet/i)).not.toBeInTheDocument();
    expect(screen.queryByText(/reviewed configuration-draft interface/i)).not.toBeInTheDocument();
    expect(screen.getByText('Assistant not connected')).toBeInTheDocument();
    expect(screen.getAllByText(/does not start a run/i).length).toBeGreaterThan(0);
    expect(screen.getByRole('button', { name: 'Connect assistant' })).toBeDisabled();
    expect(screen.getByRole('button', { name: 'Validate draft' })).toBeDisabled();
    expect(screen.getByRole('button', { name: 'New draft' })).toBeDisabled();
    expect(screen.queryByRole('button', { name: 'Adopt suggestion' })).not.toBeInTheDocument();
    expect(screen.queryByRole('link', { name: /create run|preflight|start run/i })).not.toBeInTheDocument();
    expect(document.querySelector('input[type="file"]')).toBeNull();

    const limit = await screen.findByLabelText('Local limit');
    expect(limit).toHaveValue(1);
    await user.clear(limit);
    await user.type(limit, '2');
    expect(await screen.findByText(/you have unsaved edits/i)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Save edits' })).toBeDisabled();
    expect(screen.queryByText(/edits saved|suggestion adopted|validation finished/i)).not.toBeInTheDocument();
    await waitFor(() => {
      expect(validateWorkflow).not.toHaveBeenCalled();
      expect(createRun).not.toHaveBeenCalled();
      expect(generatedMocks.openConnection).not.toHaveBeenCalled();
      expect(generatedMocks.createSession).not.toHaveBeenCalled();
      expect(generatedMocks.saveDraft).not.toHaveBeenCalled();
      expect(generatedMocks.validateDraft).not.toHaveBeenCalled();
    });
    expect(getWorkflowSchema).toHaveBeenCalledWith(WORKFLOW_ID);
  });
});
