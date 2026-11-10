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
}));

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: generatedMocks.getWorkflowSchema,
  validateWorkflow: generatedMocks.validateWorkflow,
}));

vi.mock('../../api/generated/runs/runs', () => ({
  createRun: generatedMocks.createRun,
}));

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
    generatedMocks.validateWorkflow.mockReset();
    generatedMocks.createRun.mockReset();
  });

  it('prepares the draft page without calling save, validation, or run APIs', async () => {
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
      await screen.findByRole('heading', { name: 'Configuration drafts' }),
    ).toBeInTheDocument();
    expect(
      screen.getByText('ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase'),
    ).toBeInTheDocument();
    expect(screen.getByText('Project list is not connected yet.')).toBeInTheDocument();
    expect(screen.getByText('Assistant not connected')).toBeInTheDocument();
    expect(screen.getByText(/draft list is not connected yet/i)).toBeInTheDocument();
    expect(screen.getAllByText(/does not start a run/i).length).toBeGreaterThan(0);
    expect(screen.getByRole('button', { name: 'Connect assistant' })).toBeDisabled();
    expect(screen.getByRole('button', { name: 'Validate draft' })).toBeDisabled();
    expect(screen.queryByRole('button', { name: 'Adopt suggestion' })).not.toBeInTheDocument();
    expect(screen.queryByRole('link', { name: /create run|preflight|start run/i })).not.toBeInTheDocument();
    expect(
      screen.getByRole('heading', { name: 'Suggestion request' }).closest('section')?.querySelector(
        'textarea, input[type="file"]',
      ),
    ).toBeNull();

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
    });
    expect(getWorkflowSchema).toHaveBeenCalledWith(WORKFLOW_ID);
  });
});
