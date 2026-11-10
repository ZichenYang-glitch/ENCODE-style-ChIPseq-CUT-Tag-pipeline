import { useState } from 'react';
import { Link, type RouteObject } from 'react-router-dom';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { renderWithRouter } from '../../test/test-utils';
import { ConfigurationDraftWorkbench } from './ConfigurationDraftWorkbench';
import type { ConfigurationDraftWorkbenchProps } from './ConfigurationDraftWorkbench';
import {
  createCapabilityFixture,
  createDraftEditorSchema,
  createDraftSummary,
  createPendingProposal,
} from './testFixtures';

vi.mock('@uiw/react-codemirror', () => ({
  default: () => <textarea aria-label="Advanced config YAML" />,
}));

const schema = createDraftEditorSchema();

function props(
  overrides: Partial<ConfigurationDraftWorkbenchProps> = {},
): ConfigurationDraftWorkbenchProps {
  return {
    workflowName: 'Example workflow',
    projectName: 'Example project',
    schema,
    schemaIdentity: 'example-schema',
    schemaMessage: null,
    list: {
      status: 'ready',
      drafts: [
        createDraftSummary('one', 'Alpha project'),
        createDraftSummary('two', 'Beta project'),
      ],
      pageLimit: 20,
      hasPrevious: false,
      hasNext: false,
    },
    selectedSid: 'one',
    draftExpired: false,
    connection: { state: 'disconnected', draftRestored: false },
    capability: createCapabilityFixture(),
    proposal: createPendingProposal(),
    explanation: {
      contentAvailable: true,
      labels: ['Example needs a sample selection.'],
    },
    validation: { status: 'unavailable' },
    requestStatus: 'idle',
    writesEnabled: true,
    ...overrides,
  };
}

function renderWorkbench(
  overrides: Partial<ConfigurationDraftWorkbenchProps> = {},
  element?: RouteObject['element'],
) {
  return renderWithRouter(
    [
      {
        path: '/drafts',
        element: element ?? <ConfigurationDraftWorkbench {...props(overrides)} />,
      },
      { path: '/elsewhere', element: <p>Elsewhere</p> },
    ],
    { initialEntries: ['/drafts'] },
  );
}

describe('configuration draft workbench', () => {
  it('does not save, adopt, or validate while the editor initializes', async () => {
    const onSave = vi.fn();
    const onAdopt = vi.fn();
    const onValidate = vi.fn();
    renderWorkbench({ onSave, onAdopt, onValidate });
    expect(await screen.findByLabelText('Local limit')).toHaveValue(1);
    expect(screen.getByRole('heading', { name: 'Configuration drafts' })).toBeInTheDocument();
    expect(screen.getByText('Workflow:').parentElement).toHaveTextContent(
      'Example workflow',
    );
    expect(screen.getByText('Project:').parentElement).toHaveTextContent(
      'Example project',
    );
    expect(screen.getByTestId('configuration-draft-layout').className).toContain(
      'flex-col',
    );
    expect(screen.getByTestId('configuration-draft-layout').className).toContain(
      'lg:flex-row',
    );
    await waitFor(() => {
      expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
    });
    expect(onSave).not.toHaveBeenCalled();
    expect(onAdopt).not.toHaveBeenCalled();
    expect(onValidate).not.toHaveBeenCalled();
  });

  it('requires save, discard, or cancel before leaving unsaved edits', async () => {
    const user = userEvent.setup();
    function Page() {
      const [selected, setSelected] = useState('one');
      const [saves, setSaves] = useState(0);
      return (
        <>
          <Link to="/elsewhere">Leave drafts</Link>
          <ConfigurationDraftWorkbench
            {...props({
              selectedSid: selected,
              onSelectDraft: setSelected,
              onSave: () => setSaves((value) => value + 1),
            })}
          />
          <p>Save count {saves}</p>
        </>
      );
    }
    const { router } = renderWorkbench({}, <Page />);
    const limit = await screen.findByLabelText('Local limit');
    await user.clear(limit);
    await user.type(limit, '3');
    expect(await screen.findByText(/you have unsaved edits/i)).toBeInTheDocument();

    await user.click(screen.getByRole('button', { name: /Beta project/ }));
    const dialog = await screen.findByRole('dialog', { name: 'Unsaved edits' });
    expect(within(dialog).getByRole('button', { name: 'Cancel' })).toHaveFocus();
    expect(dialog.className).toContain('max-w-lg');
    await user.keyboard('{Escape}');
    expect(screen.queryByRole('dialog')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Alpha project/ })).toHaveAttribute(
      'aria-pressed',
      'true',
    );

    await user.click(screen.getByRole('button', { name: /Beta project/ }));
    const saveDialog = await screen.findByRole('dialog', { name: 'Unsaved edits' });
    await user.click(within(saveDialog).getByRole('button', { name: 'Save edits' }));
    expect(screen.getByText('Save count 1')).toBeInTheDocument();
    expect(screen.getAllByText(/has not been confirmed/i).length).toBeGreaterThan(0);
    expect(screen.queryByText('Suggestion adopted.')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Alpha project/ })).toHaveAttribute(
      'aria-pressed',
      'true',
    );

    await user.click(within(saveDialog).getByRole('button', { name: 'Discard edits' }));
    await waitFor(() => {
      expect(screen.queryByRole('dialog')).not.toBeInTheDocument();
    });
    expect(screen.getByText('Unsaved edits discarded.')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Beta project/ })).toHaveAttribute(
      'aria-pressed',
      'true',
    );
    expect(screen.getByText('Save count 1')).toBeInTheDocument();

    await user.clear(await screen.findByLabelText('Local limit'));
    await user.type(screen.getByLabelText('Local limit'), '5');
    expect(await screen.findByText(/you have unsaved edits/i)).toBeInTheDocument();
    await user.click(screen.getByRole('link', { name: 'Leave drafts' }));
    const leaveDialog = await screen.findByRole('dialog', { name: 'Unsaved edits' });
    expect(router.state.location.pathname).toBe('/drafts');
    await user.click(within(leaveDialog).getByRole('button', { name: 'Cancel' }));
    expect(router.state.location.pathname).toBe('/drafts');
    await user.click(screen.getByRole('link', { name: 'Leave drafts' }));
    await user.click(
      within(await screen.findByRole('dialog', { name: 'Unsaved edits' })).getByRole(
        'button',
        { name: 'Discard edits' },
      ),
    );
    await waitFor(() => {
      expect(router.state.location.pathname).toBe('/elsewhere');
    });
    expect(screen.getByText('Elsewhere')).toBeInTheDocument();
  });

  it('shows a restored draft without restoring the assistant', () => {
    renderWorkbench({
      writesEnabled: false,
      connection: { state: 'disconnected', draftRestored: true },
    });
    expect(screen.getByText('Assistant not connected')).toBeInTheDocument();
    expect(screen.getByText(/connect it again/i)).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Connect assistant' })).toBeDisabled();
  });

  it('keeps manual editing available when the assistant is unavailable', async () => {
    renderWorkbench({
      writesEnabled: false,
      connection: { state: 'bridge-unavailable', draftRestored: false },
    });
    expect(await screen.findByLabelText('Local limit')).toBeInTheDocument();
    expect(screen.getByText('Assistant unavailable')).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Connect assistant' })).toBeDisabled();
    expect(screen.getByRole('button', { name: 'Validate draft' })).toBeDisabled();
  });

  it('keeps explanation separate from an adoptable suggestion', () => {
    renderWorkbench();
    expect(screen.getByRole('button', { name: 'Adopt suggestion' })).toBeInTheDocument();
    expect(screen.getByRole('heading', { name: 'Explanation' })).toBeInTheDocument();
    const explanation = screen.getByRole('heading', { name: 'Explanation' }).closest('section');
    expect(explanation?.querySelector('button')).toBeNull();
    expect(screen.getByText(/does not change the draft/i)).toBeInTheDocument();
  });
});
