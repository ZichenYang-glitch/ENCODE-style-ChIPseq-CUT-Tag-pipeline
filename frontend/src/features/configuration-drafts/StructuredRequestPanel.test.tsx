import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { StructuredRequestPanel } from './StructuredRequestPanel';
import { createCapabilityFixture } from './testFixtures';

function renderRequest(
  overrides: Partial<Parameters<typeof StructuredRequestPanel>[0]> = {},
) {
  const onRequest = vi.fn();
  render(
    <StructuredRequestPanel
      capability={createCapabilityFixture()}
      connectionAvailable
      draftExpired={false}
      writesEnabled
      dirty={false}
      requestStatus="idle"
      onRequest={onRequest}
      {...overrides}
    />,
  );
  return onRequest;
}

describe('structured suggestion request', () => {
  it('does not submit during initialization and only offers declared fields', async () => {
    const onRequest = renderRequest();
    expect(screen.queryByRole('textbox', { name: /message|chat|prompt|log/i })).not.toBeInTheDocument();
    expect(document.querySelector('textarea, input[type="file"]')).toBeNull();
    expect(screen.queryByText(/threads|cores/i)).not.toBeInTheDocument();
    expect(onRequest).not.toHaveBeenCalled();

    const user = userEvent.setup();
    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    await user.click(screen.getByRole('checkbox', { name: 'Example limit' }));
    await user.type(screen.getByRole('textbox', { name: /Preferred value/ }), '0');
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));

    expect(onRequest).toHaveBeenCalledTimes(1);
    expect(onRequest).toHaveBeenCalledWith({
      intent: 'suggest',
      fieldTokens: ['resource.limit'],
      preferences: [{ fieldToken: 'resource.limit', value: 0 }],
      issueCode: null,
    });
    expect(screen.queryByText(/suggestion adopted|request succeeded/i)).not.toBeInTheDocument();
  });

  it('sends an unchecked boolean preference as false', async () => {
    const onRequest = renderRequest();
    const user = userEvent.setup();
    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    await user.click(screen.getByRole('checkbox', { name: 'Example switch' }));
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(onRequest).toHaveBeenCalledWith({
      intent: 'suggest',
      fieldTokens: ['resource.enabled'],
      preferences: [{ fieldToken: 'resource.enabled', value: false }],
      issueCode: null,
    });
  });

  it('explains stored issues without asking for one to adopt', async () => {
    const onRequest = renderRequest();
    const user = userEvent.setup();
    await user.click(screen.getByRole('radio', { name: 'Explain a validation issue' }));
    expect(screen.queryByRole('checkbox', { name: 'Example limit' })).not.toBeInTheDocument();
    expect(screen.queryByRole('radio', { name: 'Example needs a sample selection.' })).not.toBeInTheDocument();
    expect(screen.getByText('Example needs a sample selection.')).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(onRequest).toHaveBeenCalledWith({
      intent: 'explain_issues',
      fieldTokens: [],
      preferences: [],
      issueCode: null,
    });
    expect(screen.queryByRole('button', { name: 'Adopt suggestion' })).not.toBeInTheDocument();
  });

  it('keeps an unconfirmed request from being sent again', async () => {
    const onRequest = renderRequest({ requestStatus: 'unconfirmed' });
    const user = userEvent.setup();
    expect(screen.getByRole('status')).toHaveTextContent(/will not send it again/i);
    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    expect(screen.getByRole('button', { name: 'Request suggestion' })).toBeDisabled();
    expect(onRequest).not.toHaveBeenCalled();
  });

  it('does not invent fields when the catalog is not connected', () => {
    render(
      <StructuredRequestPanel
        capability={null}
        connectionAvailable={false}
        draftExpired={false}
        writesEnabled={false}
        dirty={false}
        requestStatus="idle"
      />,
    );
    expect(
      screen.getByText(/suggestion fields are not connected yet/i),
    ).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Request suggestion' })).not.toBeInTheDocument();
    expect(document.querySelector('textarea, input[type="file"]')).toBeNull();
  });

  it('keeps the strings omit and invalid as preference values', async () => {
    const onRequest = renderRequest({
      capability: {
        available: true,
        unavailableReason: null,
        intents: ['suggest'],
        fields: [{ token: 'resource.label', label: 'Example label', type: 'string' }],
        issues: [],
      },
    });
    const user = userEvent.setup();
    await user.click(screen.getByRole('radio', { name: 'Suggest values' }));
    await user.click(screen.getByRole('checkbox', { name: 'Example label' }));
    const preference = screen.getByRole('textbox', { name: 'Preferred value' });
    await user.type(preference, 'omit');
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(onRequest).toHaveBeenLastCalledWith({
      intent: 'suggest',
      fieldTokens: ['resource.label'],
      preferences: [{ fieldToken: 'resource.label', value: 'omit' }],
      issueCode: null,
    });
    await user.clear(preference);
    await user.type(preference, 'invalid');
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(onRequest).toHaveBeenLastCalledWith({
      intent: 'suggest',
      fieldTokens: ['resource.label'],
      preferences: [{ fieldToken: 'resource.label', value: 'invalid' }],
      issueCode: null,
    });
  });

  it('clears a request type removed by the current capability', async () => {
    const onRequest = vi.fn();
    const initial = createCapabilityFixture();
    const { rerender } = render(
      <StructuredRequestPanel
        capability={initial}
        connectionAvailable
        draftExpired={false}
        writesEnabled
        dirty={false}
        requestStatus="idle"
        onRequest={onRequest}
      />,
    );
    const user = userEvent.setup();
    await user.click(screen.getByRole('radio', { name: 'Explain a validation issue' }));
    rerender(
      <StructuredRequestPanel
        capability={{ ...initial, intents: ['suggest'], issues: [] }}
        connectionAvailable
        draftExpired={false}
        writesEnabled
        dirty={false}
        requestStatus="idle"
        onRequest={onRequest}
      />,
    );
    await waitFor(() => {
      expect(
        screen.queryByRole('radio', { name: 'Explain a validation issue' }),
      ).not.toBeInTheDocument();
    });
    expect(screen.getByRole('button', { name: 'Request suggestion' })).toBeDisabled();
    await user.click(screen.getByRole('button', { name: 'Request suggestion' }));
    expect(onRequest).not.toHaveBeenCalled();
  });

  it('explains a workflow with no suggestion fields', () => {
    render(
      <StructuredRequestPanel
        capability={{
          available: false,
          unavailableReason: 'Suggestions are not available for this workflow.',
          intents: [],
          fields: [],
          issues: [],
        }}
        connectionAvailable={false}
        draftExpired={false}
        writesEnabled={false}
        dirty={false}
        requestStatus="idle"
      />,
    );
    expect(
      screen.getByText('Suggestions are not available for this workflow.'),
    ).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'Request suggestion' })).not.toBeInTheDocument();
    expect(document.querySelector('textarea, input[type="file"]')).toBeNull();
  });
});
