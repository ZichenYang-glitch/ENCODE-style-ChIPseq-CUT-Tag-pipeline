import { describe, it, expect, vi } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';

vi.mock('./api/agentClient', () => ({
  createAgentApiClient: vi.fn().mockReturnValue({
    chat: vi.fn().mockResolvedValue({
      ok: true,
      session_id: null,
      message: 'Mock agent reply.',
      suggestions: [],
      tool_calls: [],
      issues: [],
    }),
  }),
}));

describe('App shell', () => {
  it('renders the workflow platform heading', async () => {
    render(<App />);
    expect(
      await screen.findByRole('heading', { name: /workflow platform/i }),
    ).toBeInTheDocument();
  });

  it('shows the stub workflow catalog', async () => {
    render(<App />);
    expect(
      await screen.findByText(/ENCODE-style ChIP-seq/i),
    ).toBeInTheDocument();
  });

  it('loads schema hints when a workflow is selected', async () => {
    const user = userEvent.setup();
    render(<App />);
    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);
    expect(await screen.findByText(/Config schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Sample schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Options schema/i)).toBeInTheDocument();
  });

  it('renders stub validation issues when Validate is clicked', async () => {
    const user = userEvent.setup();
    render(<App />);
    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);
    const validateButton = await screen.findByTestId('validate-button');
    await user.click(validateButton);
    expect(await screen.findByText(/Sample sheet is invalid/i)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Validate/i })).toBeInTheDocument();
  });

  it('renders a frontend parse error for invalid JSON and does not require a backend', async () => {
    const user = userEvent.setup();
    render(<App />);
    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);

    const configInput = screen.getByLabelText(/Config \(JSON\)/i);
    await user.clear(configInput);
    await user.type(configInput, 'not json');

    const validateButton = screen.getByTestId('validate-button');
    await user.click(validateButton);

    expect(await screen.findByText(/FRONTEND_INPUT_INVALID/i)).toBeInTheDocument();
    expect(screen.getByText(/Invalid JSON for config/i)).toBeInTheDocument();
  });

  it('renders the read-only agent sidebar label', async () => {
    render(<App />);
    expect(
      await screen.findByText(/Validation Assistant — Read Only/i),
    ).toBeInTheDocument();
  });

  it('renders the Run progress panel and disables Create run before validation', async () => {
    render(<App />);

    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await userEvent.click(workflowButton);

    expect(await screen.findByText(/Run progress/i)).toBeInTheDocument();
    expect(
      screen.getByText(/Validate inputs before creating a run record/i),
    ).toBeInTheDocument();
    expect(screen.getByTestId('create-run-button')).toBeDisabled();
  });

  it('enables Create run after successful validation', async () => {
    const user = userEvent.setup();
    render(<App />);

    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);

    const samplesInput = screen.getByLabelText(/Samples \(path string\)/i);
    await user.type(samplesInput, 'samples.tsv');

    const validateButton = await screen.findByTestId('validate-button');
    await user.click(validateButton);

    expect(await screen.findByTestId('create-run-button')).toBeEnabled();
  });

  it('prefills the agent sidebar when Ask Agent is clicked on an issue', async () => {
    const user = userEvent.setup();
    render(<App />);

    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);

    const validateButton = await screen.findByTestId('validate-button');
    await user.click(validateButton);

    const askButton = await screen.findByRole('button', {
      name: /Ask Agent about ENCODE_SAMPLES_INVALID/i,
    });
    await user.click(askButton);

    expect(
      await screen.findByDisplayValue('Explain issue ENCODE_SAMPLES_INVALID.'),
    ).toBeInTheDocument();
  });

  it('sends only the clicked issue in current_issues when Ask Agent is clicked and sent', async () => {
    const { createAgentApiClient } = await import('./api/agentClient');
    const agentClient = vi.mocked(createAgentApiClient).mock.results[0].value as {
      chat: ReturnType<typeof vi.fn>;
    };
    const chat = vi.fn().mockResolvedValue({
      ok: true,
      session_id: null,
      message: 'Mock agent reply.',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });
    agentClient.chat = chat;

    const user = userEvent.setup();
    render(<App />);

    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);

    const validateButton = await screen.findByTestId('validate-button');
    await user.click(validateButton);

    const askButton = await screen.findByRole('button', {
      name: /Ask Agent about ENCODE_SAMPLES_INVALID/i,
    });
    await user.click(askButton);

    expect(
      await screen.findByDisplayValue('Explain issue ENCODE_SAMPLES_INVALID.'),
    ).toBeInTheDocument();

    const sendButton = screen.getByRole('button', { name: /Send message/i });
    await user.click(sendButton);

    await waitFor(() => {
      expect(chat).toHaveBeenCalledTimes(1);
    });

    const [, request] = chat.mock.calls[0];
    expect(request.context.current_issues).toHaveLength(1);
    expect(request.context.current_issues[0]).toMatchObject({
      code: 'ENCODE_SAMPLES_INVALID',
    });
  });
});
