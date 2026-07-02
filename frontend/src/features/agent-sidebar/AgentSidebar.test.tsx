import { describe, it, expect, vi } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { AgentSidebar } from './AgentSidebar';
import type { AgentResponse, Issue } from '../../api/types';
import type { AgentApiClient } from '../../api/agentClient';

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

const sampleIssue: Issue = {
  code: 'ENCODE_SAMPLES_INVALID',
  message: 'Sample sheet is invalid.',
  severity: 'error',
  path: 'samples',
  source: 'samples',
  technical_message: 'Missing samples field.',
  hint: 'Provide a path to a TSV sample sheet.',
  context: {},
};

function createMockClient(
  response: AgentResponse,
): { client: AgentApiClient; chat: ReturnType<typeof vi.fn> } {
  const chat = vi.fn().mockResolvedValue(response);
  return {
    client: { chat },
    chat,
  };
}

describe('AgentSidebar', () => {
  it('renders the read-only label and disclaimer', () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: '',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[]}
        agentClient={client}
      />,
    );

    expect(
      screen.getByText(/Validation Assistant — Read Only/i),
    ).toBeInTheDocument();
    expect(
      screen.getByText(/I cannot run, modify, or delete workflows/i),
    ).toBeInTheDocument();
  });

  it('sends a message using agentClient.chat with current_issues in context', async () => {
    const { client, chat } = createMockClient({
      ok: true,
      session_id: null,
      message: 'Deterministic mock explanation.',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    const user = userEvent.setup();
    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[sampleIssue]}
        agentClient={client}
      />,
    );

    const input = screen.getByLabelText(/Agent message input/i);
    await user.type(input, 'Explain the issues.');

    const sendButton = screen.getByRole('button', { name: /Send message/i });
    await user.click(sendButton);

    await waitFor(() => {
      expect(chat).toHaveBeenCalledTimes(1);
    });

    expect(chat).toHaveBeenCalledWith(
      WORKFLOW_ID,
      expect.objectContaining({
        session_id: null,
        message: 'Explain the issues.',
        context: expect.objectContaining({
          current_issues: [sampleIssue],
          current_config: {},
          current_schema: {},
        }),
      }),
    );
  });

  it('sends only the focused issue when one is provided', async () => {
    const focusedIssue: Issue = {
      code: 'FRONTEND_INPUT_INVALID',
      message: 'Invalid JSON for config.',
      severity: 'error',
      path: 'config',
      source: 'ui',
      technical_message: 'Unexpected token.',
      hint: 'Use valid JSON.',
      context: {},
    };

    const { client, chat } = createMockClient({
      ok: true,
      session_id: null,
      message: 'Focused issue explanation.',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    const onFocusedIssueConsumed = vi.fn();
    const user = userEvent.setup();
    const { rerender } = render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[sampleIssue]}
        agentClient={client}
        draftMessage="Explain issue FRONTEND_INPUT_INVALID."
        focusedIssue={focusedIssue}
        onFocusedIssueConsumed={onFocusedIssueConsumed}
      />,
    );

    expect(
      await screen.findByDisplayValue('Explain issue FRONTEND_INPUT_INVALID.'),
    ).toBeInTheDocument();
    expect(onFocusedIssueConsumed).toHaveBeenCalled();

    rerender(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[sampleIssue]}
        agentClient={client}
        draftMessage=""
        focusedIssue={null}
        onFocusedIssueConsumed={onFocusedIssueConsumed}
      />,
    );

    const sendButton = screen.getByRole('button', { name: /Send message/i });
    await user.click(sendButton);

    await waitFor(() => {
      expect(chat).toHaveBeenCalledTimes(1);
    });

    expect(chat).toHaveBeenCalledWith(
      WORKFLOW_ID,
      expect.objectContaining({
        session_id: null,
        message: 'Explain issue FRONTEND_INPUT_INVALID.',
        context: expect.objectContaining({
          current_issues: [focusedIssue],
          current_config: {},
          current_schema: {},
        }),
      }),
    );
  });

  it('renders assistant response messages', async () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: 'This is the assistant reply.',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    const user = userEvent.setup();
    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[]}
        agentClient={client}
      />,
    );

    const input = screen.getByLabelText(/Agent message input/i);
    await user.type(input, 'Hello');
    await user.click(screen.getByRole('button', { name: /Send message/i }));

    expect(
      await screen.findByText('This is the assistant reply.'),
    ).toBeInTheDocument();
  });

  it('renders tool calls with read_only badge', async () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: 'Here is the explanation.',
      suggestions: [],
      tool_calls: [
        {
          tool_name: 'explain_issues',
          input_summary: { issues: [sampleIssue] },
          output_summary: 'Explanation text',
          read_only: true,
        },
      ],
      issues: [],
    });

    const user = userEvent.setup();
    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[sampleIssue]}
        agentClient={client}
      />,
    );

    const input = screen.getByLabelText(/Agent message input/i);
    await user.type(input, 'Explain');
    await user.click(screen.getByRole('button', { name: /Send message/i }));

    expect(await screen.findByText('explain_issues')).toBeInTheDocument();
    expect(screen.getByText('read_only')).toBeInTheDocument();
    expect(screen.getByText('true')).toBeInTheDocument();
  });

  it('renders suggestion cards without Apply/Save/Submit/Run buttons', async () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: 'Consider this suggestion.',
      suggestions: [
        {
          type: 'config_edit',
          description: 'Add the samples field.',
          target_path: 'samples',
          current_value: null,
          proposed_value: 'samples.tsv',
          rationale: 'The workflow requires a sample sheet.',
          disclaimer: 'Verify against ENCODE guidelines.',
        },
      ],
      tool_calls: [],
      issues: [],
    });

    const user = userEvent.setup();
    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[]}
        agentClient={client}
      />,
    );

    const input = screen.getByLabelText(/Agent message input/i);
    await user.type(input, 'Fix');
    await user.click(screen.getByRole('button', { name: /Send message/i }));

    expect(
      await screen.findByText('Add the samples field.'),
    ).toBeInTheDocument();
    expect(screen.getByText(/Verify against ENCODE guidelines/i)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Copy/i })).toBeInTheDocument();

    const forbiddenButtons = ['Apply', 'Save', 'Submit', 'Run'];
    for (const label of forbiddenButtons) {
      expect(
        screen.queryByRole('button', { name: new RegExp(label, 'i') }),
      ).not.toBeInTheDocument();
    }
  });

  it('shows an error banner when chat rejects', async () => {
    const chat = vi.fn().mockRejectedValue(new Error('Network error'));
    const client: AgentApiClient = { chat };

    const user = userEvent.setup();
    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[]}
        agentClient={client}
      />,
    );

    const input = screen.getByLabelText(/Agent message input/i);
    await user.type(input, 'Explain');
    await user.click(screen.getByRole('button', { name: /Send message/i }));

    expect(await screen.findByRole('alert')).toHaveTextContent('Network error');
  });

  it('disables input and send when no workflow is selected', async () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: '',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    render(
      <AgentSidebar workflowId={null} issues={[]} agentClient={client} />,
    );

    expect(screen.getByLabelText(/Agent message input/i)).toBeDisabled();
    expect(screen.getByRole('button', { name: /Send message/i })).toBeDisabled();
  });

  it('disables send when input is empty', async () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: '',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    render(
      <AgentSidebar workflowId={WORKFLOW_ID} issues={[]} agentClient={client} />,
    );

    expect(screen.getByRole('button', { name: /Send message/i })).toBeDisabled();
  });

  it('prefills input when draftMessage is provided', async () => {
    const { client } = createMockClient({
      ok: true,
      session_id: null,
      message: '',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });

    const onDraftConsumed = vi.fn();
    render(
      <AgentSidebar
        workflowId={WORKFLOW_ID}
        issues={[]}
        agentClient={client}
        draftMessage="Explain issue ENCODE_SAMPLES_INVALID."
        onDraftConsumed={onDraftConsumed}
      />,
    );

    expect(
      await screen.findByDisplayValue('Explain issue ENCODE_SAMPLES_INVALID.'),
    ).toBeInTheDocument();
    expect(onDraftConsumed).toHaveBeenCalled();
  });
});
