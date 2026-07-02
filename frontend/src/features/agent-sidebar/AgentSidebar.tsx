import { useEffect, useState } from 'react';
import type { AgentApiClient } from '../../api/agentClient';
import type {
  AgentRequest,
  AgentResponse,
  AgentSuggestion,
  AgentToolCall,
  Issue,
} from '../../api/types';
import { Button } from '../../components/Button';

interface Message {
  role: 'user' | 'assistant';
  content: string;
  tool_calls?: AgentToolCall[];
  suggestions?: AgentSuggestion[];
}

interface AgentSidebarProps {
  workflowId: string | null;
  issues: Issue[];
  agentClient: AgentApiClient;
  draftMessage?: string;
  onDraftConsumed?: () => void;
  focusedIssue?: Issue | null;
  onFocusedIssueConsumed?: () => void;
}

export function AgentSidebar({
  workflowId,
  issues,
  agentClient,
  draftMessage,
  onDraftConsumed,
  focusedIssue,
  onFocusedIssueConsumed,
}: AgentSidebarProps) {
  const [input, setInput] = useState(draftMessage ?? '');
  const [messages, setMessages] = useState<Message[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [pendingIssue, setPendingIssue] = useState<Issue | null>(null);

  // Prefill input and capture the focused issue when a new draft arrives from
  // the issue panel.
  useEffect(() => {
    if (draftMessage) {
      setInput(draftMessage);
      onDraftConsumed?.();
    }
    if (focusedIssue) {
      setPendingIssue(focusedIssue);
      onFocusedIssueConsumed?.();
    }
  }, [draftMessage, focusedIssue, onDraftConsumed, onFocusedIssueConsumed]);

  async function handleSend() {
    if (!workflowId || input.trim().length === 0) return;

    const userMessage: Message = { role: 'user', content: input.trim() };
    setMessages((prev) => [...prev, userMessage]);
    setInput('');
    setLoading(true);
    setError(null);

    const currentIssues = pendingIssue ? [pendingIssue] : issues;

    const request: AgentRequest = {
      session_id: null,
      message: userMessage.content,
      context: {
        current_issues: currentIssues,
        current_config: {},
        current_schema: {},
      },
    };

    try {
      const response: AgentResponse = await agentClient.chat(
        workflowId,
        request,
      );
      setMessages((prev) => [
        ...prev,
        {
          role: 'assistant',
          content: response.message,
          tool_calls: response.tool_calls,
          suggestions: response.suggestions,
        },
      ]);
    } catch (err) {
      const message =
        err instanceof Error ? err.message : 'Agent request failed.';
      setError(message);
    } finally {
      setLoading(false);
      setPendingIssue(null);
    }
  }

  function handleKeyDown(event: React.KeyboardEvent<HTMLTextAreaElement>) {
    if (event.key === 'Enter' && !event.shiftKey) {
      event.preventDefault();
      void handleSend();
    }
  }

  return (
    <div className="flex h-full flex-col rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] p-3 shadow-sm">
      <header className="mb-3 border-b border-[var(--color-border)] pb-2">
        <h3 className="text-sm font-semibold text-[var(--color-text)]">
          Validation Assistant — Read Only
        </h3>
        <p className="mt-1 text-xs text-[var(--color-text-muted)]">
          I can help you understand validation errors and schema requirements. I
          cannot run, modify, or delete workflows.
        </p>
      </header>

      <div className="flex-1 space-y-3 overflow-y-auto text-sm">
        {messages.length === 0 && (
          <p className="text-xs text-[var(--color-text-muted)]">
            Ask a question about the selected workflow.
          </p>
        )}

        {messages.map((message, index) => (
          <div
            key={index}
            className={`space-y-2 ${
              message.role === 'user' ? 'text-right' : 'text-left'
            }`}
          >
            <div
              className={`inline-block max-w-full rounded border px-2.5 py-2 text-left ${
                message.role === 'user'
                  ? 'border-[var(--color-accent)] bg-[var(--color-accent-bg)] text-[var(--color-text)]'
                  : 'border-[var(--color-border)] bg-[var(--color-bg)] text-[var(--color-text)]'
              }`}
            >
              <p className="text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
                {message.role === 'user' ? 'You' : 'Assistant'}
              </p>
              <p className="mt-1 whitespace-pre-wrap">{message.content}</p>

              {message.tool_calls && message.tool_calls.length > 0 && (
                <div className="mt-2 space-y-2">
                  {message.tool_calls.map((toolCall, toolIndex) => (
                    <ToolCallCard key={toolIndex} toolCall={toolCall} />
                  ))}
                </div>
              )}

              {message.suggestions && message.suggestions.length > 0 && (
                <div className="mt-2 space-y-2">
                  {message.suggestions.map((suggestion, suggestionIndex) => (
                    <SuggestionCard
                      key={suggestionIndex}
                      suggestion={suggestion}
                    />
                  ))}
                </div>
              )}
            </div>
          </div>
        ))}

        {loading && (
          <div className="flex items-center gap-2 text-xs text-[var(--color-text-muted)]">
            <span
              className="inline-block h-4 w-4 animate-spin rounded-full border-2 border-[var(--color-border)] border-t-[var(--color-accent)]"
              role="status"
              aria-label="Loading"
            />
            Assistant is thinking…
          </div>
        )}

        {error && (
          <div
            className="rounded border border-[var(--color-error)] bg-[var(--color-error-bg)] p-2 text-xs text-[var(--color-error)]"
            role="alert"
          >
            {error}
          </div>
        )}
      </div>

      <div className="mt-3 space-y-2 border-t border-[var(--color-border)] pt-3">
        <textarea
          value={input}
          onChange={(event) => setInput(event.target.value)}
          onKeyDown={handleKeyDown}
          placeholder={
            workflowId
              ? 'Ask about this workflow…'
              : 'Select a workflow to ask the assistant.'
          }
          disabled={!workflowId || loading}
          className="min-h-[4.5rem] w-full resize-y rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2.5 py-2 text-sm text-[var(--color-text)] placeholder:text-[var(--color-text-muted)] focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)] disabled:opacity-60"
          aria-label="Agent message input"
        />
        <Button
          variant="primary"
          onClick={handleSend}
          disabled={!workflowId || input.trim().length === 0 || loading}
          className="w-full"
          aria-label="Send message"
        >
          Send
        </Button>
      </div>
    </div>
  );
}

function ToolCallCard({ toolCall }: { toolCall: AgentToolCall }) {
  return (
    <div className="rounded border border-[var(--color-border)] bg-[var(--color-surface)] p-2 text-left text-xs">
      <div className="flex flex-wrap items-center gap-2">
        <span className="font-semibold text-[var(--color-text)]">
          {toolCall.tool_name}
        </span>
        <span className="rounded bg-[var(--color-info-bg)] px-1.5 py-0.5 text-[var(--color-info)]">
          read_only
        </span>
        {toolCall.read_only === true && (
          <span className="text-[var(--color-text-muted)]">true</span>
        )}
      </div>
      {Object.keys(toolCall.input_summary).length > 0 && (
        <details className="mt-1">
          <summary className="cursor-pointer text-[var(--color-accent)] hover:underline">
            Input summary
          </summary>
          <pre className="mt-1 overflow-auto rounded bg-[var(--color-bg)] p-1.5 text-[var(--color-text)]">
            {JSON.stringify(toolCall.input_summary, null, 2)}
          </pre>
        </details>
      )}
      {toolCall.output_summary && (
        <details className="mt-1">
          <summary className="cursor-pointer text-[var(--color-accent)] hover:underline">
            Output summary
          </summary>
          <pre className="mt-1 overflow-auto rounded bg-[var(--color-bg)] p-1.5 text-[var(--color-text)]">
            {toolCall.output_summary}
          </pre>
        </details>
      )}
    </div>
  );
}

function SuggestionCard({ suggestion }: { suggestion: AgentSuggestion }) {
  async function handleCopy() {
    await navigator.clipboard.writeText(suggestion.description);
  }

  return (
    <div className="rounded border border-[var(--color-border)] bg-[var(--color-surface)] p-2 text-left text-xs">
      <p className="font-semibold text-[var(--color-text)]">
        {suggestion.description}
      </p>
      {suggestion.rationale && (
        <p className="mt-1 text-[var(--color-text-muted)]">
          {suggestion.rationale}
        </p>
      )}
      {suggestion.disclaimer && (
        <p className="mt-1 italic text-[var(--color-text-muted)]">
          {suggestion.disclaimer}
        </p>
      )}
      <div className="mt-2 flex flex-wrap gap-2">
        <Button
          variant="secondary"
          onClick={() => void handleCopy()}
          className="text-xs"
          aria-label="Copy suggestion description"
        >
          Copy
        </Button>
      </div>
    </div>
  );
}
