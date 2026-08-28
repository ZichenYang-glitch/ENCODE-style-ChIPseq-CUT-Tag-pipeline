import { useId, type KeyboardEvent } from 'react';
import type { RunLogChunkResponse } from '../../api/runTypes';

interface RunLogPanelProps {
  stdoutChunks: RunLogChunkResponse[];
  stderrChunks: RunLogChunkResponse[];
  activeStream: 'stdout' | 'stderr';
  onStreamChange: (stream: 'stdout' | 'stderr') => void;
}

export function RunLogPanel({
  stdoutChunks,
  stderrChunks,
  activeStream,
  onStreamChange,
}: RunLogPanelProps) {
  const baseId = useId();
  const stdoutTabId = `${baseId}-stdout-tab`;
  const stderrTabId = `${baseId}-stderr-tab`;
  const panelId = `${baseId}-panel`;
  const chunks = activeStream === 'stdout' ? stdoutChunks : stderrChunks;

  function handleTabKey(event: KeyboardEvent<HTMLButtonElement>) {
    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    const next =
      event.key === 'ArrowLeft' || event.key === 'Home' ? 'stdout' : 'stderr';
    onStreamChange(next);
    document.getElementById(next === 'stdout' ? stdoutTabId : stderrTabId)?.focus();
  }

  return (
    <div className="space-y-2" data-testid="run-log-panel">
      <div className="flex gap-1" role="tablist" aria-label="Log streams">
        <button
          type="button"
          id={stdoutTabId}
          role="tab"
          aria-selected={activeStream === 'stdout'}
          aria-controls={panelId}
          tabIndex={activeStream === 'stdout' ? 0 : -1}
          className={`min-h-9 rounded-[4px] px-3 py-1 text-xs font-medium ${
            activeStream === 'stdout'
              ? 'bg-[var(--color-accent)] text-white'
              : 'border border-[var(--color-border)] bg-[var(--color-bg)] text-[var(--color-text)]'
          }`}
          onClick={() => onStreamChange('stdout')}
          onKeyDown={handleTabKey}
          data-testid="stdout-tab"
        >
          stdout
        </button>
        <button
          type="button"
          id={stderrTabId}
          role="tab"
          aria-selected={activeStream === 'stderr'}
          aria-controls={panelId}
          tabIndex={activeStream === 'stderr' ? 0 : -1}
          className={`min-h-9 rounded-[4px] px-3 py-1 text-xs font-medium ${
            activeStream === 'stderr'
              ? 'bg-[var(--color-accent)] text-white'
              : 'border border-[var(--color-border)] bg-[var(--color-bg)] text-[var(--color-text)]'
          }`}
          onClick={() => onStreamChange('stderr')}
          onKeyDown={handleTabKey}
          data-testid="stderr-tab"
        >
          stderr
        </button>
      </div>

      <div
        id={panelId}
        role="tabpanel"
        aria-labelledby={activeStream === 'stdout' ? stdoutTabId : stderrTabId}
      >
        {chunks.length === 0 ? (
          <p
            className="text-sm text-[var(--color-text-muted)]"
            data-testid="run-log-panel-empty"
          >
            No log entries yet.
          </p>
        ) : (
          <div
            className="max-h-[28rem] min-w-0 overflow-auto rounded-[4px] border border-[var(--color-border)] bg-[var(--color-surface-subtle)]"
            data-testid="run-log-chunks"
          >
            {chunks.map((chunk) => (
              <div
                key={chunk.chunk_id}
                className="border-b border-[var(--color-border)] p-3 last:border-b-0"
              >
                <div className="mb-1 text-xs text-[var(--color-text-faint)] tabular-nums">
                  {new Date(chunk.timestamp).toLocaleString()} — {chunk.stream_name}
                </div>
                <pre className="max-w-full whitespace-pre-wrap break-words font-mono text-xs leading-[1.125rem] text-[var(--color-text)]">
                  {chunk.lines.join('\n')}
                </pre>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
