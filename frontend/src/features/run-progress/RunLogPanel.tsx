import type { RunLogChunkResponse } from '../../api/runTypes';

interface RunLogPanelProps {
  chunks: RunLogChunkResponse[];
}

export function RunLogPanel({ chunks }: RunLogPanelProps) {
  if (chunks.length === 0) {
    return (
      <p
        className="text-sm text-[var(--color-text-muted)]"
        data-testid="run-log-panel-empty"
      >
        No log entries yet.
      </p>
    );
  }

  return (
    <div className="space-y-2" data-testid="run-log-panel">
      {chunks.map((chunk) => (
        <div
          key={chunk.chunk_id}
          className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-2"
        >
          <div className="mb-1 text-xs text-[var(--color-text-muted)]">
            {new Date(chunk.timestamp).toLocaleString()} — {chunk.stream_name}
          </div>
          <pre className="overflow-auto font-mono text-xs text-[var(--color-text)]">
            {chunk.lines.join('\n')}
          </pre>
        </div>
      ))}
    </div>
  );
}
