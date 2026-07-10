import { describe, it, expect, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { RunLogPanel } from './RunLogPanel';
import type { RunLogChunkResponse } from '../../api/runTypes';

function makeChunk(overrides: Partial<RunLogChunkResponse> = {}): RunLogChunkResponse {
  return {
    chunk_id: 'log-1',
    run_id: 'run-1',
    stream_name: 'stdout',
    sequence: 1,
    timestamp: '2026-07-04T12:00:00.000Z',
    lines: ['line 1'],
    ...overrides,
  };
}

describe('RunLogPanel', () => {
  it('renders stdout chunks by default', () => {
    render(
      <RunLogPanel
        stdoutChunks={[makeChunk({ lines: ['stdout line'] })]}
        stderrChunks={[]}
        activeStream="stdout"
        onStreamChange={vi.fn()}
      />,
    );

    expect(screen.getByText('stdout line')).toBeInTheDocument();
    expect(screen.queryByText('stderr line')).not.toBeInTheDocument();
  });

  it('switches to stderr chunks when stderr tab is selected', async () => {
    const user = userEvent.setup();
    const onStreamChange = vi.fn();

    render(
      <RunLogPanel
        stdoutChunks={[makeChunk({ lines: ['stdout line'] })]}
        stderrChunks={[makeChunk({ stream_name: 'stderr', lines: ['stderr line'] })]}
        activeStream="stdout"
        onStreamChange={onStreamChange}
      />,
    );

    await user.click(screen.getByTestId('stderr-tab'));
    expect(onStreamChange).toHaveBeenCalledWith('stderr');
  });

  it('renders empty state when active stream has no chunks', () => {
    render(
      <RunLogPanel
        stdoutChunks={[]}
        stderrChunks={[]}
        activeStream="stdout"
        onStreamChange={vi.fn()}
      />,
    );

    expect(screen.getByTestId('run-log-panel-empty')).toHaveTextContent(
      /No log entries yet/i,
    );
  });

  it('does not duplicate rendered chunks when given identical arrays', () => {
    const { rerender } = render(
      <RunLogPanel
        stdoutChunks={[makeChunk({ chunk_id: 'log-1', lines: ['line'] })]}
        stderrChunks={[]}
        activeStream="stdout"
        onStreamChange={vi.fn()}
      />,
    );

    expect(screen.getAllByText('line')).toHaveLength(1);

    rerender(
      <RunLogPanel
        stdoutChunks={[makeChunk({ chunk_id: 'log-1', lines: ['line'] })]}
        stderrChunks={[]}
        activeStream="stdout"
        onStreamChange={vi.fn()}
      />,
    );

    expect(screen.getAllByText('line')).toHaveLength(1);
  });
});
