import { act, fireEvent, render, screen } from '@testing-library/react';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { closeViewerMessage } from '../../viewer/contract';
import { VIEWER_HANDSHAKE_TIMEOUT_MS, ViewerHost } from './ViewerHost';

function currentFrame() {
  const frame = screen.getByTitle('Genome browser') as HTMLIFrameElement;
  const params = new URLSearchParams((frame.getAttribute('src') ?? '').split('#')[1] ?? '');
  const instanceId = params.get('instance') ?? '';
  const target = frame.contentWindow;
  if (target === null) throw new Error('iframe has no content window');
  return { frame, instanceId, target };
}

function spyFrame() {
  const current = currentFrame();
  const postMessage = vi.spyOn(current.target, 'postMessage').mockImplementation(() => {});
  return { ...current, postMessage };
}

function deliver(data: unknown, source: Window | null, origin = window.location.origin) {
  act(() => {
    window.dispatchEvent(new MessageEvent('message', { data, origin, source }));
  });
}

const EPOCH_A = '33333333-3333-4333-8333-333333333333';
const EPOCH_B = '44444444-4444-4444-8444-444444444444';

function listening(instanceId: string, epoch = EPOCH_A) {
  return {
    version: 1,
    type: 'helixweave.viewer.listening',
    instanceId,
    payload: { epoch },
  };
}

function accepted(instanceId: string, runId = 'run-1') {
  return {
    version: 1,
    type: 'helixweave.viewer.accepted',
    instanceId,
    payload: { runId },
  };
}

function ready(instanceId: string, runId = 'run-1') {
  return {
    version: 1,
    type: 'helixweave.viewer.ready',
    instanceId,
    payload: { runId, referenceStatus: 'bound', artifactIds: ['artifact-peak'] },
  };
}

function configureCalls(postMessage: { mock: { calls: unknown[][] } }, instanceId?: string) {
  return postMessage.mock.calls.filter((call) => {
    const message = call[0] as { type?: string; instanceId?: string };
    return message?.type === 'helixweave.viewer.configure' && (instanceId === undefined || message.instanceId === instanceId);
  });
}

describe('ViewerHost', () => {
  afterEach(() => {
    vi.useRealTimers();
    vi.restoreAllMocks();
  });

  it('embeds the same-origin viewer and opens the full page without a wildcard target', () => {
    render(<ViewerHost runId="run-1" />);
    const frame = screen.getByTitle('Genome browser');
    expect(frame).toHaveAttribute('src', expect.stringMatching(/^\/viewer\.html#instance=[0-9a-f-]{36}$/i));
    expect(frame.getAttribute('src')).not.toContain('*');
    const fullPage = screen.getByRole('link', { name: 'Open genome browser full page' });
    expect(fullPage).toHaveAttribute('href', '/viewer.html?run=run-1');
    expect(fullPage).toHaveAttribute('rel', 'noopener noreferrer');
  });

  it('ignores a ready notice from a retired instance', () => {
    render(<ViewerHost runId="run-1" />);
    const { target } = currentFrame();
    const retired = '22222222-2222-4222-8222-222222222222';
    deliver(ready(retired), target);
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
  });

  it('holds configure through a slow document until the child listener announces', () => {
    vi.useFakeTimers();
    render(<ViewerHost runId="run-1" />);
    const { instanceId, target, postMessage } = spyFrame();
    act(() => {
      vi.advanceTimersByTime(3_500);
    });
    expect(configureCalls(postMessage)).toHaveLength(0);
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');

    deliver(listening(instanceId), target);
    const sent = configureCalls(postMessage, instanceId);
    expect(sent).toHaveLength(1);
    expect(sent[0]?.[0]).toEqual(
      expect.objectContaining({
        version: 1,
        type: 'helixweave.viewer.configure',
        instanceId,
        payload: {
          runId: 'run-1',
          region: null,
          artifactIds: null,
          display: { showToolbar: true },
        },
      }),
    );
    expect(sent[0]?.[1]).toBe(window.location.origin);
    act(() => {
      vi.advanceTimersByTime(1_000);
    });
    expect(configureCalls(postMessage, instanceId)).toHaveLength(1);
  });

  it('ignores a repeated announcement and configures again when the frame reloads', () => {
    vi.useFakeTimers();
    render(<ViewerHost runId="run-1" />);
    const { instanceId, target, postMessage } = spyFrame();
    deliver(listening(instanceId), target);
    deliver(accepted(instanceId), target);
    deliver(ready(instanceId), target);
    expect(screen.getByRole('status')).toHaveTextContent(
      'Genome browser is using the authorized track API.',
    );
    act(() => {
      vi.advanceTimersByTime(3_000);
    });
    expect(configureCalls(postMessage, instanceId)).toHaveLength(1);

    deliver(listening(instanceId, EPOCH_A), target);
    expect(configureCalls(postMessage, instanceId)).toHaveLength(1);
    expect(screen.getByRole('status')).toHaveTextContent(
      'Genome browser is using the authorized track API.',
    );

    deliver(listening(instanceId, EPOCH_B), target);
    expect(configureCalls(postMessage, instanceId)).toHaveLength(2);
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
  });

  it('closes a destroyed instance and does not configure it after the run changes', () => {
    const first = render(<ViewerHost runId="run-1" />);
    const old = spyFrame();
    deliver(listening(old.instanceId), old.target);
    expect(configureCalls(old.postMessage, old.instanceId)).toHaveLength(1);

    first.rerender(<ViewerHost runId="run-2" />);
    expect(old.postMessage.mock.calls.some((call) => {
      const message = call[0] as { type?: string; instanceId?: string };
      return message?.type === 'helixweave.viewer.close' && message.instanceId === old.instanceId;
    })).toBe(true);
    expect(old.postMessage).toHaveBeenCalledWith(
      closeViewerMessage(old.instanceId),
      window.location.origin,
    );

    const callsBeforeLate = old.postMessage.mock.calls.length;
    deliver(listening(old.instanceId), old.target);
    expect(old.postMessage.mock.calls.length).toBe(callsBeforeLate);

    const next = spyFrame();
    expect(next.instanceId).not.toBe(old.instanceId);
    next.postMessage.mockClear();
    deliver(listening(old.instanceId), next.target);
    expect(configureCalls(next.postMessage)).toHaveLength(0);
    deliver(listening(next.instanceId), next.target);
    expect(configureCalls(next.postMessage, next.instanceId)).toHaveLength(1);
    expect(configureCalls(next.postMessage, next.instanceId)[0]?.[0]).toEqual(
      expect.objectContaining({
        payload: expect.objectContaining({ runId: 'run-2' }),
      }),
    );

    first.unmount();
    const afterUnmount = next.postMessage.mock.calls.length;
    deliver(listening(next.instanceId), next.target);
    expect(next.postMessage.mock.calls.length).toBe(afterUnmount);
    expect(configureCalls(next.postMessage, next.instanceId)).toHaveLength(1);
  });

  it('shows a bounded retry state when the listener never announces, then opens on retry', () => {
    vi.useFakeTimers();
    render(<ViewerHost runId="run-1" />);
    const first = spyFrame();
    act(() => {
      vi.advanceTimersByTime(VIEWER_HANDSHAKE_TIMEOUT_MS - 1);
    });
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
    expect(screen.queryByRole('button', { name: 'Retry genome browser' })).not.toBeInTheDocument();

    act(() => {
      vi.advanceTimersByTime(1);
    });
    expect(screen.getByRole('status')).toHaveTextContent('Genome browser did not finish opening.');
    fireEvent.click(screen.getByRole('button', { name: 'Retry genome browser' }));
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
    const second = spyFrame();
    expect(second.instanceId).not.toBe(first.instanceId);
    act(() => {
      vi.advanceTimersByTime(3_500);
    });
    expect(configureCalls(second.postMessage)).toHaveLength(0);
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
    deliver(listening(second.instanceId), second.target);
    expect(configureCalls(second.postMessage, second.instanceId)).toHaveLength(1);
  });

  it('ignores a listening message from another origin or source', () => {
    render(<ViewerHost runId="run-1" />);
    const { instanceId, target, postMessage } = spyFrame();
    const other = document.createElement('iframe');
    document.body.appendChild(other);
    deliver(listening(instanceId), target, 'https://evil.test');
    deliver(listening(instanceId), other.contentWindow);
    expect(configureCalls(postMessage)).toHaveLength(0);
    other.remove();
  });

  it('recovers from a handshake timeout once configure is accepted and still waits for ready', () => {
    vi.useFakeTimers();
    render(<ViewerHost runId="run-1" />);
    const { instanceId, target, postMessage } = spyFrame();
    act(() => {
      vi.advanceTimersByTime(VIEWER_HANDSHAKE_TIMEOUT_MS);
    });
    expect(screen.getByRole('status')).toHaveTextContent('Genome browser did not finish opening.');

    deliver(listening(instanceId), target);
    deliver(accepted(instanceId), target);
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
    act(() => {
      vi.advanceTimersByTime(VIEWER_HANDSHAKE_TIMEOUT_MS);
    });
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
    expect(configureCalls(postMessage, instanceId)).toHaveLength(1);

    deliver(ready(instanceId), target);
    expect(screen.getByRole('status')).toHaveTextContent(
      'Genome browser is using the authorized track API.',
    );
  });
});
