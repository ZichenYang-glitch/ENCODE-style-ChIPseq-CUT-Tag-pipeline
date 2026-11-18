import { describe, expect, it, vi } from 'vitest';
import {
  configureViewerMessage,
  listeningViewerMessage,
  type ViewerConfigurePayload,
} from './contract';
import { bindFramedViewer, VIEWER_LISTENING_INTERVAL_MS } from './framedSession';

const INSTANCE = '11111111-1111-4111-8111-111111111111';
const OTHER = '22222222-2222-4222-8222-222222222222';
const ORIGIN = 'http://platform.test';

function payload(overrides: Partial<ViewerConfigurePayload> = {}): ViewerConfigurePayload {
  return {
    runId: 'run-1',
    region: null,
    artifactIds: null,
    display: { showToolbar: true },
    ...overrides,
  };
}

function fakeWindow() {
  const listeners = new Set<(event: MessageEvent) => void>();
  let interval: (() => void) | undefined;
  const page = {
    location: { origin: ORIGIN },
    postMessage: vi.fn(),
    addEventListener: (type: string, listener: (event: MessageEvent) => void) => {
      if (type === 'message') listeners.add(listener);
    },
    removeEventListener: (type: string, listener: (event: MessageEvent) => void) => {
      if (type === 'message') listeners.delete(listener);
    },
    setInterval: (handler: () => void, timeout: number) => {
      expect(timeout).toBe(VIEWER_LISTENING_INTERVAL_MS);
      interval = handler;
      return 1;
    },
    clearInterval: () => {
      interval = undefined;
    },
    dispatch(event: { data: unknown; origin: string; source: object }) {
      for (const listener of [...listeners]) listener(event as MessageEvent);
    },
    tick() {
      interval?.();
    },
    listening: () => listeners.size,
  };
  return page;
}

describe('framed viewer handshake', () => {
  it('announces only after the handler is installed and repeats until configure', () => {
    const page = fakeWindow();
    const host = fakeWindow();
    const onConfigure = vi.fn();
    expect(host.postMessage).not.toHaveBeenCalled();
    const stop = bindFramedViewer({
      page: page as unknown as Window,
      host: host as unknown as Window,
      instanceId: INSTANCE,
      onConfigure,
      onClose: vi.fn(),
    });
    expect(page.listening()).toBe(1);
    expect(host.postMessage).toHaveBeenCalledTimes(1);
    const announced = host.postMessage.mock.calls[0]?.[0] as { payload?: { epoch?: string } };
    expect(announced.payload?.epoch).toMatch(
      /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/i,
    );
    expect(host.postMessage).toHaveBeenCalledWith(
      listeningViewerMessage(INSTANCE, announced.payload?.epoch ?? ''),
      ORIGIN,
    );
    page.tick();
    expect(host.postMessage).toHaveBeenCalledTimes(2);
    expect(host.postMessage.mock.calls[1]?.[0]).toEqual(
      listeningViewerMessage(INSTANCE, announced.payload?.epoch ?? ''),
    );

    page.dispatch({
      data: configureViewerMessage(INSTANCE, payload()),
      origin: ORIGIN,
      source: host,
    });
    expect(onConfigure).toHaveBeenCalledWith(payload());
    expect(host.postMessage).toHaveBeenLastCalledWith(
      expect.objectContaining({
        type: 'helixweave.viewer.accepted',
        instanceId: INSTANCE,
        payload: { runId: 'run-1' },
      }),
      ORIGIN,
    );
    const afterConfigure = host.postMessage.mock.calls.length;
    page.tick();
    expect(host.postMessage).toHaveBeenCalledTimes(afterConfigure);
    stop();
    expect(page.listening()).toBe(0);
  });

  it('keeps announcing when configure fails the origin or instance check', () => {
    const page = fakeWindow();
    const host = fakeWindow();
    const onConfigure = vi.fn();
    bindFramedViewer({
      page: page as unknown as Window,
      host: host as unknown as Window,
      instanceId: INSTANCE,
      onConfigure,
      onClose: vi.fn(),
    });
    page.dispatch({
      data: configureViewerMessage(INSTANCE, payload()),
      origin: 'https://evil.test',
      source: host,
    });
    page.dispatch({
      data: configureViewerMessage(OTHER, payload()),
      origin: ORIGIN,
      source: host,
    });
    expect(onConfigure).not.toHaveBeenCalled();
    const announced = host.postMessage.mock.calls.length;
    page.tick();
    expect(host.postMessage.mock.calls.length).toBe(announced + 1);
  });

  it('stops a destroyed session and ignores it after another instance binds', () => {
    const page = fakeWindow();
    const host = fakeWindow();
    const onConfigure = vi.fn();
    const onClose = vi.fn();
    const stop = bindFramedViewer({
      page: page as unknown as Window,
      host: host as unknown as Window,
      instanceId: INSTANCE,
      onConfigure,
      onClose,
    });
    stop();
    page.tick();
    page.dispatch({
      data: configureViewerMessage(INSTANCE, payload()),
      origin: ORIGIN,
      source: host,
    });
    expect(onConfigure).not.toHaveBeenCalled();

    const nextConfigure = vi.fn();
    bindFramedViewer({
      page: page as unknown as Window,
      host: host as unknown as Window,
      instanceId: OTHER,
      onConfigure: nextConfigure,
      onClose: vi.fn(),
    });
    page.dispatch({
      data: configureViewerMessage(INSTANCE, payload()),
      origin: ORIGIN,
      source: host,
    });
    expect(nextConfigure).not.toHaveBeenCalled();
    page.dispatch({
      data: configureViewerMessage(OTHER, payload({ runId: 'run-2' })),
      origin: ORIGIN,
      source: host,
    });
    expect(nextConfigure).toHaveBeenCalledWith(payload({ runId: 'run-2' }));
  });

  it('stops announcements when the host closes the instance', () => {
    const page = fakeWindow();
    const host = fakeWindow();
    const onClose = vi.fn();
    bindFramedViewer({
      page: page as unknown as Window,
      host: host as unknown as Window,
      instanceId: INSTANCE,
      onConfigure: vi.fn(),
      onClose,
    });
    page.dispatch({
      data: {
        version: 1,
        type: 'helixweave.viewer.close',
        instanceId: INSTANCE,
        payload: {},
      },
      origin: ORIGIN,
      source: host,
    });
    expect(onClose).toHaveBeenCalledTimes(1);
    const announced = host.postMessage.mock.calls.length;
    page.tick();
    expect(host.postMessage).toHaveBeenCalledTimes(announced);
  });

  it('does not announce from a standalone page', () => {
    const page = fakeWindow();
    const stop = bindFramedViewer({
      page: page as unknown as Window,
      host: page as unknown as Window,
      instanceId: INSTANCE,
      onConfigure: vi.fn(),
      onClose: vi.fn(),
    });
    expect(page.postMessage).not.toHaveBeenCalled();
    expect(page.listening()).toBe(0);
    stop();
  });
});
