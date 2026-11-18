import { describe, expect, it, vi } from 'vitest';
import {
  acceptedViewerMessage,
  closeViewerMessage,
  configureViewerMessage,
  framedInstanceId,
  listeningViewerMessage,
  newViewerInstanceId,
  postViewerMessage,
  readHostCommand,
  readViewerNotice,
  standaloneRunId,
  type ViewerConfigurePayload,
} from './contract';

const INSTANCE = '11111111-1111-4111-8111-111111111111';
const OTHER = '22222222-2222-4222-8222-222222222222';
const ORIGIN = 'http://platform.test';

function payload(overrides: Partial<ViewerConfigurePayload> = {}): ViewerConfigurePayload {
  return {
    runId: 'run-1',
    region: 'chr1:500-2700',
    artifactIds: ['artifact-peak'],
    display: { showToolbar: true },
    ...overrides,
  };
}

function channel(source: object) {
  return {
    event: { origin: ORIGIN, source },
    expected: { origin: ORIGIN, source },
  };
}

describe('viewer message contract', () => {
  it('accepts a same-origin configure for the current instance', () => {
    const source = {};
    const message = configureViewerMessage(INSTANCE, payload());
    const { event, expected } = channel(source);
    expect(readHostCommand(message, INSTANCE, event, expected)).toEqual(message);
  });

  it('rejects a different origin, source, version, or instance', () => {
    const source = {};
    const message = configureViewerMessage(INSTANCE, payload());
    expect(
      readHostCommand(message, INSTANCE, { origin: 'https://other.test', source }, {
        origin: ORIGIN,
        source,
      }),
    ).toBeNull();
    expect(
      readHostCommand(message, INSTANCE, { origin: ORIGIN, source: {} }, {
        origin: ORIGIN,
        source,
      }),
    ).toBeNull();
    expect(
      readHostCommand({ ...message, version: 2 }, INSTANCE, channel(source).event, channel(source).expected),
    ).toBeNull();
    expect(readHostCommand(message, OTHER, channel(source).event, channel(source).expected)).toBeNull();
  });

  it('rejects wildcard delivery, urls, and unbounded track lists', () => {
    const target = { postMessage: vi.fn() };
    postViewerMessage(target as unknown as Window, closeViewerMessage(INSTANCE), '*');
    postViewerMessage(target as unknown as Window, closeViewerMessage(INSTANCE), '');
    expect(target.postMessage).not.toHaveBeenCalled();
    postViewerMessage(target as unknown as Window, closeViewerMessage(INSTANCE), ORIGIN);
    expect(target.postMessage).toHaveBeenCalledWith(closeViewerMessage(INSTANCE), ORIGIN);

    const source = {};
    const { event, expected } = channel(source);
    expect(
      readHostCommand(
        configureViewerMessage(INSTANCE, payload({ runId: 'https://evil.test/read' })),
        INSTANCE,
        event,
        expected,
      ),
    ).toBeNull();
    expect(
      readHostCommand(
        configureViewerMessage(INSTANCE, payload({ artifactIds: Array.from({ length: 33 }, (_, index) => `track-${index}`) })),
        INSTANCE,
        event,
        expected,
      ),
    ).toBeNull();
    expect(
      readHostCommand(
        {
          ...configureViewerMessage(INSTANCE, payload()),
          payload: { ...payload(), contentUrl: 'https://evil.test/track' },
        },
        INSTANCE,
        event,
        expected,
      ),
    ).toBeNull();
  });

  it('keeps a late notice from a retired instance from being read', () => {
    const source = {};
    const notice = {
      version: 1 as const,
      type: 'helixweave.viewer.ready' as const,
      instanceId: OTHER,
      payload: {
        runId: 'run-1',
        referenceStatus: 'bound' as const,
        artifactIds: ['artifact-peak'],
      },
    };
    expect(readViewerNotice(notice, INSTANCE, channel(source).event, channel(source).expected)).toBeNull();
    expect(readViewerNotice(notice, OTHER, channel(source).event, channel(source).expected)?.type).toBe(
      'helixweave.viewer.ready',
    );
  });

  it('accepts a listening announcement and a configure receipt on the same channel', () => {
    const source = {};
    const { event, expected } = channel(source);
    const epoch = '33333333-3333-4333-8333-333333333333';
    expect(readViewerNotice(listeningViewerMessage(INSTANCE, epoch), INSTANCE, event, expected)).toEqual(
      listeningViewerMessage(INSTANCE, epoch),
    );
    expect(readViewerNotice(acceptedViewerMessage(INSTANCE, 'run-1'), INSTANCE, event, expected)).toEqual(
      acceptedViewerMessage(INSTANCE, 'run-1'),
    );
    expect(
      readViewerNotice(listeningViewerMessage(INSTANCE, epoch), INSTANCE, { origin: 'https://other.test', source }, expected),
    ).toBeNull();
    expect(
      readViewerNotice(listeningViewerMessage(INSTANCE, epoch), OTHER, event, expected),
    ).toBeNull();
    expect(
      readViewerNotice(
        { ...listeningViewerMessage(INSTANCE, epoch), payload: { runId: 'run-1' } },
        INSTANCE,
        event,
        expected,
      ),
    ).toBeNull();
    expect(
      readViewerNotice(listeningViewerMessage(INSTANCE, 'not-an-epoch'), INSTANCE, event, expected),
    ).toBeNull();
    expect(
      readViewerNotice(
        { ...acceptedViewerMessage(INSTANCE, 'run-1'), payload: { runId: 'https://evil.test/read' } },
        INSTANCE,
        event,
        expected,
      ),
    ).toBeNull();
  });

  it('reads only a bounded run id and instance id from the page address', () => {
    expect(standaloneRunId('?run=run-1&url=https://evil.test')).toBe('run-1');
    expect(standaloneRunId('?run=../secret')).toBeNull();
    expect(framedInstanceId(`#instance=${newViewerInstanceId()}`)).toMatch(
      /^[0-9a-f-]{36}$/i,
    );
    expect(framedInstanceId('#instance=session-cookie')).toBeNull();
  });
});
