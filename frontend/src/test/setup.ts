import '@testing-library/jest-dom/vitest';
import { vi } from 'vitest';

// Node's Request accepts only its own AbortSignal. jsdom's signal is a
// different class, so React Router navigations throw before the location
// changes. Forward the DOM signal onto the controller captured before jsdom
// replaced it. Dropping the signal would let a navigation ignore cancellation.
const nodeAbort = (
  globalThis as typeof globalThis & {
    __nodeAbort?: {
      AbortController: typeof AbortController;
      AbortSignal: typeof AbortSignal;
    };
  }
).__nodeAbort;

const OriginalRequest = globalThis.Request;

if (nodeAbort !== undefined) {
  const NodeAbortController = nodeAbort.AbortController;
  const NodeAbortSignal = nodeAbort.AbortSignal;

  function isNodeAbortSignal(signal: AbortSignal): boolean {
    return signal instanceof NodeAbortSignal;
  }

  function nodeAbortSignal(signal: AbortSignal): AbortSignal {
    if (isNodeAbortSignal(signal)) return signal;
    const controller = new NodeAbortController();
    const abort = () => {
      controller.abort(signal.reason);
    };
    if (signal.aborted) abort();
    else signal.addEventListener('abort', abort, { once: true });
    return controller.signal;
  }

  globalThis.Request = class Request extends OriginalRequest {
    constructor(input: RequestInfo | URL, init?: RequestInit) {
      const signal = init?.signal;
      if (signal != null && !(signal instanceof NodeAbortSignal)) {
        super(input, { ...init, signal: nodeAbortSignal(signal) });
        return;
      }
      super(input, init);
    }
  };
}

vi.mock('../api/authClient', async (importOriginal) => {
  const actual =
    await importOriginal<typeof import('../api/authClient')>();
  return {
    ...actual,
    fetchSessionState: vi.fn(async () => ({
      ok: true,
      setup_required: false,
      authenticated: true,
      principal: {
        user_id: 'usr_00000000000000000000000000000000',
        username: 'test-admin',
        role: 'administrator',
      },
    })),
    login: vi.fn(async () => ({
      user_id: 'usr_00000000000000000000000000000000',
      username: 'test-admin',
      role: 'administrator',
    })),
    logout: vi.fn(async () => undefined),
  };
});
