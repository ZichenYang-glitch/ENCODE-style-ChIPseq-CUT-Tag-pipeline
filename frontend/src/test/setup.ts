import '@testing-library/jest-dom/vitest';
import { vi } from 'vitest';

// Node's Request rejects a DOM AbortSignal. React Router builds one for every
// client navigation, so the location never changes under jsdom on this runtime.
const OriginalRequest = globalThis.Request;
globalThis.Request = class Request extends OriginalRequest {
  constructor(input: RequestInfo | URL, init?: RequestInit) {
    if (init && 'signal' in init) {
      const { signal: _signal, ...rest } = init;
      super(input, rest);
      return;
    }
    super(input, init);
  }
};

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
