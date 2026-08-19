import '@testing-library/jest-dom/vitest';
import { vi } from 'vitest';

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
