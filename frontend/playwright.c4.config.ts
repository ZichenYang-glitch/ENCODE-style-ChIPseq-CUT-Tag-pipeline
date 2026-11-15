import { defineConfig, devices } from '@playwright/test';

import { AUTH_STORAGE_PATH } from './e2e/auth';

const noProxy = [process.env.NO_PROXY, '127.0.0.1', 'localhost', '::1']
  .filter(Boolean)
  .join(',');
process.env.NO_PROXY = noProxy;
process.env.no_proxy = noProxy;

if (
  !process.env.ENCODE_PIPELINE_E2E_ROOT ||
  !process.env.ENCODE_PIPELINE_E2E_OWNER ||
  !process.env.HELIXWEAVE_C4_JUNIT
) {
  throw new Error('Run the C4 browser check through scripts/run-c4-playwright.mjs.');
}

export default defineConfig({
  testDir: './e2e',
  testMatch: ['**/auth.setup.ts', '**/*c4-native.spec.ts'],
  fullyParallel: false,
  workers: 1,
  retries: 0,
  timeout: 180_000,
  expect: { timeout: 30_000 },
  reporter: [
    ['list'],
    ['junit', { outputFile: process.env.HELIXWEAVE_C4_JUNIT }],
  ],
  use: {
    baseURL: 'http://127.0.0.1:4173',
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
    video: 'off',
  },
  webServer: {
    command: 'exec python3 ../test/browser/platform_runtime.py',
    url: 'http://127.0.0.1:4173',
    timeout: 180_000,
    reuseExistingServer: false,
    stdout: 'pipe',
    stderr: 'pipe',
  },
  projects: [
    {
      name: 'setup',
      testMatch: /auth\.setup\.ts/,
    },
    {
      name: 'c4-desktop',
      testMatch: /c4-native\.spec\.ts/,
      dependencies: ['setup'],
      use: {
        ...devices['Desktop Chrome'],
        viewport: { width: 1440, height: 900 },
        storageState: AUTH_STORAGE_PATH,
      },
    },
  ],
});
