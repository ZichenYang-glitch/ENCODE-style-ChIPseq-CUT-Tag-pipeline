import { test as setup } from '@playwright/test';

import { AUTH_STORAGE_PATH, login } from './auth';

// Login is rate-limited per subject, so authenticate once per invocation and
// share the storage state instead of signing in before every test.
setup('authenticate as the seeded e2e administrator', async ({ page }) => {
  await login(page);
  await page.context().storageState({ path: AUTH_STORAGE_PATH });
});
