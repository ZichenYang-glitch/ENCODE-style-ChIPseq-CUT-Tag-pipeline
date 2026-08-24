import { expect, type Page } from '@playwright/test';

export const AUTH_STORAGE_PATH = 'e2e/.auth/e2e-admin.json';

export const e2eAdmin = {
  username: 'e2e-admin',
  password: 'e2e playwright admin password',
};

export async function login(page: Page): Promise<void> {
  await page.goto('/login');
  await page.getByLabel('Username').fill(e2eAdmin.username);
  await page.getByLabel('Password').fill(e2eAdmin.password);
  await page.getByRole('button', { name: 'Sign in' }).click();
  await expect(page).toHaveURL(/\/workflows/);
}
