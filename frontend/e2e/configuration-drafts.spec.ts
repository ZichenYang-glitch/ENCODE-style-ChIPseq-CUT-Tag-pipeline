import { expect, test, type Page } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';

const runtimeRoot = process.env.ENCODE_PIPELINE_E2E_ROOT!;

function workflowId(): string {
  const manifest = JSON.parse(
    readFileSync(join(runtimeRoot, 'runtime.json'), 'utf8'),
  ) as { workflowId: string };
  return manifest.workflowId;
}

async function expectNoHorizontalOverflow(page: Page) {
  await expect
    .poll(() =>
      page.evaluate(
        () =>
          document.documentElement.scrollWidth <=
          document.documentElement.clientWidth,
      ),
    )
    .toBe(true);
}

async function openDrafts(page: Page) {
  await page.goto(`/workflows/${workflowId()}/configuration-drafts`);
  await expect(
    page.getByRole('heading', { name: 'Configuration drafts', exact: true }),
  ).toBeVisible();
  await expect(page.getByRole('combobox', { name: 'Project' })).toBeVisible();
  await expect(page.getByRole('button', { name: 'New draft' })).toBeVisible();
  await expect(page.getByText(/reviewed configuration-draft interface/i)).toHaveCount(0);
  await expect(page.locator('input[type="file"]')).toHaveCount(0);
  await expect(page.getByText('Assistant available')).toHaveCount(0);
}

test('configuration drafts load from the server on desktop @desktop', async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await openDrafts(page);
  await expectNoHorizontalOverflow(page);
  await exerciseDraft(page);
});

test('configuration drafts stay usable on a phone @mobile', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await openDrafts(page);
  await expectNoHorizontalOverflow(page);
  const project = page.getByRole('combobox', { name: 'Project' });
  await expect(project).toBeVisible();
  const box = await project.boundingBox();
  expect(box).not.toBeNull();
  expect(box!.x).toBeGreaterThanOrEqual(0);
  expect(box!.x + box!.width).toBeLessThanOrEqual(390);
  await exerciseDraft(page);
});

async function exerciseDraft(page: Page) {
  const create = page.getByRole('button', { name: 'New draft' });
  if (await create.isDisabled()) {
    await expect(page.getByText(/no projects are available/i)).toBeVisible();
    return;
  }
  const responsePromise = page.waitForResponse(
    (response) =>
      response.url().includes('/api/v1/agent/authoring-sessions') &&
      response.request().method() === 'POST',
  );
  await create.click();
  const response = await responsePromise;
  if (!response.ok()) {
    const failure = page.getByRole('alert');
    await expect(failure).toBeVisible();
    await expect(failure).not.toContainText(/traceback|\/home\/|exception/i);
    return;
  }
  await expect(page.getByRole('button', { name: /revision /i }).first()).toBeVisible();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  const connect = page.getByRole('button', { name: 'Connect assistant' });
  await expect(connect).toBeEnabled();
  const connectionPromise = page.waitForResponse(
    (next) =>
      next.url().includes('/connection') && next.request().method() === 'POST',
  );
  await connect.click();
  const connection = await connectionPromise;
  expect(connection.status()).toBe(503);
  await expect(page.getByText('Assistant unavailable')).toBeVisible();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  await expect(page.getByRole('button', { name: 'Save edits' })).toBeVisible();
  await expect(page.getByRole('button', { name: 'Validate draft' })).toBeVisible();
}
