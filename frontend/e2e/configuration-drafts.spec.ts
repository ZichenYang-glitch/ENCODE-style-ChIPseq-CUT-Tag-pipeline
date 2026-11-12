import { execFileSync } from 'node:child_process';
import { expect, test, type Page, type Response } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';

const runtimeRoot = process.env.ENCODE_PIPELINE_E2E_ROOT!;

function workflowId(): string {
  const manifest = JSON.parse(
    readFileSync(join(runtimeRoot, 'runtime.json'), 'utf8'),
  ) as { workflowId: string };
  return manifest.workflowId;
}

function seedProject(): void {
  execFileSync('python3', [join(process.cwd(), 'e2e/seed_authoring_fixture.py')], {
    env: process.env,
    stdio: 'pipe',
  });
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
    page.getByRole('heading', { level: 1, name: 'Configuration drafts' }),
  ).toBeVisible();
}

test('an account without a project cannot create a draft @desktop', async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await openDrafts(page);
  await expect(page.getByText(/no projects are available/i)).toBeVisible();
  await expect(page.getByRole('button', { name: 'New draft' })).toBeDisabled();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
});

test('a controlled draft can be saved, reread, and validated on desktop @desktop', async ({
  page,
}) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await exerciseControlledDraft(page);
  await expectNoHorizontalOverflow(page);
});

test('a controlled draft stays operable on a phone @mobile', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await exerciseControlledDraft(page);
  await expectNoHorizontalOverflow(page);
  const project = page.getByRole('combobox', { name: 'Project' });
  const box = await project.boundingBox();
  expect(box).not.toBeNull();
  expect(box!.x).toBeGreaterThanOrEqual(0);
  expect(box!.x + box!.width).toBeLessThanOrEqual(390);
});

async function exerciseControlledDraft(page: Page) {
  seedProject();
  await openDrafts(page);
  await expect(page.getByRole('combobox', { name: 'Project' })).toContainText(
    'Browser authoring',
  );
  const created = page.waitForResponse(
    (response) =>
      response.url().includes('/api/v1/agent/authoring-sessions') &&
      response.request().method() === 'POST' &&
      !response.url().includes('/draft') &&
      !response.url().includes('/validate') &&
      !response.url().includes('/connection'),
  );
  await page.getByRole('button', { name: 'New draft' }).click();
  const createResponse = await created;
  expect(createResponse.ok()).toBeTruthy();
  const createdBody = (await createResponse.json()) as { draft_revision: number };
  expect(createdBody.draft_revision).toBe(0);
  await expect(page.getByText('Assistant available')).toHaveCount(0);

  const threads = page.getByRole('spinbutton', { name: /thread/i });
  await expect(threads).toBeVisible();
  await threads.click();
  await threads.fill('');
  await threads.pressSequentially('4');
  await expect(page.getByRole('button', { name: 'Save edits' })).toBeEnabled();
  const saved = page.waitForResponse(
    (response) =>
      response.url().includes('/draft') && response.request().method() === 'PUT',
  );
  await page.getByRole('button', { name: 'Save edits' }).click();
  const saveResponse = await saved;
  expect(saveResponse.ok()).toBeTruthy();
  const saveBody = (await saveResponse.json()) as {
    draft_revision: number;
    draft: { config: { threads?: number } };
  };
  expect(saveBody.draft_revision).toBeGreaterThan(0);
  expect(saveBody.draft.config.threads).toBe(4);

  await page.reload();
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  await expect(page.getByRole('spinbutton', { name: /thread/i })).toHaveValue('4');

  const sample = page.getByRole('checkbox', { name: 'Browser sample' });
  if (await sample.isVisible()) await sample.check();
  const reference = page.getByRole('combobox', { name: 'Reference profile' });
  if (await reference.isEnabled()) {
    const option = reference.locator('option').nth(1);
    const value = await option.getAttribute('value');
    if (value) await reference.selectOption(value);
  }
  if (await page.getByRole('button', { name: 'Save edits' }).isEnabled()) {
    const selectionSaved = page.waitForResponse(
      (response) =>
        response.url().includes('/draft') && response.request().method() === 'PUT',
    );
    await page.getByRole('button', { name: 'Save edits' }).click();
    expect((await selectionSaved).ok()).toBeTruthy();
  }

  const validated = page.waitForResponse(
    (response) =>
      response.url().includes('/validate') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Validate draft' }).click();
  await assertRealValidation(page, await validated);

  const connection = page.waitForResponse(
    (response) =>
      response.url().includes('/connection') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Connect assistant' }).click();
  expect((await connection).status()).toBe(503);
  await expect(page.getByText('Assistant unavailable')).toBeVisible();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  await expect(page.getByRole('button', { name: 'Save edits' })).toBeVisible();
  await expect(page.getByRole('button', { name: 'Validate draft' })).toBeVisible();
}

async function assertRealValidation(page: Page, response: Response) {
  const body = (await response.json()) as {
    ok: boolean;
    draft_revision: number;
    snapshot: { snapshot_id?: string; expires_at?: string } | null;
    issues?: { code: string; message: string }[];
  };
  expect(body.draft_revision).toBeGreaterThan(0);
  if (body.ok && body.snapshot?.snapshot_id) {
    await expect(page.getByText(body.snapshot.snapshot_id)).toBeVisible();
    await expect(page.getByText(`Revision ${body.draft_revision}`)).toBeVisible();
    return;
  }
  expect(body.snapshot).toBeNull();
  const issue = body.issues?.[0];
  expect(issue?.code).toBeTruthy();
  const message = issue?.message ?? '';
  const safe = message.trim().length > 0 && !/[\\/\n]|traceback/i.test(message);
  if (safe) {
    await expect(page.getByText(message, { exact: false })).toBeVisible();
  } else {
    await expect(page.getByText(/validation reported problems|not accepted/i)).toBeVisible();
  }
  await expect(page.getByText(/validation finished/i)).toHaveCount(0);
}
