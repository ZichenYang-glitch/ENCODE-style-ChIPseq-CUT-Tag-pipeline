// Invoked by the isolated real-HTTP fixture in test/api/peri_run_actions_browser.py.
import { readFileSync } from 'node:fs';
import { chromium, expect } from '@playwright/test';

const fixture = JSON.parse(readFileSync(0, 'utf8'));
const browser = await chromium.launch();

try {
  const context = await browser.newContext({ baseURL: fixture.url });
  await context.addCookies([fixture.session, fixture.csrf].map(([name, value]) => ({
    name, value, url: fixture.url,
  })));
  const page = await context.newPage();
  page.setDefaultTimeout(8000);
  const actionRoot = `/api/v1/agent/authoring-sessions/${fixture.sid}/run-actions`;
  const runPath = (id) => `/api/v1/runs/${id}`;
  const apiResponse = (path, method) => page.waitForResponse(
    (response) => new URL(response.url()).pathname === path &&
      response.request().method() === method,
  );

  await page.goto(`/workflows/${fixture.workflow}/configuration-drafts`);
  await page.getByRole('combobox', { name: 'Project' }).selectOption(fixture.projectId);
  await page.getByRole('button', { name: /Peri run action fixture/ }).click();
  const connect = page.getByRole('button', { name: 'Connect assistant' });
  await connect.click();
  await expect(page.getByText('Assistant available')).toBeVisible();
  const validatePath = actionRoot.replace('/run-actions', '/validate');
  const validation = apiResponse(validatePath, 'POST');
  await page.getByRole('button', { name: 'Validate draft' }).click();
  const validated = (await (await validation).json()).snapshot.snapshot_id;
  await page.getByLabel('Tags (optional, one key=value per line)').fill('purpose=browser');
  await page.getByRole('button', { name: 'Review run action' }).click();

  const review = page.getByRole('group', { name: 'Review Peri run action' });
  await expect(review.getByText(validated)).toBeVisible();
  await expect(review.getByText('{"purpose":"browser"}')).toBeVisible();
  await expect(page.getByRole('link', { name: /^Open run / })).toHaveCount(0);

  // Allow the actual write and lose only its response at the browser boundary.
  let first;
  await page.route(`**${actionRoot}/confirm`, async (route) => {
    const response = await route.fetch();
    first = { status: response.status(), body: await response.json() };
    await route.abort('failed');
  }, { times: 1 });
  await page.getByRole('button', { name: 'Confirm run creation' }).click();
  const retry = page.getByRole('button', { name: 'Retry confirmation' });
  await expect(retry).toBeVisible();
  expect(first.status).toBe(201);

  const readCreated = apiResponse(runPath(first.body.run_id), 'GET');
  const replay = apiResponse(`${actionRoot}/confirm`, 'POST');
  await retry.click();
  const replayed = await replay;
  expect(replayed.status()).toBe(200);
  expect(await replayed.json()).toEqual(first.body);
  expect((await (await readCreated).json()).run.status).toBe('created');
  await expect(page.getByRole('link', { name: `Open run ${first.body.run_id}` }))
    .toHaveAttribute('href', `/runs/${first.body.run_id}`);

  await page.getByLabel('Cancel a specific run').check();
  await page.getByLabel('Run ID', { exact: true }).fill(fixture.runB);
  await page.getByRole('button', { name: 'Review run action' }).click();
  await expect(review.getByText(fixture.runB)).toBeVisible();
  await page.getByRole('button', { name: 'Change request' }).click();
  await page.getByLabel('Run ID', { exact: true }).fill(fixture.runA);
  await page.getByRole('button', { name: 'Review run action' }).click();
  await expect(review.getByText(fixture.runA)).toBeVisible();
  await expect(review.getByText(fixture.runB)).toHaveCount(0);

  const stop = apiResponse(`${actionRoot}/confirm`, 'POST');
  const readRunning = apiResponse(runPath(fixture.runA), 'GET');
  await page.getByRole('button', { name: 'Confirm cancellation' }).click();
  expect((await stop).status()).toBe(202);
  expect((await (await readRunning).json()).run.status).toBe('running');
  await expect(page.getByText('Stop requested.')).toBeVisible();
  await expect(page.getByText('Run status: running.')).toBeVisible();
  const untouched = await context.request.get(runPath(fixture.runB));
  expect((await untouched.json()).run.status).toBe('created');

  expect((await context.request.post('/fixture/acknowledge-stop')).ok()).toBeTruthy();
  await expect(page.getByText('Run status: cancelled.')).toBeVisible();
  await expect(page.getByText('Stop requested.')).toHaveCount(0);
} finally {
  await browser.close();
}
