import { expect, test, type Page } from '@playwright/test';
import { existsSync, readFileSync } from 'node:fs';
import { join } from 'node:path';

interface RuntimeManifest {
  workflowId: string;
  samplesPath: string;
  successConfig: Record<string, unknown>;
  cancelConfig: Record<string, unknown>;
  runtimeRoot: string;
  workspaceRoot: string;
  markerRoot: string;
  queueName: string;
}

const runtimeRoot = process.env.ENCODE_PIPELINE_E2E_ROOT!;

function manifest(): RuntimeManifest {
  return JSON.parse(
    readFileSync(join(runtimeRoot, 'runtime.json'), 'utf8'),
  ) as RuntimeManifest;
}

async function createPlannedRun(
  page: Page,
  config: Record<string, unknown>,
): Promise<string> {
  const runtime = manifest();
  await page.goto(`/workflows/${runtime.workflowId}`);
  await page.getByLabel('Config (JSON)').fill(JSON.stringify(config));
  await page.getByLabel('Samples (path string)').fill(runtime.samplesPath);
  await page.getByRole('button', { name: 'Validate' }).click();
  await expect(page.getByTestId('create-run-button')).toBeEnabled();

  const preflightRequests: string[] = [];
  const recordPreflight = (request: { method(): string; url(): string }) => {
    if (request.method() === 'POST' && request.url().endsWith('/preflight')) {
      preflightRequests.push(request.url());
    }
  };
  page.on('request', recordPreflight);
  await page.getByRole('button', { name: 'Create run' }).click();
  await expect(page).toHaveURL(/\/runs\/[^/]+$/);
  const runId = page.url().split('/').at(-1);
  expect(runId).toBeTruthy();
  await expect(page.getByTestId('run-status-badge')).toHaveText('planned', {
    timeout: 60_000,
  });
  page.off('request', recordPreflight);
  expect(preflightRequests).toHaveLength(1);
  return runId!;
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

function processExists(pid: number): boolean {
  try {
    process.kill(pid, 0);
  } catch (error) {
    return (error as NodeJS.ErrnoException).code !== 'ESRCH';
  }
  return true;
}

test('real durable execution succeeds, persists logs, and survives reload @desktop', async ({
  page,
}) => {
  const runtime = manifest();
  const runId = await createPlannedRun(page, runtime.successConfig);

  const startButton = page.getByRole('button', { name: 'Start run' });
  await expect(startButton).toBeEnabled();
  await startButton.click();
  await expect(page.getByTestId('run-status-badge')).toHaveText(
    /queued|running|succeeded/,
  );
  await expect(page.getByTestId('run-status-badge')).toHaveText('succeeded', {
    timeout: 60_000,
  });
  await expect(page.getByText('browser-e2e-success')).toBeVisible();
  await expect(
    page.getByText('Local workflow execution completed successfully.'),
  ).toBeVisible();
  await expectNoHorizontalOverflow(page);

  await page.reload();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}$`));
  await expect(page.getByTestId('run-status-badge')).toHaveText('succeeded');
  await expect(page.getByText('browser-e2e-success')).toBeVisible();
});

test('running cancellation stays requested until worker acknowledgement @mobile', async ({
  page,
}) => {
  const runtime = manifest();
  const runId = await createPlannedRun(page, runtime.cancelConfig);

  const startButton = page.getByRole('button', { name: 'Start run' });
  await expect(startButton).toBeEnabled();
  await startButton.click();
  await expect(page.getByTestId('run-status-badge')).toHaveText('running', {
    timeout: 60_000,
  });
  await expect(page.getByText('browser-e2e-cancel-entered')).toBeVisible({
    timeout: 60_000,
  });

  const cancelResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response.url().endsWith(`/api/v1/runs/${runId}/cancel`),
  );
  await page.getByRole('button', { name: 'Cancel run' }).click();
  const cancelResponse = await cancelResponsePromise;
  expect(cancelResponse.status()).toBe(202);
  const cancelBody = (await cancelResponse.json()) as { run: { status: string } };
  expect(cancelBody.run.status).toBe('running');

  await expect(page.getByTestId('cancellation-requested')).toBeVisible();
  await expect(page.getByTestId('run-status-badge')).toHaveText('running');
  await expect(page.getByTestId('run-status-badge')).not.toHaveText('cancelling');
  await expectNoHorizontalOverflow(page);

  await expect(page.getByTestId('run-status-badge')).toHaveText('cancelled', {
    timeout: 60_000,
  });
  await expect(page.getByTestId('cancellation-requested')).not.toBeVisible();
  await expect(page.getByText('cancellation_acknowledged')).toBeVisible();

  const processEvidencePath = join(runtime.markerRoot, 'browser-processes.json');
  await expect.poll(() => existsSync(processEvidencePath)).toBe(true);
  const evidence = JSON.parse(readFileSync(processEvidencePath, 'utf8')) as {
    helper_pid: number;
    child_pid: number;
  };
  await expect
    .poll(() => [evidence.helper_pid, evidence.child_pid].some(processExists))
    .toBe(false);
  expect(existsSync(join(runtime.markerRoot, 'browser-helper-completed'))).toBe(false);
  expect(
    existsSync(join(runtime.workspaceRoot, runId, 'result', 'complete.txt')),
  ).toBe(false);

  await page.reload();
  await expect(page.getByTestId('run-status-badge')).toHaveText('cancelled');
  await expect(page.getByText('browser-e2e-cancel-entered')).toBeVisible();
});
