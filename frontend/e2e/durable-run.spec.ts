import { expect, test, type Page, type TestInfo } from '@playwright/test';
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
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

async function captureArtifactViewport(
  page: Page,
  testInfo: TestInfo,
  name: string,
  width: number,
  height: number,
) {
  await page.setViewportSize({ width, height });
  if (width < 1280) {
    await page.reload();
  }
  await expect(page.getByRole('tab', { name: 'Artifacts' })).toBeVisible();
  await expect(page.getByRole('button', { name: 'Refresh run progress' })).toBeVisible();
  await expect(page.getByTestId('artifact-inspector')).toBeVisible();
  await expect(
    page.getByRole('button', { name: 'Download result_manifest.tsv' }),
  ).toBeVisible();
  if (width < 1280) {
    await expect(page.getByTestId('artifact-inspector')).toBeInViewport();
    await expect(page.getByRole('button', { name: 'Back to artifact list' })).toBeVisible();
  }
  await expectNoHorizontalOverflow(page);

  const list = await page.getByTestId('artifact-list').boundingBox();
  const inspector = await page.getByTestId('artifact-inspector').boundingBox();
  expect(list).not.toBeNull();
  expect(inspector).not.toBeNull();
  if (width >= 1280) {
    expect(inspector!.x).toBeGreaterThan(list!.x + list!.width - 2);
  } else {
    expect(inspector!.y).toBeGreaterThanOrEqual(list!.y + list!.height - 2);
  }

  const screenshot = await page.screenshot({ fullPage: true });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `${name}.png`), screenshot);
  }
  await testInfo.attach(`artifact-browser-${name}`, {
    body: screenshot,
    contentType: 'image/png',
  });
}

async function captureQcViewport(
  page: Page,
  testInfo: TestInfo,
  name: string,
  width: number,
  height: number,
) {
  await page.setViewportSize({ width, height });
  await page.reload();
  await expect(page.getByRole('tab', { name: 'QC' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await expect(page.getByRole('button', { name: 'Refresh run progress' })).toBeVisible();
  const emptyState = page.getByRole('heading', { name: 'No indexed QC metrics' });
  await expect(emptyState).toBeVisible({ timeout: 60_000 });
  await expect(
    page.getByText(
      'QC metric indexing completed successfully and produced an empty result set.',
    ),
  ).toBeVisible();
  await expectNoHorizontalOverflow(page);

  const statusBox = await emptyState.boundingBox();
  expect(statusBox).not.toBeNull();
  expect(statusBox!.x).toBeGreaterThanOrEqual(0);
  expect(statusBox!.x + statusBox!.width).toBeLessThanOrEqual(width + 1);

  const screenshot = await page.screenshot({ fullPage: true });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `qc-${name}.png`), screenshot);
  }
  await testInfo.attach(`qc-workbench-${name}`, {
    body: screenshot,
    contentType: 'image/png',
  });
}

function processExists(pid: number): boolean {
  try {
    process.kill(pid, 0);
  } catch (error) {
    return (error as NodeJS.ErrnoException).code !== 'ESRCH';
  }
  return true;
}

test('real durable execution succeeds, exposes artifacts, and survives deep-link reload @desktop', async ({
  page,
}, testInfo) => {
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

  await page.getByRole('tab', { name: 'Artifacts' }).click();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}\\?view=artifacts$`));
  const manifestArtifact = page.getByRole('button', {
    name: 'Open artifact result_manifest.tsv',
  }).first();
  await expect(manifestArtifact).toBeVisible({ timeout: 60_000 });
  await manifestArtifact.click();
  await expect(page).toHaveURL(
    new RegExp(`/runs/${runId}\\?view=artifacts&artifact=[A-Za-z0-9_.-]+$`),
  );
  await expect(page.getByRole('heading', { name: 'Artifact details' })).toBeVisible();
  await expect(
    page
      .getByTestId('artifact-inspector')
      .getByText('results/multiqc/result_manifest.tsv', { exact: true }),
  ).toBeVisible();
  await expect(page.getByTestId('artifact-inspector').getByText(/^run:\/\/runs\//)).toBeVisible();
  await expect(page.getByRole('button', { name: 'Copy opaque uri' })).toBeVisible();
  const artifactId = new URL(page.url()).searchParams.get('artifact');
  expect(artifactId).toBeTruthy();
  const rangeResponse = await page.request.get(
    `/api/v1/runs/${runId}/artifacts/${artifactId}/download`,
    { headers: { Range: 'bytes=0-3' } },
  );
  expect(rangeResponse.status()).toBe(200);
  expect(rangeResponse.headers()['content-range']).toBeUndefined();
  expect(rangeResponse.headers()['accept-ranges']).toBe('none');
  expect((await rangeResponse.body()).toString('utf8')).toBe(
    'output_type\tstatus\tpath\n' +
      'result_manifest\tpresent\tresults/multiqc/result_manifest.tsv\n',
  );
  const downloadPromise = page.waitForEvent('download');
  await page
    .getByRole('button', { name: 'Download result_manifest.tsv' })
    .click();
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toBe('result_manifest.tsv');
  const downloadedPath = await download.path();
  expect(downloadedPath).not.toBeNull();
  expect(readFileSync(downloadedPath!, 'utf8')).toBe(
    'output_type\tstatus\tpath\n' +
      'result_manifest\tpresent\tresults/multiqc/result_manifest.tsv\n',
  );
  await expect(page.getByText('Download prepared successfully.')).toBeVisible();

  await page.reload();
  await expect(page.getByRole('tab', { name: 'Artifacts' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await expect(
    page
      .getByTestId('artifact-inspector')
      .getByText('results/multiqc/result_manifest.tsv', { exact: true }),
  ).toBeVisible();

  await captureArtifactViewport(page, testInfo, 'desktop-1440x900', 1440, 900);
  await captureArtifactViewport(page, testInfo, 'tablet-1024x768', 1024, 768);
  await captureArtifactViewport(page, testInfo, 'mobile-390x844', 390, 844);
  await captureArtifactViewport(page, testInfo, 'mobile-360x800', 360, 800);

  await page.getByRole('tab', { name: 'QC' }).click();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}\\?view=qc$`));
  await expect(
    page.getByRole('heading', { name: 'No indexed QC metrics' }),
  ).toBeVisible({ timeout: 60_000 });
  await page.reload();
  await expect(page.getByRole('tab', { name: 'QC' })).toHaveAttribute(
    'aria-selected',
    'true',
  );

  await captureQcViewport(page, testInfo, 'desktop-1440x900', 1440, 900);
  await captureQcViewport(page, testInfo, 'tablet-1024x768', 1024, 768);
  await captureQcViewport(page, testInfo, 'mobile-390x844', 390, 844);
  await captureQcViewport(page, testInfo, 'mobile-360x800', 360, 800);

  await page.getByRole('tab', { name: 'Activity' }).click();
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
  expect(
    existsSync(
      join(
        runtime.workspaceRoot,
        runId,
        'results',
        'multiqc',
        'result_manifest.tsv',
      ),
    ),
  ).toBe(false);

  await page.reload();
  await expect(page.getByTestId('run-status-badge')).toHaveText('cancelled');
  await expect(page.getByText('browser-e2e-cancel-entered')).toBeVisible();
});
