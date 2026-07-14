import { expect, test, type Page, type TestInfo } from '@playwright/test';
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

interface RuntimeManifest {
  workflowId: string;
  samplesPath: string;
  resultsConfig: Record<string, unknown>;
  cancelConfig: Record<string, unknown>;
  emptyConfig: Record<string, unknown>;
  malformedConfig: Record<string, unknown>;
  expectedQcSummary: string;
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

async function replaceCodeMirror(
  editor: ReturnType<Page['getByRole']>,
  text: string,
) {
  await editor.click();
  await editor.press('Control+A');
  await editor.press('Backspace');
  await editor.fill(text);
}

async function setWorkbenchConfig(
  page: Page,
  config: Record<string, unknown>,
) {
  await page.getByRole('tab', { name: 'Config' }).click();
  const yamlEditor = page.getByRole('textbox', { name: 'Advanced config YAML' });
  if (!(await yamlEditor.isVisible())) {
    await page.getByRole('button', { name: 'YAML mode' }).click();
  }
  await expect(yamlEditor).toBeVisible();
  await replaceCodeMirror(yamlEditor, JSON.stringify(config));
}

async function createPlannedRun(
  page: Page,
  config: Record<string, unknown>,
  options: { proveBackendRejection?: boolean } = {},
): Promise<string> {
  const runtime = manifest();
  await page.goto(`/workflows/${runtime.workflowId}/new-run`);
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();
  await setWorkbenchConfig(
    page,
    options.proveBackendRejection
      ? { ...config, stage4b: false, stage5: true }
      : config,
  );
  await page.getByRole('tab', { name: 'Samples' }).click();
  await page.getByLabel('Import samples TSV').setInputFiles(runtime.samplesPath);
  await expect(
    page.locator('[data-testid="sample-desktop-table"], [data-testid="sample-mobile-list"]')
      .filter({ visible: true }),
  ).toBeVisible();
  await page.getByRole('tab', { name: 'Review' }).click();
  await expect(
    page.getByText(/structurally ready for backend validation/i),
  ).toBeVisible();

  if (options.proveBackendRejection) {
    const rejectedResponsePromise = page.waitForResponse(
      (response) =>
        response.request().method() === 'POST' &&
        response.url().endsWith(
          `/api/v1/workflows/${runtime.workflowId}/validate`,
        ),
    );
    await page.getByRole('button', { name: 'Validate current inputs' }).click();
    const rejectedResponse = await rejectedResponsePromise;
    expect(rejectedResponse.status()).toBe(200);
    const rejected = (await rejectedResponse.json()) as {
      ok: boolean;
      snapshot: unknown;
      issues: Array<{ code: string }>;
    };
    expect(rejected.ok).toBe(false);
    expect(rejected.snapshot).toBeNull();
    expect(rejected.issues.map((issue) => issue.code)).toContain(
      'ENCODE_CONFIG_INVALID',
    );
    await expect(page.getByRole('alert')).toContainText(
      'ENCODE_CONFIG_INVALID',
    );
    await expect(
      page.getByRole('button', { name: 'Create run from validated inputs' }),
    ).toBeDisabled();
    await setWorkbenchConfig(page, config);
    await page.getByRole('tab', { name: 'Review' }).click();
  }

  const validateRequestPromise = page.waitForRequest(
    (request) =>
      request.method() === 'POST' &&
      request.url().endsWith(
        `/api/v1/workflows/${runtime.workflowId}/validate`,
      ),
  );
  const validateResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response.url().endsWith(
        `/api/v1/workflows/${runtime.workflowId}/validate`,
      ),
  );
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  const validateRequest = await validateRequestPromise;
  const validatePayload = validateRequest.postDataJSON() as {
    config: Record<string, unknown>;
    samples: Array<Record<string, string>>;
    options: Record<string, unknown>;
  };
  expect(validatePayload.config).not.toHaveProperty('samples');
  expect(JSON.stringify(validatePayload.config)).not.toContain(
    runtime.samplesPath,
  );
  expect(validatePayload.samples).toHaveLength(1);
  expect(validatePayload.samples[0]?.sample).toBe('C1');
  expect(validatePayload.options).toEqual({ strict_inputs: false });
  const validateResponse = await validateResponsePromise;
  expect(validateResponse.status()).toBe(200);
  const validation = (await validateResponse.json()) as {
    ok: boolean;
    snapshot: null | { snapshot_id: string; payload_digest: string };
  };
  expect(validation.ok).toBe(true);
  expect(validation.snapshot?.snapshot_id).toMatch(/^vsnap_[0-9a-f]{32}$/);
  expect(validation.snapshot?.payload_digest).toMatch(/^[0-9a-f]{64}$/);
  await expect(
    page.getByText(/This exact draft can create one run/i),
  ).toBeVisible();

  const preflightRequests: string[] = [];
  const recordPreflight = (request: { method(): string; url(): string }) => {
    if (request.method() === 'POST' && request.url().endsWith('/preflight')) {
      preflightRequests.push(request.url());
    }
  };
  page.on('request', recordPreflight);
  const createRequestPromise = page.waitForRequest(
    (request) =>
      request.method() === 'POST' &&
      request.url().endsWith(`/api/v1/workflows/${runtime.workflowId}/runs`),
  );
  await page
    .getByRole('button', { name: 'Create run from validated inputs' })
    .click();
  const createRequest = await createRequestPromise;
  expect(createRequest.postDataJSON()).toEqual({
    snapshot_id: validation.snapshot?.snapshot_id,
  });
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

async function executeToSucceeded(
  page: Page,
  config: Record<string, unknown>,
  logMarker: string,
  options: { proveBackendRejection?: boolean; requireRunning?: boolean } = {},
): Promise<string> {
  const runId = await createPlannedRun(page, config, options);
  const startButton = page.getByRole('button', { name: 'Start run' });
  await expect(startButton).toBeEnabled();
  await startButton.click();
  if (options.requireRunning) {
    await expect(page.getByTestId('run-status-badge')).toHaveText('running', {
      timeout: 60_000,
    });
  } else {
    await expect(page.getByTestId('run-status-badge')).toHaveText(
      /queued|running|succeeded/,
    );
  }
  await expect(page.getByTestId('run-status-badge')).toHaveText('succeeded', {
    timeout: 60_000,
  });
  await expect(page.getByText(logMarker)).toBeVisible({ timeout: 60_000 });
  await expect(
    page.getByText('Local workflow execution completed successfully.'),
  ).toBeVisible();
  return runId;
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

async function attachScreenshot(
  page: Page,
  testInfo: TestInfo,
  name: string,
) {
  const screenshot = await page.screenshot({ fullPage: true });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `${name}.png`), screenshot);
  }
  await testInfo.attach(name, {
    body: screenshot,
    contentType: 'image/png',
  });
}

async function captureArtifactViewport(
  page: Page,
  testInfo: TestInfo,
  name: string,
  width: number,
  height: number,
) {
  await page.setViewportSize({ width, height });
  if (width < 1280) await page.reload();
  await expect(page.getByRole('tab', { name: 'Artifacts' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await expect(
    page.getByRole('button', { name: 'Refresh run progress' }),
  ).toBeVisible();
  await expect(page.getByTestId('artifact-inspector')).toBeVisible();
  await expect(
    page.getByRole('button', { name: 'Download C1.qc_summary.tsv' }),
  ).toBeVisible();
  if (width < 1280) {
    await expect(page.getByTestId('artifact-inspector')).toBeInViewport();
    await expect(
      page.getByRole('button', { name: 'Back to artifact list' }),
    ).toBeVisible();
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
  await attachScreenshot(page, testInfo, `artifact-browser-${name}`);
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
  await expect(
    page.getByRole('button', { name: 'Refresh run progress' }),
  ).toBeVisible();
  await expect(page.getByTestId('qc-metric-list')).toBeVisible({
    timeout: 60_000,
  });
  await expect(
    page.getByRole('button', {
      name: 'Open source artifact for Total reads',
    }),
  ).toBeVisible();
  if (width >= 768) {
    const totalReadsRow = page
      .getByRole('table', { name: 'Indexed run QC metrics' })
      .getByRole('row')
      .filter({ hasText: 'Total reads' });
    await expect(totalReadsRow.getByText('1000', { exact: true })).toBeVisible();
  } else {
    const totalReadsItem = page
      .getByRole('list', { name: 'Indexed run QC metrics' })
      .getByRole('listitem')
      .filter({ hasText: 'Total reads' });
    await expect(totalReadsItem.getByText('1000', { exact: true })).toBeVisible();
  }
  await expectNoHorizontalOverflow(page);
  const listBox = await page.getByTestId('qc-metric-list').boundingBox();
  expect(listBox).not.toBeNull();
  expect(listBox!.x).toBeGreaterThanOrEqual(0);
  expect(listBox!.x + listBox!.width).toBeLessThanOrEqual(width + 1);
  await attachScreenshot(page, testInfo, `qc-workbench-${name}`);
}

async function captureRunHistoryViewport(
  page: Page,
  testInfo: TestInfo,
  runId: string,
  name: string,
  width: number,
  height: number,
) {
  await page.setViewportSize({ width, height });
  await page.goto('/runs');
  await expect(
    page.getByRole('navigation', { name: 'Primary' }).getByRole('link', {
      name: 'Runs',
    }),
  ).toHaveAttribute('aria-current', 'page');
  await expect(
    page.getByRole('link', { name: runId, exact: true }).filter({ visible: true }),
  ).toBeVisible({ timeout: 60_000 });
  await expect(page.getByTestId('run-history-list')).toBeVisible();
  await expectNoHorizontalOverflow(page);
  const historyBox = await page.getByTestId('run-history-list').boundingBox();
  expect(historyBox).not.toBeNull();
  expect(historyBox!.x).toBeGreaterThanOrEqual(0);
  expect(historyBox!.x + historyBox!.width).toBeLessThanOrEqual(width + 1);
  await attachScreenshot(page, testInfo, `run-history-${name}`);
}

async function verifyRunHistoryResultLinks(page: Page, runId: string) {
  const primaryNavigation = page.getByRole('navigation', { name: 'Primary' });
  await primaryNavigation.getByRole('link', { name: 'Runs' }).click();
  await expect(page).toHaveURL('/runs');
  const runLink = page
    .getByRole('link', { name: runId, exact: true })
    .filter({ visible: true });
  await expect(runLink).toBeVisible({ timeout: 60_000 });

  await page
    .getByRole('link', { name: `Open QC for run ${runId}` })
    .filter({ visible: true })
    .click();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}\\?view=qc$`));
  await expect(page.getByRole('tab', { name: 'QC' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await primaryNavigation.getByRole('link', { name: 'Runs' }).click();

  await page
    .getByRole('link', { name: `Open artifacts for run ${runId}` })
    .filter({ visible: true })
    .click();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}\\?view=artifacts$`));
  await expect(page.getByRole('tab', { name: 'Artifacts' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await primaryNavigation.getByRole('link', { name: 'Runs' }).click();

  await runLink.click();
  await expect(page).toHaveURL(`/runs/${runId}`);
  await expect(page.getByRole('tab', { name: 'Activity' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await expect(primaryNavigation.getByRole('link', { name: 'Runs' })).toHaveAttribute(
    'aria-current',
    'page',
  );
}

function processExists(pid: number): boolean {
  try {
    process.kill(pid, 0);
  } catch (error) {
    return (error as NodeJS.ErrnoException).code !== 'ESRCH';
  }
  return true;
}

test('real run exposes QC source artifact and exact download @desktop', async ({
  page,
}, testInfo) => {
  test.setTimeout(180_000);
  const runtime = manifest();
  const runId = await executeToSucceeded(
    page,
    runtime.resultsConfig,
    'browser-e2e-results',
    { proveBackendRejection: true, requireRunning: true },
  );
  await expectNoHorizontalOverflow(page);

  await verifyRunHistoryResultLinks(page, runId);
  await captureRunHistoryViewport(
    page,
    testInfo,
    runId,
    'desktop-1440x900',
    1440,
    900,
  );
  await captureRunHistoryViewport(
    page,
    testInfo,
    runId,
    'tablet-1024x768',
    1024,
    768,
  );
  await page.goto(`/runs/${runId}`);

  await page.getByRole('tab', { name: 'QC' }).click();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}\\?view=qc$`));
  await expect(
    page.getByRole('heading', { name: 'Indexed QC metrics' }),
  ).toBeVisible({ timeout: 60_000 });
  await expect(page.getByText('8 loaded of 8')).toBeVisible();

  const qcResponse = await page.request.get(
    `/api/v1/runs/${runId}/qc-metrics?limit=50`,
  );
  expect(qcResponse.status()).toBe(200);
  const qcBody = (await qcResponse.json()) as {
    qc_metrics: Array<{ metric_key: string; value: string }>;
  };
  expect(qcBody.qc_metrics).toHaveLength(8);
  expect(
    qcBody.qc_metrics.find(
      (metric) => metric.metric_key === 'sequencing.total_reads',
    )?.value,
  ).toBe('1000');

  await captureQcViewport(page, testInfo, 'desktop-1440x900', 1440, 900);
  await captureQcViewport(page, testInfo, 'tablet-1024x768', 1024, 768);
  await captureQcViewport(page, testInfo, 'mobile-390x844', 390, 844);
  await captureQcViewport(page, testInfo, 'mobile-360x800', 360, 800);

  await page
    .getByRole('button', { name: 'Open source artifact for Total reads' })
    .click();
  await expect(page).toHaveURL(
    new RegExp(`/runs/${runId}\\?view=artifacts&artifact=[A-Za-z0-9_.-]+$`),
  );
  await expect(page.getByRole('heading', { name: 'Artifact details' })).toBeVisible();
  await expect(
    page
      .getByTestId('artifact-inspector')
      .getByText('results/C1/01_qc/C1.qc_summary.tsv', { exact: true }),
  ).toBeVisible();
  await expect(
    page.getByTestId('artifact-inspector').getByText(/^run:\/\/runs\//),
  ).toBeVisible();
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
  const workspaceArtifact = readFileSync(
    join(
      runtime.workspaceRoot,
      runId,
      'results',
      'C1',
      '01_qc',
      'C1.qc_summary.tsv',
    ),
  );
  expect(workspaceArtifact.toString('utf8')).toBe(runtime.expectedQcSummary);
  expect(await rangeResponse.body()).toEqual(workspaceArtifact);

  const downloadPromise = page.waitForEvent('download');
  await page
    .getByRole('button', { name: 'Download C1.qc_summary.tsv' })
    .click();
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toBe('C1.qc_summary.tsv');
  const downloadedPath = await download.path();
  expect(downloadedPath).not.toBeNull();
  expect(readFileSync(downloadedPath!)).toEqual(workspaceArtifact);
  await expect(page.getByText('Download prepared successfully.')).toBeVisible();

  await page.reload();
  await expect(page.getByRole('tab', { name: 'Artifacts' })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await expect(
    page
      .getByTestId('artifact-inspector')
      .getByText('results/C1/01_qc/C1.qc_summary.tsv', { exact: true }),
  ).toBeVisible();

  await captureArtifactViewport(page, testInfo, 'desktop-1440x900', 1440, 900);
  await captureArtifactViewport(page, testInfo, 'tablet-1024x768', 1024, 768);
  await captureArtifactViewport(page, testInfo, 'mobile-390x844', 390, 844);
  await captureArtifactViewport(page, testInfo, 'mobile-360x800', 360, 800);

  await page.getByRole('tab', { name: 'QC' }).click();
  await page.reload();
  await expect(page).toHaveURL(new RegExp(`/runs/${runId}\\?view=qc$`));
  await expect(page.getByTestId('qc-metric-list')).toBeVisible();

  const emptyRunId = await executeToSucceeded(
    page,
    runtime.emptyConfig,
    'browser-e2e-empty',
  );
  await page.getByRole('tab', { name: 'QC' }).click();
  await expect(page).toHaveURL(new RegExp(`/runs/${emptyRunId}\\?view=qc$`));
  await expect(
    page.getByRole('heading', { name: 'No indexed QC metrics' }),
  ).toBeVisible({ timeout: 60_000 });
  await expect(
    page.getByText(
      'QC metric indexing completed successfully and produced an empty result set.',
    ),
  ).toBeVisible();
  await page.setViewportSize({ width: 390, height: 844 });
  await expectNoHorizontalOverflow(page);
  await attachScreenshot(page, testInfo, 'qc-confirmed-empty-mobile-390x844');

  await page.setViewportSize({ width: 1440, height: 900 });
  const malformedRunId = await executeToSucceeded(
    page,
    runtime.malformedConfig,
    'browser-e2e-malformed',
  );
  await page.getByRole('tab', { name: 'QC' }).click();
  await expect(page).toHaveURL(
    new RegExp(`/runs/${malformedRunId}\\?view=qc$`),
  );
  await expect(
    page.getByRole('heading', { name: 'QC metric indexing failed' }),
  ).toBeVisible({ timeout: 60_000 });
  await expect(page.getByTestId('run-status-badge')).toHaveText('succeeded');
  await expect(page.getByText('sample\tassay')).not.toBeVisible();
  await expectNoHorizontalOverflow(page);
  await attachScreenshot(page, testInfo, 'qc-indexing-failed-desktop-1440x900');
});

test('real successful run is discoverable from mobile run history @mobile', async ({
  page,
}, testInfo) => {
  test.setTimeout(180_000);
  const runtime = manifest();
  const runId = await executeToSucceeded(
    page,
    runtime.resultsConfig,
    'browser-e2e-results',
  );

  await verifyRunHistoryResultLinks(page, runId);
  await captureRunHistoryViewport(
    page,
    testInfo,
    runId,
    'mobile-390x844',
    390,
    844,
  );
  await captureRunHistoryViewport(
    page,
    testInfo,
    runId,
    'mobile-360x800',
    360,
    800,
  );
});

test('running cancellation stays requested until worker acknowledgement @mobile', async ({
  page,
}) => {
  test.setTimeout(180_000);
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
  expect(
    existsSync(
      join(
        runtime.workspaceRoot,
        runId,
        'results',
        'C1',
        '01_qc',
        'C1.qc_summary.tsv',
      ),
    ),
  ).toBe(false);

  await page.reload();
  await expect(page.getByTestId('run-status-badge')).toHaveText('cancelled');
  await expect(page.getByText('browser-e2e-cancel-entered')).toBeVisible();
});
