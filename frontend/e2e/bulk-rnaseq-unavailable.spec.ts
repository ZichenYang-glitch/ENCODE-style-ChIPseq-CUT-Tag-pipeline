import {
  expect,
  test,
  type Locator,
  type Page,
  type TestInfo,
} from '@playwright/test';
import { createHash } from 'node:crypto';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

interface RuntimeManifest {
  bulkWorkflowId: string;
  bulkSamplesPath: string;
  bulkConfig: Record<string, unknown>;
  bulkOptions: Record<string, unknown>;
  bulkExpectedExecution: 'available' | 'not_configured';
  bulkRequiredArtifactOutputTypes: string[];
  bulkRequiredQcMetricKeys: string[];
}

const runtimeRoot = process.env.ENCODE_PIPELINE_E2E_ROOT!;

function manifest(): RuntimeManifest {
  return JSON.parse(
    readFileSync(join(runtimeRoot, 'runtime.json'), 'utf8'),
  ) as RuntimeManifest;
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

async function capture(page: Page, testInfo: TestInfo, name: string) {
  await expectNoHorizontalOverflow(page);
  const screenshot = await page.screenshot({ fullPage: true, scale: 'css' });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `${name}.png`), screenshot);
  }
  await testInfo.attach(name, { body: screenshot, contentType: 'image/png' });
}

async function captureElement(
  page: Page,
  locator: Locator,
  testInfo: TestInfo,
  name: string,
) {
  await expectNoHorizontalOverflow(page);
  await expect(locator).toBeVisible();
  const screenshot = await locator.screenshot({ scale: 'css' });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `${name}.png`), screenshot);
  }
  await testInfo.attach(name, { body: screenshot, contentType: 'image/png' });
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

async function openBulkAuthoring(page: Page, fixture: RuntimeManifest) {
  await page.goto('/workflows');
  await expect(page.getByRole('heading', { name: 'Workflows' })).toBeVisible();
  const bulkCard = page.getByRole('button', { name: /Bulk RNA-seq/ });
  await expect(bulkCard).toContainText('bulk-rnaseq');
  await expect(bulkCard).toContainText(
    fixture.bulkExpectedExecution === 'available' ? 'Runnable' : 'Not configured',
  );
  await bulkCard.click();
  await expect(page).toHaveURL(`/workflows/${fixture.bulkWorkflowId}`);
  await expect(
    page.getByText(/^Upstream: nf-core\/rnaseq 3\.26\.0/),
  ).toBeVisible();
  await expect(
    page.getByText(
      fixture.bulkExpectedExecution === 'available'
        ? 'WORKFLOW_EXECUTION_READY'
        : 'WORKFLOW_EXECUTION_NOT_CONFIGURED',
    ),
  ).toBeVisible();
  const authorInputs = page.getByRole('link', { name: 'Author inputs' });
  const newAnalysis = page.getByRole('link', { name: 'New analysis' });
  await expect(authorInputs).toHaveAttribute(
    'href',
    `/workflows/${fixture.bulkWorkflowId}/new-run`,
  );
  await expect(newAnalysis).toHaveAttribute(
    'href',
    `/workflows/${fixture.bulkWorkflowId}/new-run`,
  );
  return newAnalysis;
}

async function populateBulkDraft(page: Page, fixture: RuntimeManifest) {
  await page.getByRole('button', { name: 'YAML mode' }).click();
  const yamlEditor = page.getByRole('textbox', { name: 'Advanced config YAML' });
  await replaceCodeMirror(yamlEditor, JSON.stringify(fixture.bulkConfig, null, 2));

  await page.getByRole('tab', { name: 'Samples' }).click();
  await page
    .getByLabel('Import samples TSV')
    .setInputFiles(fixture.bulkSamplesPath);
  await expect(page.getByText('3 rows', { exact: true })).toBeVisible();
  await expect(page.getByTestId('sample-transport-issue')).toHaveCount(0);
}

interface ArtifactRecord {
  artifact_id: string;
  name: string;
  revision: string;
  output_type: string;
  size_bytes: number;
}

interface QcMetricRecord {
  metric_id: string;
  metric_key: string;
  value: string;
  unit: string;
}

interface PublicWorkflowIdentity {
  schemaVersion: string;
  availabilityReason: string;
  upstreamName: string;
  upstreamVersion: string;
  upstreamRevision: string;
}

async function loadPublicWorkflowIdentity(
  page: Page,
  workflowId: string,
): Promise<PublicWorkflowIdentity> {
  const response = await page.request.get(`/api/v1/workflows/${workflowId}`);
  expect(response.status()).toBe(200);
  const body = (await response.json()) as {
    ok: boolean;
    workflow_id: string;
    workflow: null | {
      schema_version: string;
      availability: { reason_code: string };
      upstream_identity: null | {
        name: string;
        version: string;
        revision: string;
      };
    };
  };
  expect(body.ok).toBe(true);
  expect(body.workflow_id).toBe(workflowId);
  expect(body.workflow).not.toBeNull();
  expect(body.workflow!.schema_version).toMatch(/^\d+\.\d+\.\d+$/);
  expect(body.workflow!.availability.reason_code).toBe('WORKFLOW_EXECUTION_READY');
  expect(body.workflow!.upstream_identity).not.toBeNull();
  expect(body.workflow!.upstream_identity).toEqual({
    name: 'nf-core/rnaseq',
    version: '3.26.0',
    revision: 'e7ca46272c8f9d5ceee3f71759f4ba551d3217a4',
  });
  return {
    schemaVersion: body.workflow!.schema_version,
    availabilityReason: body.workflow!.availability.reason_code,
    upstreamName: body.workflow!.upstream_identity!.name,
    upstreamVersion: body.workflow!.upstream_identity!.version,
    upstreamRevision: body.workflow!.upstream_identity!.revision,
  };
}

async function loadArtifacts(page: Page, runId: string) {
  const artifacts: ArtifactRecord[] = [];
  const cursors = new Set<string>();
  let generation: string | null = null;
  let after: string | null = null;
  for (let pageNumber = 0; pageNumber < 20; pageNumber += 1) {
    const parameters = new URLSearchParams({ limit: '100' });
    if (generation !== null) parameters.set('generation', generation);
    if (after !== null) parameters.set('after', after);
    const response = await page.request.get(
      `/api/v1/runs/${runId}/artifacts?${parameters.toString()}`,
    );
    expect(response.status()).toBe(200);
    const body = (await response.json()) as {
      ok: boolean;
      run_id: string;
      artifact_generation: string;
      artifacts?: ArtifactRecord[];
      next_cursor?: string | null;
    };
    expect(body.ok).toBe(true);
    expect(body.run_id).toBe(runId);
    expect(body.artifact_generation).toMatch(/^artifactgen-[0-9a-f]{64}$/);
    if (generation === null) generation = body.artifact_generation;
    expect(body.artifact_generation).toBe(generation);
    artifacts.push(...(body.artifacts ?? []));
    if (!body.next_cursor) return { artifacts, generation };
    expect(cursors.has(body.next_cursor)).toBe(false);
    cursors.add(body.next_cursor);
    after = body.next_cursor;
  }
  throw new Error('Artifact pagination exceeded the product Gate bound.');
}

async function loadQcMetrics(page: Page, runId: string) {
  const metrics: QcMetricRecord[] = [];
  const cursors = new Set<string>();
  let generation: string | null = null;
  let after: string | null = null;
  for (let pageNumber = 0; pageNumber < 20; pageNumber += 1) {
    const parameters = new URLSearchParams({ limit: '100' });
    if (generation !== null) parameters.set('generation', generation);
    if (after !== null) parameters.set('after', after);
    const response = await page.request.get(
      `/api/v1/runs/${runId}/qc-metrics?${parameters.toString()}`,
    );
    expect(response.status()).toBe(200);
    const body = (await response.json()) as {
      ok: boolean;
      run_id: string;
      qc_generation: string;
      qc_metrics?: QcMetricRecord[];
      next_cursor?: string | null;
    };
    expect(body.ok).toBe(true);
    expect(body.run_id).toBe(runId);
    expect(body.qc_generation).toMatch(/^qcgen-[0-9a-f]{64}$/);
    if (generation === null) generation = body.qc_generation;
    expect(body.qc_generation).toBe(generation);
    metrics.push(...(body.qc_metrics ?? []));
    if (!body.next_cursor) return { metrics, generation };
    expect(cursors.has(body.next_cursor)).toBe(false);
    cursors.add(body.next_cursor);
    after = body.next_cursor;
  }
  throw new Error('QC pagination exceeded the product Gate bound.');
}

async function persistPathFreeEvidence(
  testInfo: TestInfo,
  values: {
    workflowId: string;
    runId: string;
    snapshotId: string;
    payloadDigest: string;
    artifactGeneration: string;
    artifacts: ArtifactRecord[];
    qcGeneration: string;
    metrics: QcMetricRecord[];
    downloadedArtifactId: string;
    downloadedArtifactRevision: string;
    downloadedSha256: string;
    downloadedByteCount: number;
    workflowIdentity: PublicWorkflowIdentity;
  },
) {
  const evidenceDirectory =
    process.env.HELIXWEAVE_BULK_RNASEQ_PRODUCT_EVIDENCE_DIR;
  const testedHead = process.env.GITHUB_SHA;
  if (
    evidenceDirectory &&
    (typeof testedHead !== 'string' || !/^[0-9a-f]{40}$/.test(testedHead))
  ) {
    throw new Error(
      'Protected product evidence requires exact lowercase GITHUB_SHA identity.',
    );
  }
  const evidence = {
    ...(evidenceDirectory ? { tested_head: testedHead! } : {}),
    workflow_id: values.workflowId,
    workflow_schema_version: values.workflowIdentity.schemaVersion,
    workflow_availability_reason: values.workflowIdentity.availabilityReason,
    upstream: {
      name: values.workflowIdentity.upstreamName,
      version: values.workflowIdentity.upstreamVersion,
      revision: values.workflowIdentity.upstreamRevision,
    },
    run_id: values.runId,
    snapshot_id: values.snapshotId,
    payload_digest: values.payloadDigest,
    artifact_generation: values.artifactGeneration,
    artifact_count: values.artifacts.length,
    artifact_output_types: [
      ...new Set(values.artifacts.map((artifact) => artifact.output_type)),
    ].sort(),
    artifacts: values.artifacts
      .map((artifact) => ({
        artifact_id: artifact.artifact_id,
        output_type: artifact.output_type,
        revision: artifact.revision,
      }))
      .sort((left, right) => left.artifact_id.localeCompare(right.artifact_id)),
    qc_generation: values.qcGeneration,
    qc_metric_count: values.metrics.length,
    qc_metrics: values.metrics
      .map((metric) => ({
        metric_id: metric.metric_id,
        metric_key: metric.metric_key,
        value: metric.value,
        unit: metric.unit,
      }))
      .sort((left, right) =>
        `${left.metric_key}\u0000${left.metric_id}`.localeCompare(
          `${right.metric_key}\u0000${right.metric_id}`,
        ),
      ),
    download: {
      artifact_id: values.downloadedArtifactId,
      revision: values.downloadedArtifactRevision,
      sha256: values.downloadedSha256,
      byte_count: values.downloadedByteCount,
    },
  };
  const serialized = `${JSON.stringify(evidence, null, 2)}\n`;
  expect(serialized).not.toContain(runtimeRoot);
  expect(serialized).not.toMatch(
    /(?:relative_path|workspace|launch_dir|work_dir|docker|socket|environment|uri|filename)/i,
  );
  if (evidenceDirectory) {
    mkdirSync(evidenceDirectory, { recursive: true });
    writeFileSync(
      join(evidenceDirectory, 'bulk-rnaseq-product-evidence.json'),
      serialized,
      'utf8',
    );
  }
  await testInfo.attach('bulk-rnaseq-product-evidence', {
    body: Buffer.from(serialized, 'utf8'),
    contentType: 'application/json',
  });
}

function readExpectedSamples(fixture: RuntimeManifest) {
  const lines = readFileSync(fixture.bulkSamplesPath, 'utf8')
    .trimEnd()
    .split(/\r?\n/);
  const headers = lines[0]!.split('\t');
  return lines.slice(1).map((line) =>
    Object.fromEntries(
      line.split('\t').map((value, index) => [headers[index]!, value]),
    ),
  );
}

test('Bulk RNA-seq product path is fail-closed or executes by declared admission @desktop', async ({
  page,
}, testInfo) => {
  const fixture = manifest();
  if (fixture.bulkExpectedExecution === 'available') {
    test.setTimeout(45 * 60 * 1_000);
  }
  expect(fixture.bulkWorkflowId).toBe('bulk-rnaseq');
  const expectedSamples = readExpectedSamples(fixture);
  expect(new Set(expectedSamples.map((sample) => sample.layout))).toEqual(
    new Set(['SE', 'PE']),
  );
  expect(
    new Set(expectedSamples.map((sample) => sample.sample)).size,
  ).toBeLessThan(expectedSamples.length);
  expect(
    expectedSamples.every(
      (sample) => sample.strandedness && sample.platform === 'ILLUMINA',
    ),
  ).toBe(true);
  expect(fixture.bulkConfig).toHaveProperty('standard.trimming');
  expect(fixture.bulkConfig).toHaveProperty(
    'standard.ribosomal_rna_removal',
  );
  if (fixture.bulkExpectedExecution === 'not_configured') {
    expect(fixture.bulkConfig).toHaveProperty('standard.umi');
  }
  let runMutationRequests = 0;
  page.on('request', (request) => {
    if (
      request.method() === 'POST' &&
      (/\/api\/v1\/workflows\/bulk-rnaseq\/runs$/.test(request.url()) ||
        /\/api\/v1\/runs\/[^/]+\/start$/.test(request.url()))
    ) {
      runMutationRequests += 1;
    }
  });

  const newAnalysis = await openBulkAuthoring(page, fixture);
  const evidencePrefix =
    fixture.bulkExpectedExecution === 'available'
      ? 'bulk-available'
      : 'bulk-unavailable';
  await capture(page, testInfo, `${evidencePrefix}-detail-desktop-1440x900`);
  await newAnalysis.click();
  await expect(page).toHaveURL(
    `/workflows/${fixture.bulkWorkflowId}/new-run`,
  );
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();

  const configForm = page.getByLabel('Workflow config form');
  await expect(configForm).toContainText('Standard bulk RNA-seq parameters');
  await expect(configForm).toContainText(/trimming/i);
  await expect(configForm).toContainText(/ribosomal[_ ]rna[_ ]removal/i);
  await expect(configForm).toContainText(/umi/i);
  await capture(
    page,
    testInfo,
    `${evidencePrefix}-config-desktop-1440x900`,
  );
  await populateBulkDraft(page, fixture);

  const table = page.getByTestId('sample-desktop-table');
  await expect(table).toBeVisible();
  for (const [index, sample] of expectedSamples.entries()) {
    const rowNumber = index + 1;
    await expect(table.getByLabel(`Sample ${rowNumber} sample`)).toHaveValue(
      sample.sample,
    );
    await expect(table.getByLabel(`Sample ${rowNumber} layout`)).toHaveValue(
      sample.layout,
    );
    await expect(
      table.getByLabel(`Sample ${rowNumber} strandedness`),
    ).toHaveValue(sample.strandedness);
    await expect(table.getByLabel(`Sample ${rowNumber} platform`)).toHaveValue(
      'ILLUMINA',
    );
    await expect(
      table.getByLabel(`Sample ${rowNumber} platform`),
    ).toHaveAttribute('readonly', '');
  }
  await expectNoHorizontalOverflow(page);

  await page.getByRole('tab', { name: 'Options' }).click();
  await expect(
    page.getByRole('heading', { name: 'Adapter options' }),
  ).toBeVisible();
  await page.getByRole('tab', { name: 'Review' }).click();
  await expect(
    page.getByText(/structurally ready for backend validation/i),
  ).toBeVisible();
  const review = page.getByTestId('draft-review-json');
  await expect(review).toContainText('trimming');
  await expect(review).toContainText('ribosomal_rna_removal');
  if (fixture.bulkExpectedExecution === 'not_configured') {
    await expect(review).toContainText('umi');
  }
  for (const sample of new Set(expectedSamples.map((row) => row.sample))) {
    await expect(review).toContainText(sample);
  }
  if (fixture.bulkExpectedExecution === 'not_configured') {
    const availabilityNotice = page.getByTestId('execution-availability-notice');
    await expect(availabilityNotice).toContainText(
      'Execution runtime is not configured',
    );
    await captureElement(
      page,
      availabilityNotice,
      testInfo,
      'bulk-unavailable-availability-desktop',
    );
  } else {
    await expect(page.getByTestId('execution-availability-notice')).toHaveCount(0);
  }
  const createRun = page.getByRole('button', {
    name: 'Create run from validated inputs',
  });
  await expect(createRun).toBeDisabled();
  await expect(page.getByRole('button', { name: 'Start run' })).toHaveCount(0);
  await expectNoHorizontalOverflow(page);

  const validateRequestPromise = page.waitForRequest(
    (request) =>
      request.method() === 'POST' &&
      request
        .url()
        .endsWith(`/api/v1/workflows/${fixture.bulkWorkflowId}/validate`),
  );
  const validateResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response
        .url()
        .endsWith(`/api/v1/workflows/${fixture.bulkWorkflowId}/validate`),
  );
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  const validateRequest = await validateRequestPromise;
  const validateResponse = await validateResponsePromise;
  expect(validateResponse.status()).toBe(200);
  const requestPayload = validateRequest.postDataJSON() as {
    config: Record<string, unknown>;
    samples: Array<Record<string, string>>;
    options: Record<string, unknown>;
  };
  expect(requestPayload.config).toEqual(fixture.bulkConfig);
  expect(requestPayload.samples).toEqual(expectedSamples);
  expect(requestPayload.options).toEqual(fixture.bulkOptions);

  const validation = (await validateResponse.json()) as {
    ok: boolean;
    snapshot: null | { snapshot_id: string; payload_digest: string };
    canonical_payload?: unknown;
  };
  expect(validation.ok).toBe(true);
  expect(validation).not.toHaveProperty('canonical_payload');
  if (fixture.bulkExpectedExecution === 'not_configured') {
    expect(validation.snapshot).toBeNull();
    await expect(
      page.getByText(
        /No runnable snapshot was issued because execution is unavailable/i,
      ),
    ).toBeVisible();
    await expect(createRun).toBeDisabled();
    await expect(page.getByRole('button', { name: 'Start run' })).toHaveCount(0);
    expect(runMutationRequests).toBe(0);
    return;
  }

  expect(validation.snapshot?.snapshot_id).toMatch(/^vsnap_[0-9a-f]{32}$/);
  expect(validation.snapshot?.payload_digest).toMatch(/^[0-9a-f]{64}$/);
  expect(fixture.bulkRequiredArtifactOutputTypes.length).toBeGreaterThan(0);
  expect(fixture.bulkRequiredQcMetricKeys.length).toBeGreaterThan(0);
  await expect(
    page.getByText(/This exact draft can create one run/i),
  ).toBeVisible();
  await expect(createRun).toBeEnabled();

  const createResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response
        .url()
        .endsWith(`/api/v1/workflows/${fixture.bulkWorkflowId}/runs`),
  );
  const preflightResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      /\/api\/v1\/runs\/[^/]+\/preflight$/.test(response.url()),
  );
  await createRun.click();
  const createResponse = await createResponsePromise;
  expect(createResponse.status()).toBe(201);
  await expect(page).toHaveURL(/\/runs\/[^/]+$/);
  const runId = page.url().split('/').at(-1);
  expect(runId).toBeTruthy();
  const preflightResponse = await preflightResponsePromise;
  expect(preflightResponse.status()).toBe(202);
  await expect(page.getByTestId('run-status-badge')).toHaveText('planned', {
    timeout: 180_000,
  });

  const startResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response.url().endsWith(`/api/v1/runs/${runId}/start`),
  );
  const startRun = page.getByRole('button', { name: 'Start run' });
  await expect(startRun).toBeEnabled();
  await startRun.click();
  const startResponse = await startResponsePromise;
  expect(startResponse.status()).toBe(202);
  await expect(page.getByTestId('run-status-badge')).toHaveText(
    /queued|running|succeeded/,
    { timeout: 180_000 },
  );
  await expect(page.getByTestId('run-status-badge')).toHaveText('succeeded', {
    timeout: 40 * 60 * 1_000,
  });
  await expect(page).toHaveURL(`/runs/${runId}`);
  await captureElement(
    page,
    page.getByTestId('run-status-badge').locator('..'),
    testInfo,
    'bulk-available-run-succeeded-desktop',
  );

  await page.getByRole('tab', { name: 'Artifacts' }).click();
  await expect(
    page.getByRole('heading', { name: 'Indexed artifacts' }),
  ).toBeVisible({ timeout: 180_000 });
  const { artifacts, generation: artifactGeneration } = await loadArtifacts(
    page,
    runId!,
  );
  expect(artifacts.length).toBeGreaterThan(0);
  const outputTypes = new Set(artifacts.map((artifact) => artifact.output_type));
  for (const outputType of fixture.bulkRequiredArtifactOutputTypes) {
    expect(outputTypes.has(outputType), `missing artifact ${outputType}`).toBe(true);
  }
  const downloadedArtifact = artifacts.find((artifact) =>
    [
      'bulk_rnaseq.star.log_final',
      'bulk_rnaseq.salmon.meta_info',
      'bulk_rnaseq.trim.fastp.json',
    ].includes(artifact.output_type),
  );
  if (!downloadedArtifact) {
    throw new Error('No bounded machine-readable artifact is available to download.');
  }
  expect(downloadedArtifact.size_bytes).toBeLessThanOrEqual(2_097_152);
  await page
    .getByRole('button', { name: `Open artifact ${downloadedArtifact.name}` })
    .first()
    .click();
  const artifactInspector = page.getByTestId('artifact-inspector');
  await expect(artifactInspector).toContainText(downloadedArtifact.output_type);
  const downloadPromise = page.waitForEvent('download');
  const downloadResponsePromise = page.waitForResponse(
    (response) =>
      response.request().method() === 'GET' &&
      response
        .url()
        .includes(
          `/api/v1/runs/${runId}/artifacts/${downloadedArtifact.artifact_id}/download?`,
        ),
  );
  await artifactInspector
    .getByRole('button', { name: `Download ${downloadedArtifact.name}` })
    .click();
  const download = await downloadPromise;
  const downloadResponse = await downloadResponsePromise;
  expect(downloadResponse.status()).toBe(200);
  const downloadUrl = new URL(downloadResponse.url());
  expect(downloadUrl.searchParams.get('generation')).toBe(artifactGeneration);
  expect(downloadUrl.searchParams.get('revision')).toBe(
    downloadedArtifact.revision,
  );
  expect(await download.failure()).toBeNull();
  const downloadedPath = await download.path();
  expect(downloadedPath).not.toBeNull();
  const downloadedBytes = readFileSync(downloadedPath!);
  expect(downloadedBytes.byteLength).toBe(downloadedArtifact.size_bytes);
  const downloadedSha256 = createHash('sha256')
    .update(downloadedBytes)
    .digest('hex');
  await expect(page.getByText('Download prepared successfully.')).toBeVisible();
  await captureElement(
    page,
    page.getByRole('heading', { name: 'Indexed artifacts' }).locator('..'),
    testInfo,
    'bulk-available-artifacts-desktop',
  );

  await page.getByRole('tab', { name: 'QC' }).click();
  await expect(
    page.getByRole('heading', { name: 'Indexed QC metrics' }),
  ).toBeVisible({ timeout: 180_000 });
  const { metrics, generation: qcGeneration } = await loadQcMetrics(page, runId!);
  expect(metrics.length).toBeGreaterThan(0);
  const metricKeys = new Set(metrics.map((metric) => metric.metric_key));
  for (const metricKey of fixture.bulkRequiredQcMetricKeys) {
    expect(metricKeys.has(metricKey), `missing QC metric ${metricKey}`).toBe(true);
  }
  const firstMetric = metrics[0]!;
  const firstMetricRow = page
    .getByRole('table', { name: 'Indexed run QC metrics' })
    .getByRole('row')
    .filter({ hasText: firstMetric.metric_key });
  await expect(firstMetricRow.getByText(firstMetric.value, { exact: true })).toBeVisible();
  await expect(firstMetricRow.getByText(firstMetric.unit, { exact: true })).toBeVisible();
  await captureElement(
    page,
    page.getByTestId('qc-metric-list'),
    testInfo,
    'bulk-available-qc-desktop',
  );
  const workflowIdentity = await loadPublicWorkflowIdentity(
    page,
    fixture.bulkWorkflowId,
  );

  await persistPathFreeEvidence(testInfo, {
    workflowId: fixture.bulkWorkflowId,
    runId: runId!,
    snapshotId: validation.snapshot!.snapshot_id,
    payloadDigest: validation.snapshot!.payload_digest,
    artifactGeneration,
    artifacts,
    qcGeneration,
    metrics,
    downloadedArtifactId: downloadedArtifact.artifact_id,
    downloadedArtifactRevision: downloadedArtifact.revision,
    downloadedSha256,
    downloadedByteCount: downloadedBytes.byteLength,
    workflowIdentity,
  });
});

test('Bulk RNA-seq authoring remains operable without mobile overflow @mobile', async ({
  page,
}, testInfo) => {
  const fixture = manifest();
  const evidencePrefix =
    fixture.bulkExpectedExecution === 'available'
      ? 'bulk-available'
      : 'bulk-unavailable';
  const newAnalysis = await openBulkAuthoring(page, fixture);
  await capture(page, testInfo, `${evidencePrefix}-detail-mobile-390x844`);
  await newAnalysis.click();
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();
  await capture(page, testInfo, `${evidencePrefix}-config-mobile-390x844`);
  await populateBulkDraft(page, fixture);

  const sampleList = page.getByTestId('sample-mobile-list');
  await expect(sampleList).toBeVisible();
  await expect(page.getByTestId('sample-desktop-table')).toBeHidden();
  const expectedSamples = readExpectedSamples(fixture);
  for (const [index, sample] of expectedSamples.entries()) {
    const rowNumber = index + 1;
    await expect(sampleList.getByLabel(`Sample ${rowNumber} layout`)).toHaveValue(
      sample.layout,
    );
  }
  await expect(sampleList.getByLabel('Sample 1 platform')).toHaveValue('ILLUMINA');
  await expectNoHorizontalOverflow(page);

  await page.getByRole('tab', { name: 'Review' }).click();
  await expect(
    page.getByText(/structurally ready for backend validation/i),
  ).toBeVisible();
  await expect(
    page.getByRole('button', { name: 'Create run from validated inputs' }),
  ).toBeDisabled();
  if (fixture.bulkExpectedExecution === 'not_configured') {
    const availabilityNotice = page.getByTestId('execution-availability-notice');
    await expect(availabilityNotice).toContainText(
      'Execution runtime is not configured',
    );
    await captureElement(
      page,
      availabilityNotice,
      testInfo,
      'bulk-unavailable-availability-mobile',
    );
  } else {
    await expect(page.getByTestId('execution-availability-notice')).toHaveCount(0);
  }
  await expectNoHorizontalOverflow(page);
});
