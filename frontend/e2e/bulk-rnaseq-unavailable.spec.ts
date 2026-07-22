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
  bulkRequiredSampleIds: string[];
  bulkRequiredArtifactSampleOutputTypes: [string, string][];
  bulkRequiredQcSampleMetricKeys: [string, string][];
  bulkRequiredQcSampleMetricValues: [string, string, string][];
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
  metadata: {
    sample_id?: string | null;
  };
}

interface QcMetricRecord {
  metric_id: string;
  metric_key: string;
  display_name: string;
  value: string;
  unit: string;
  sample_id: string | null;
  source_artifact_id: string;
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
    requiredSampleIds: string[];
    requiredArtifactSampleOutputTypes: [string, string][];
    qcGeneration: string;
    metrics: QcMetricRecord[];
    requiredQcSampleMetricKeys: [string, string][];
    requiredQcSampleMetricValues: [string, string, string][];
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
    fixture_oracle: {
      required_sample_ids: [...values.requiredSampleIds].sort(),
      artifact_sample_output_types: [...values.requiredArtifactSampleOutputTypes]
        .map(([sampleId, outputType]) => ({
          sample_id: sampleId,
          output_type: outputType,
        }))
        .sort((left, right) =>
          `${left.sample_id}\u0000${left.output_type}`.localeCompare(
            `${right.sample_id}\u0000${right.output_type}`,
          ),
        ),
      qc_sample_metric_keys: [...values.requiredQcSampleMetricKeys]
        .map(([sampleId, metricKey]) => ({
          sample_id: sampleId,
          metric_key: metricKey,
        }))
        .sort((left, right) =>
          `${left.sample_id}\u0000${left.metric_key}`.localeCompare(
            `${right.sample_id}\u0000${right.metric_key}`,
          ),
        ),
      qc_sample_metric_values: [...values.requiredQcSampleMetricValues]
        .map(([sampleId, metricKey, value]) => ({
          sample_id: sampleId,
          metric_key: metricKey,
          value,
        }))
        .sort((left, right) => {
          const leftKey = `${left.sample_id}\u0000${left.metric_key}\u0000${left.value}`;
          const rightKey = `${right.sample_id}\u0000${right.metric_key}\u0000${right.value}`;
          return leftKey.localeCompare(rightKey);
        }),
    },
    artifacts: values.artifacts
      .map((artifact) => ({
        artifact_id: artifact.artifact_id,
        output_type: artifact.output_type,
        revision: artifact.revision,
        sample_id: artifact.metadata.sample_id ?? null,
      }))
      .sort((left, right) => left.artifact_id.localeCompare(right.artifact_id)),
    qc_generation: values.qcGeneration,
    qc_metric_count: values.metrics.length,
    qc_metrics: values.metrics
      .map((metric) => ({
        metric_id: metric.metric_id,
        metric_key: metric.metric_key,
        display_name: metric.display_name,
        value: metric.value,
        unit: metric.unit,
        sample_id: metric.sample_id,
        source_artifact_id: metric.source_artifact_id,
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
  expect(fixture.bulkRequiredArtifactSampleOutputTypes.length).toBeGreaterThan(0);
  expect(fixture.bulkRequiredQcSampleMetricValues.length).toBeGreaterThan(0);
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
  const artifactIds = new Set(artifacts.map((artifact) => artifact.artifact_id));
  const outputTypes = new Set(artifacts.map((artifact) => artifact.output_type));
  for (const outputType of fixture.bulkRequiredArtifactOutputTypes) {
    expect(outputTypes.has(outputType), `missing artifact ${outputType}`).toBe(true);
  }
  const artifactSampleOutputTypes = new Set(
    artifacts.flatMap((artifact) =>
      artifact.metadata.sample_id
        ? [`${artifact.metadata.sample_id}\u0000${artifact.output_type}`]
        : [],
    ),
  );
  for (const [sampleId, outputType] of
    fixture.bulkRequiredArtifactSampleOutputTypes) {
    expect(
      artifactSampleOutputTypes.has(`${sampleId}\u0000${outputType}`),
      `missing artifact ${sampleId}/${outputType}`,
    ).toBe(true);
  }
  for (const sampleId of fixture.bulkRequiredSampleIds) {
    expect(
      artifacts.some((artifact) => artifact.metadata.sample_id === sampleId),
      `missing artifact ownership for ${sampleId}`,
    ).toBe(true);
  }
  const artifactTable = page.getByRole('table', {
    name: 'Indexed run artifacts',
  });
  for (let pageNumber = 0; pageNumber < 20; pageNumber += 1) {
    const loadMore = page.getByRole('button', { name: 'Load more artifacts' });
    if ((await loadMore.count()) === 0) break;
    const rowsBefore = await artifactTable.getByRole('row').count();
    await expect(loadMore).toBeEnabled();
    await loadMore.click();
    await expect
      .poll(() => artifactTable.getByRole('row').count())
      .toBeGreaterThan(rowsBefore);
  }
  await expect(
    page.getByRole('button', { name: 'Load more artifacts' }),
  ).toHaveCount(0);
  for (const [sampleId, outputType] of
    fixture.bulkRequiredArtifactSampleOutputTypes) {
    const row = artifactTable
      .getByRole('row')
      .filter({ hasText: outputType })
      .filter({ hasText: sampleId });
    await expect(
      row.first(),
      `artifact UI omitted ${sampleId}/${outputType}`,
    ).toBeVisible();
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
  const qcSampleMetricKeys = new Set(
    metrics.flatMap((metric) =>
      metric.sample_id
        ? [`${metric.sample_id}\u0000${metric.metric_key}`]
        : [],
    ),
  );
  const qcSampleMetricValues = new Set(
    metrics.flatMap((metric) =>
      metric.sample_id
        ? [`${metric.sample_id}\u0000${metric.metric_key}\u0000${metric.value}`]
        : [],
    ),
  );
  for (const [sampleId, metricKey] of fixture.bulkRequiredQcSampleMetricKeys) {
    expect(
      qcSampleMetricKeys.has(`${sampleId}\u0000${metricKey}`),
      `missing QC metric ${sampleId}/${metricKey}`,
    ).toBe(true);
  }
  for (const [sampleId, metricKey, value] of
    fixture.bulkRequiredQcSampleMetricValues) {
    expect(
      qcSampleMetricValues.has(`${sampleId}\u0000${metricKey}\u0000${value}`),
      `QC value drifted for ${sampleId}/${metricKey}`,
    ).toBe(true);
  }
  for (const metric of metrics) {
    expect(
      artifactIds.has(metric.source_artifact_id),
      `QC source artifact is missing for ${metric.metric_key}`,
    ).toBe(true);
  }
  for (const sampleId of fixture.bulkRequiredSampleIds) {
    expect(
      metrics.some((metric) => metric.sample_id === sampleId),
      `missing QC ownership for ${sampleId}`,
    ).toBe(true);
  }
  const qcTable = page.getByRole('table', {
    name: 'Indexed run QC metrics',
  });
  for (let pageNumber = 0; pageNumber < 20; pageNumber += 1) {
    const loadMore = page.getByRole('button', { name: 'Load more QC metrics' });
    if ((await loadMore.count()) === 0) break;
    const rowsBefore = await qcTable.getByRole('row').count();
    await expect(loadMore).toBeEnabled();
    await loadMore.click();
    await expect
      .poll(() => qcTable.getByRole('row').count())
      .toBeGreaterThan(rowsBefore);
  }
  await expect(
    page.getByRole('button', { name: 'Load more QC metrics' }),
  ).toHaveCount(0);
  for (const [sampleId, metricKey, value] of
    fixture.bulkRequiredQcSampleMetricValues) {
    const metric = metrics.find(
      (candidate) =>
        candidate.sample_id === sampleId &&
        candidate.metric_key === metricKey &&
        candidate.value === value,
    );
    expect(
      metric,
      `QC API omitted ${sampleId}/${metricKey}/${value}`,
    ).toBeDefined();
    const row = qcTable
      .getByRole('row')
      .filter({ hasText: metricKey })
      .filter({ hasText: sampleId });
    await expect(
      row.first(),
      `QC UI omitted ${sampleId}/${metricKey}`,
    ).toBeVisible();
    await expect(row.first().getByText(value, { exact: true })).toBeVisible();
    await expect(row.first()).toContainText(metric!.source_artifact_id);
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
    requiredSampleIds: fixture.bulkRequiredSampleIds,
    requiredArtifactSampleOutputTypes:
      fixture.bulkRequiredArtifactSampleOutputTypes,
    qcGeneration,
    metrics,
    requiredQcSampleMetricKeys: fixture.bulkRequiredQcSampleMetricKeys,
    requiredQcSampleMetricValues: fixture.bulkRequiredQcSampleMetricValues,
    downloadedArtifactId: downloadedArtifact.artifact_id,
    downloadedArtifactRevision: downloadedArtifact.revision,
    downloadedSha256,
    downloadedByteCount: downloadedBytes.byteLength,
    workflowIdentity,
  });

  await page.setViewportSize({ width: 390, height: 844 });
  await page.getByRole('tab', { name: 'Artifacts' }).click();
  await expect(artifactTable).toBeHidden();
  const artifactMobileList = page.getByTestId('artifact-list').locator('ul');
  await expect(artifactMobileList).toBeVisible();
  await expectNoHorizontalOverflow(page);
  const [mobileArtifactSampleId, mobileArtifactOutputType] =
    fixture.bulkRequiredArtifactSampleOutputTypes[0]!;
  const mobileArtifact = artifacts.find(
    (artifact) =>
      artifact.metadata.sample_id === mobileArtifactSampleId &&
      artifact.output_type === mobileArtifactOutputType,
  );
  expect(mobileArtifact).toBeDefined();
  const mobileArtifactCard = artifactMobileList
    .getByRole('listitem')
    .filter({ hasText: mobileArtifactSampleId })
    .filter({ hasText: mobileArtifactOutputType })
    .first();
  await expect(mobileArtifactCard).toBeVisible();
  const openMobileArtifact = mobileArtifactCard.getByRole('button', {
    name: `Open artifact ${mobileArtifact!.name}`,
  });
  await expect(openMobileArtifact).toBeEnabled();
  await captureElement(
    page,
    mobileArtifactCard,
    testInfo,
    'bulk-available-artifact-mobile-390x844',
  );
  await openMobileArtifact.click();
  await expect(page).toHaveURL(
    `/runs/${runId}?view=artifacts&artifact=${mobileArtifact!.artifact_id}`,
  );
  await expect(artifactInspector).toContainText(mobileArtifactOutputType);
  await expectNoHorizontalOverflow(page);
  await page.getByRole('tab', { name: 'QC' }).click();
  await expect(qcTable).toBeHidden();
  const qcMobileList = page.getByRole('list', {
    name: 'Indexed run QC metrics',
  });
  await expect(qcMobileList).toBeVisible();
  await expectNoHorizontalOverflow(page);
  const [mobileQcSampleId, mobileQcMetricKey, mobileQcValue] =
    fixture.bulkRequiredQcSampleMetricValues[0]!;
  const mobileQcMetric = metrics.find(
    (metric) =>
      metric.sample_id === mobileQcSampleId &&
      metric.metric_key === mobileQcMetricKey &&
      metric.value === mobileQcValue,
  );
  expect(mobileQcMetric).toBeDefined();
  const mobileQcSourceArtifact = artifacts.find(
    (artifact) =>
      artifact.artifact_id === mobileQcMetric!.source_artifact_id,
  );
  expect(mobileQcSourceArtifact).toBeDefined();
  const mobileQcCard = qcMobileList
    .getByRole('listitem')
    .filter({ hasText: mobileQcSampleId })
    .filter({ hasText: mobileQcMetricKey })
    .filter({ hasText: mobileQcValue })
    .filter({ hasText: mobileQcMetric!.source_artifact_id })
    .first();
  await expect(mobileQcCard).toBeVisible();
  const openMobileQcSource = mobileQcCard.getByRole('button', {
    name: `Open source artifact for ${mobileQcMetric!.display_name}`,
  });
  await expect(openMobileQcSource).toBeEnabled();
  await captureElement(
    page,
    mobileQcCard,
    testInfo,
    'bulk-available-qc-mobile-390x844',
  );
  await openMobileQcSource.click();
  await expect(page).toHaveURL(
    `/runs/${runId}?view=artifacts&artifact=${mobileQcMetric!.source_artifact_id}`,
  );
  await expect(artifactInspector).toContainText(
    mobileQcSourceArtifact!.output_type,
  );
  await expectNoHorizontalOverflow(page);
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
