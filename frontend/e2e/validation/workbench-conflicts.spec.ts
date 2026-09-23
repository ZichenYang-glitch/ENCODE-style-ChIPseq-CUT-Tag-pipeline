// Run against the validation-only real API fixture in test/browser/workbench_validation_runtime.py.
// No API response is synthesized. A delayed-response case only holds a real response.
import { expect, test, type Page, type BrowserContext } from '@playwright/test';
import { login } from '../auth';

// Reuse a real authenticated session within each viewport project; do not exhaust login throttling.
let sessionCookies: Awaited<ReturnType<BrowserContext['cookies']>>;
test.beforeAll(async ({ browser, baseURL }) => {
  const context = await browser.newContext({ baseURL });
  try {
    await login(await context.newPage());
    sessionCookies = await context.cookies();
  } finally {
    await context.close();
  }
});
test.beforeEach(async ({ context }) => { await context.addCookies(sessionCookies); });

async function configure(page: Page, config: object) {
  await page.getByRole('tab', { name: 'Config', exact: true }).click();
  await page.getByRole('button', { name: 'YAML mode' }).click();
  const editor = page.getByRole('textbox', { name: 'Advanced config YAML' });
  // CodeMirror owns selection; use its keyboard path before inserting replacement text.
  await editor.click();
  await editor.press('Control+A');
  await editor.press('Backspace');
  await editor.fill(JSON.stringify(config));
  await expect(editor).toHaveText(JSON.stringify(config));
  await page.getByRole('button', { name: 'Form mode' }).click();
}

async function start(page: Page, config: object) {
  await page.goto('/workflows/bulk-rnaseq/new-run');
  await expect(page.getByRole('combobox', { name: 'Reference profile', exact: true })).toContainText('Tiny validation reference');
  await configure(page, config);
  await page.getByRole('tab', { name: 'Samples', exact: true }).click();
  await page.getByLabel('Import samples TSV').setInputFiles({
    name: 'tiny.tsv', mimeType: 'text/tab-separated-values',
    buffer: Buffer.from('sample\tlibrary\tlane\tlayout\tfastq_1\tstrandedness\tplatform\nS1\tlib1\tL001\tSE\t/synthetic/reads.fastq.gz\treverse\tILLUMINA\n'),
  });
}

async function validate(page: Page) {
  await page.getByRole('tab', { name: 'Review', exact: true }).click();
  const response = page.waitForResponse(r => r.url().endsWith('/bulk-rnaseq/validate') && r.request().method() === 'POST');
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  const result = await response;
  expect(result.status()).toBe(200);
  return { response: await result.json(), payload: result.request().postDataJSON() };
}

async function noOverflow(page: Page) {
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= document.documentElement.clientWidth)).toBe(true);
}

test('QC conflict uses the real HTTP hint, preserves values, and recovers by keyboard @validation', async ({ page }, info) => {
  await start(page, { standard: { qc: { enabled: true, rseqc: true } }, advanced: { rseqc_modules: 'bam_stat' } });
  await page.getByRole('tab', { name: 'Config', exact: true }).click();
  await page.locator('#root_standard_qc_enabled').uncheck();
  const { response, payload } = await validate(page);
  expect(response.ok).toBe(false);
  expect(response.issues[0].code).toBe('BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT');
  expect(payload.config.advanced.rseqc_modules).toBe('bam_stat');
  const errors = page.getByTestId('validation-errors');
  await expect(errors).toContainText('enable both standard.qc.enabled and standard.qc.rseqc');
  await noOverflow(page);
  await info.attach('conflict', { body: await errors.screenshot(), contentType: 'image/png' });
  await page.getByRole('tab', { name: 'Options', exact: true }).click();
  const locate = page.getByRole('button', { name: 'Go to config.advanced.rseqc_modules', exact: true });
  await locate.focus();
  await page.keyboard.press('Enter');
  const target = page.locator('#root_advanced_rseqc_modules');
  await expect(target).toBeFocused();
  await expect(target).toHaveValue('bam_stat');
  await expect(page.locator('#root_advanced_rseqc_modules__error')).toContainText('Remove advanced.rseqc_modules');
  await page.getByRole('tab', { name: 'Review', exact: true }).click();
  expect(JSON.parse(await page.getByTestId('draft-review-json').innerText()).config).toEqual(payload.config);
  await page.getByRole('button', { name: 'Go to config.advanced.rseqc_modules' }).click();
  await page.locator('#root_standard_qc_enabled').check();
  await expect(errors).toHaveCount(0);
  // Reopening the master alone does not restore the dependent sub-flag.
  expect((await validate(page)).response.ok).toBe(false);
  await page.getByRole('button', { name: 'Go to config.advanced.rseqc_modules' }).click();
  await page.locator('#root_standard_qc_rseqc').focus();
  await page.keyboard.press('Space');
  await expect(errors).toHaveCount(0);
  const corrected = await validate(page);
  expect(corrected.response.ok).toBe(true);
  expect(corrected.response.snapshot).toBeNull();
  expect(corrected.payload.config.advanced).toEqual(payload.config.advanced);
  await expect(page.getByText(/Backend validation succeeded/)).toBeVisible();
  await noOverflow(page);
});

test('UMI output conflict locates the outputs section and leaves unrelated values intact @validation', async ({ page }, info) => {
  await start(page, { standard: { umi: { enabled: true, mode: 'read_name' }, outputs: { umi_intermediates: true } } });
  await page.getByRole('tab', { name: 'Config', exact: true }).click();
  await page.locator('#root_standard_umi_enabled').uncheck();
  const initial = await validate(page);
  expect(initial.response.issues[0].code).toBe('BULK_RNASEQ_OUTPUT_CONFLICT');
  expect(initial.payload.config.standard.umi).toEqual({ enabled: false });
  expect(initial.payload.config.standard.outputs.umi_intermediates).toBe(true);
  await expect(page.getByTestId('validation-errors')).toContainText('complete its required settings');
  await page.getByRole('button', { name: 'Go to config.standard.outputs' }).focus();
  await page.keyboard.press('Enter');
  await expect(page.locator('#root_standard_outputs > legend')).toBeFocused();
  await noOverflow(page);
  await info.attach('outputs', { body: await page.locator('#root_standard_outputs').screenshot(), contentType: 'image/png' });
  await page.locator('#root_standard_outputs_umi_intermediates').focus();
  await page.keyboard.press('Space');
  const corrected = await validate(page);
  expect(corrected.response.ok).toBe(true);
  const expected = structuredClone(initial.payload.config);
  expected.standard.outputs.umi_intermediates = false;
  expect(corrected.payload.config).toEqual(expected);
});

test('a delayed real conflict response cannot overwrite an edited draft @validation', async ({ page }) => {
  await start(page, { standard: { qc: { enabled: false } }, advanced: { deseq2_vst: false } });
  let release!: () => void;
  let received!: () => void;
  const hold = new Promise<void>(resolve => { release = resolve; });
  const ready = new Promise<void>(resolve => { received = resolve; });
  await page.route('**/bulk-rnaseq/validate', async route => {
    const response = await route.fetch();
    expect((await response.json()).issues[0].code).toBe('BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT');
    received();
    await hold;
    await route.fulfill({ response });
  }, { times: 1 });
  await page.getByRole('tab', { name: 'Review', exact: true }).click();
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  await ready;
  // Explicitly remove the optional key; false was already present and is not deletion.
  await configure(page, { standard: { qc: { enabled: false } }, advanced: {} });
  release();
  await expect(page.getByText(/Inputs changed while validation was running/)).toBeVisible();
  await expect(page.getByTestId('validation-errors')).toHaveCount(0);
  expect((await validate(page)).response.ok).toBe(true);
  await noOverflow(page);
});
