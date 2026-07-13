import { expect, test, type Page, type TestInfo } from '@playwright/test';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

interface RuntimeManifest {
  workflowId: string;
  samplesPath: string;
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

async function capture(
  page: Page,
  testInfo: TestInfo,
  name: string,
  width: number,
  height: number,
) {
  await page.setViewportSize({ width, height });
  await expectNoHorizontalOverflow(page);
  const screenshot = await page.screenshot({ fullPage: true });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `${name}.png`), screenshot);
  }
  await testInfo.attach(name, { body: screenshot, contentType: 'image/png' });
}

async function replaceCodeMirror(editor: ReturnType<Page['getByRole']>, text: string) {
  await editor.click();
  await editor.press('Control+A');
  await editor.press('Backspace');
  await editor.fill(text);
}

test('schema workbench preserves one draft across form YAML samples review and history @desktop', async ({
  page,
}, testInfo) => {
  const runtime = manifest();
  await page.goto(`/workflows/${runtime.workflowId}/new-run`);
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();
  await expect(page.getByText(/Draft only/)).toBeVisible();

  await page.getByRole('button', { name: 'YAML mode' }).click();
  const yamlEditor = page.getByRole('textbox', { name: 'Advanced config YAML' });
  await expect(yamlEditor).toBeVisible();
  await expect(yamlEditor).toContainText('outdir: results');
  await replaceCodeMirror(
    yamlEditor,
    'advanced_browser_field: retained\noutdir: results\nthreads: 8\n',
  );
  await page.getByRole('button', { name: 'Form mode' }).click();
  await page.getByRole('spinbutton', { name: /threads/i }).fill('12');
  await page.getByRole('button', { name: 'YAML mode' }).click();
  await expect(yamlEditor).toContainText('advanced_browser_field: retained');
  await expect(yamlEditor).toContainText('threads: 12');

  await page.getByRole('tab', { name: 'Samples' }).click();
  await page.getByLabel('Import samples TSV').setInputFiles(runtime.samplesPath);
  const desktopTable = page.getByTestId('sample-desktop-table');
  await expect(desktopTable).toBeVisible();
  const sampleInput = desktopTable.getByLabel('Sample 1 sample');
  await expect(sampleInput).not.toHaveValue('');
  await sampleInput.fill('browser-edited-sample');
  await capture(
    page,
    testInfo,
    'input-workbench-samples-desktop-1440x900',
    1440,
    900,
  );
  await capture(
    page,
    testInfo,
    'input-workbench-samples-tablet-1024x768',
    1024,
    768,
  );

  await page.getByRole('tab', { name: 'Options' }).click();
  await expect(page.getByLabel(/strict_inputs/i)).toBeVisible();
  await page.getByRole('tab', { name: 'Review' }).click();
  await expect(page).toHaveURL(/step=review/);
  const review = page.getByTestId('draft-review-json');
  await expect(review).toContainText('advanced_browser_field');
  await expect(review).toContainText('browser-edited-sample');
  await expect(page.getByText(/structurally ready for backend validation/i)).toBeVisible();
  await expect(
    page.getByRole('button', { name: /Create run|Start run|Validate/i }),
  ).toHaveCount(0);

  await page.goBack();
  await expect(page.getByRole('tab', { name: 'Options' })).toHaveAttribute(
    'data-state',
    'active',
  );
  await page.goBack();
  await expect(page.getByRole('tab', { name: 'Samples' })).toHaveAttribute(
    'data-state',
    'active',
  );
  await expect(desktopTable.getByLabel('Sample 1 sample')).toHaveValue(
    'browser-edited-sample',
  );
  await page.goForward();
  await page.goForward();
  await expect(review).toContainText('browser-edited-sample');

  await capture(page, testInfo, 'input-workbench-desktop-1440x900', 1440, 900);
  await capture(page, testInfo, 'input-workbench-tablet-1024x768', 1024, 768);

  await page.reload();
  await expect(page).toHaveURL(/step=review/);
  await expect(page.getByText(/Samples do not satisfy the adapter row schema/)).toBeVisible();
  await expect(review).not.toContainText('browser-edited-sample');
});

test('schema workbench sample records and YAML remain operable on mobile @mobile', async ({
  page,
}, testInfo) => {
  const runtime = manifest();
  await page.goto(
    `/workflows/${runtime.workflowId}/new-run?step=samples`,
  );
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();
  await page.getByLabel('Import samples TSV').setInputFiles(runtime.samplesPath);

  const mobileList = page.getByTestId('sample-mobile-list');
  await expect(mobileList).toBeVisible();
  await expect(page.getByTestId('sample-desktop-table')).toBeHidden();
  const sampleInput = mobileList.getByLabel('Sample 1 sample');
  await expect(sampleInput).not.toHaveValue('');
  await sampleInput.fill('mobile-edited-sample');
  await expect(sampleInput).toHaveValue('mobile-edited-sample');
  await capture(
    page,
    testInfo,
    'input-workbench-samples-mobile-390x844',
    390,
    844,
  );
  await capture(
    page,
    testInfo,
    'input-workbench-samples-mobile-360x800',
    360,
    800,
  );

  await page.getByRole('tab', { name: 'Config' }).click();
  await page.getByRole('button', { name: 'YAML mode' }).click();
  const yamlEditor = page.getByRole('textbox', { name: 'Advanced config YAML' });
  await expect(yamlEditor).toBeVisible();
  await replaceCodeMirror(
    yamlEditor,
    'outdir: results\nthreads: 6\nmobile_field: true\n',
  );
  await expect(yamlEditor).toContainText('mobile_field: true');

  await expect(page.getByRole('tab', { name: 'Review' })).toBeEnabled();
  await page.getByRole('tab', { name: 'Review' }).click();
  await expect(page.getByTestId('draft-review-json')).toContainText(
    'mobile-edited-sample',
  );
  await capture(page, testInfo, 'input-workbench-mobile-390x844', 390, 844);
  await capture(page, testInfo, 'input-workbench-mobile-360x800', 360, 800);
});
