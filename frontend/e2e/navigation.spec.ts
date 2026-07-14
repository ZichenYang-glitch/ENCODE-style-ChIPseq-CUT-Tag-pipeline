import { expect, test, type Page, type TestInfo } from '@playwright/test';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

interface RuntimeManifest {
  workflowId: string;
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
  const screenshot = await page.screenshot({ fullPage: false, scale: 'css' });
  const screenshotDirectory = process.env.ENCODE_PIPELINE_E2E_SCREENSHOT_DIR;
  if (screenshotDirectory) {
    mkdirSync(screenshotDirectory, { recursive: true });
    writeFileSync(join(screenshotDirectory, `${name}.png`), screenshot);
  }
  await testInfo.attach(name, { body: screenshot, contentType: 'image/png' });
}

test('global navigation exposes authoring before developer schemas on desktop @desktop', async ({
  page,
}, testInfo) => {
  const { workflowId } = manifest();
  await page.goto(`/workflows/${workflowId}`);

  await expect(page).toHaveTitle('HelixWeave');
  await expect(
    page.getByRole('link', { name: 'HelixWeave' }),
  ).toHaveAttribute('href', '/workflows');

  const primaryNavigation = page.getByRole('navigation', { name: 'Primary' });
  const workflowsLink = primaryNavigation.getByRole('link', {
    name: 'Workflows',
  });
  const runsLink = primaryNavigation.getByRole('link', { name: 'Runs' });
  const newAnalysisLink = primaryNavigation.getByRole('link', {
    name: 'New analysis',
  });
  await expect(primaryNavigation.getByRole('link')).toHaveText([
    'Workflows',
    'Runs',
    'New analysis',
  ]);
  await expect(workflowsLink).toHaveAttribute('aria-current', 'page');
  await expect(runsLink).not.toHaveAttribute('aria-current');
  await expect(runsLink).toHaveAttribute('href', '/runs');
  await expect(newAnalysisLink).toHaveAttribute(
    'href',
    `/workflows/${workflowId}/new-run`,
  );

  const authorInputs = page.getByRole('link', { name: 'Author inputs' });
  const schemaDetails = page.getByTestId('developer-schema-details');
  await expect(authorInputs).toBeVisible();
  await expect(authorInputs).toBeInViewport();
  await expect(schemaDetails).not.toHaveAttribute('open');
  const authorBox = await authorInputs.boundingBox();
  const schemaBox = await schemaDetails.boundingBox();
  expect(authorBox).not.toBeNull();
  expect(schemaBox).not.toBeNull();
  expect(authorBox!.y).toBeLessThan(schemaBox!.y);

  await capture(
    page,
    testInfo,
    'navigation-workflow-detail-desktop-1440x900',
    1440,
    900,
  );
  await capture(
    page,
    testInfo,
    'navigation-workflow-detail-tablet-1024x768',
    1024,
    768,
  );

  const schemaSummary = page.getByText('Developer schema', { exact: true });
  await schemaSummary.focus();
  await schemaSummary.press('Enter');
  await expect(schemaDetails).toHaveAttribute('open', '');
  await expect(page.getByText('Config schema', { exact: true })).toBeVisible();
  await expect(page.getByText('Sample schema', { exact: true })).toBeVisible();
  await expect(page.getByText('Options schema', { exact: true })).toBeVisible();
  await expectNoHorizontalOverflow(page);

  await newAnalysisLink.click();
  await expect(page).toHaveURL(`/workflows/${workflowId}/new-run`);
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();
  await expect(
    primaryNavigation.getByRole('link', { name: 'New analysis' }),
  ).toHaveAttribute('aria-current', 'page');
  await primaryNavigation.getByRole('link', { name: 'Workflows' }).click();
  await expect(page).toHaveURL('/workflows');
  await expect(page.getByRole('heading', { name: 'Workflows' })).toBeVisible();
});

test('navigation and workflow authoring remain operable without mobile overflow @mobile', async ({
  page,
}, testInfo) => {
  const { workflowId } = manifest();
  await page.goto(`/workflows/${workflowId}`);

  const primaryNavigation = page.getByRole('navigation', { name: 'Primary' });
  await expect(primaryNavigation).toBeVisible();
  await expect(
    primaryNavigation.getByRole('link', { name: 'Workflows' }),
  ).toHaveAttribute('aria-current', 'page');
  const authorInputs = page.getByRole('link', { name: 'Author inputs' });
  await expect(authorInputs).toBeVisible();
  await expect(authorInputs).toBeInViewport();
  await expect(page.getByTestId('developer-schema-details')).not.toHaveAttribute(
    'open',
  );

  await capture(
    page,
    testInfo,
    'navigation-workflow-detail-mobile-390x844',
    390,
    844,
  );
  await capture(
    page,
    testInfo,
    'navigation-workflow-detail-mobile-360x800',
    360,
    800,
  );

  await primaryNavigation.getByRole('link', { name: 'New analysis' }).click();
  await expect(page).toHaveURL(`/workflows/${workflowId}/new-run`);
  await expect(
    page.getByRole('heading', { name: 'Input workbench' }),
  ).toBeVisible();
  await expectNoHorizontalOverflow(page);
  await primaryNavigation.getByRole('link', { name: 'Workflows' }).click();
  await expect(page).toHaveURL('/workflows');
  await expectNoHorizontalOverflow(page);
});
