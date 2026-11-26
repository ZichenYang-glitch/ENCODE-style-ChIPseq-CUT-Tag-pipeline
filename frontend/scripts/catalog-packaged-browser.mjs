import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { chromium } from '@playwright/test';

const repoRoot = fileURLToPath(new URL('../..', import.meta.url));
const serverPath = fileURLToPath(
  new URL('../../test/browser/catalog_page_server.py', import.meta.url),
);
const python = process.env.HELIXWEAVE_PYTHON ?? 'python3';

for (const key of [
  'HTTP_PROXY',
  'HTTPS_PROXY',
  'ALL_PROXY',
  'http_proxy',
  'https_proxy',
  'all_proxy',
]) {
  delete process.env[key];
}
process.env.NO_PROXY = '*';
process.env.no_proxy = '*';

function waitForReady(child) {
  return new Promise((resolve, reject) => {
    let settled = false;
    let buffer = '';
    const finish = (error, value) => {
      if (settled) return;
      settled = true;
      if (error) reject(error);
      else resolve(value);
    };
    child.on('exit', (code) => {
      finish(new Error(`catalog API exited before it was ready (${code})\n${buffer}`));
    });
    child.stderr.on('data', (chunk) => {
      process.stderr.write(chunk);
    });
    child.stdout.on('data', (chunk) => {
      buffer += chunk.toString();
      const line = buffer
        .split('\n')
        .find((item) => item.startsWith('CATALOG_PAGE_READY '));
      if (!line) return;
      try {
        finish(null, JSON.parse(line.slice('CATALOG_PAGE_READY '.length)));
      } catch (error) {
        finish(error instanceof Error ? error : new Error(String(error)));
      }
    });
  });
}

const api = spawn(python, [serverPath, '--long-history', '--packaged'], {
  cwd: repoRoot,
  env: {
    ...process.env,
    PYTHONPATH: fileURLToPath(new URL('../../src', import.meta.url)),
  },
  stdio: ['ignore', 'pipe', 'pipe'],
});

try {
  const manifest = await waitForReady(api);
  const browser = await chromium.launch({
    headless: true,
    args: ['--no-proxy-server'],
  });
  const page = await browser.newPage({ baseURL: manifest.api_base, locale: 'zh-CN' });
  page.setDefaultTimeout(20_000);

  const catalogDocument = page.waitForResponse(
    (response) =>
      response.request().method() === 'GET' &&
      new URL(response.url()).pathname === '/catalog',
  );
  await page.goto('/');
  await page.locator('#login-username').fill(manifest.admin_username);
  await page.locator('#login-password').fill(manifest.admin_password);
  await page.getByRole('button', { name: 'Sign in' }).click();
  await page.waitForURL(/\/workflows$/);
  await page.goto('/catalog');
  const catalogResponse = await catalogDocument;
  if (catalogResponse.status() !== 200) {
    throw new Error(`packaged /catalog returned ${catalogResponse.status()}`);
  }
  const contentType = catalogResponse.headers()['content-type'] ?? '';
  if (!contentType.includes('text/html')) {
    throw new Error(`packaged /catalog content type was ${contentType}`);
  }

  const search = page.getByRole('form', { name: 'Find experiments' });
  await search.getByLabel('Project').selectOption(manifest.project_id);
  await search.getByLabel('Assay').selectOption([manifest.assay_term_id]);
  await search.getByLabel('Cell type').selectOption([manifest.cell_type_term_id]);
  await search.getByRole('button', { name: 'Find experiments' }).click();
  await page
    .getByRole('list', { name: 'Experiments' })
    .getByRole('button', { name: 'History 24' })
    .click();
  const detail = page.getByRole('region', { name: 'Experiment revision' });
  await detail.getByRole('heading', { name: 'History 24', exact: true }).waitFor();
  await detail.locator('p', { hasText: '当前 24' }).waitFor();
  await page.getByRole('button', { name: 'Next revision page' }).click();
  await page.getByText('Revision page 2').waitFor();
  await detail.locator('option', { hasText: /^旧版 21:/ }).waitFor({ state: 'attached' });

  await browser.close();
  console.log(
    JSON.stringify(
      {
        packaged_catalog_document: true,
        opened_current_24: true,
        opened_revision_page_2: true,
      },
      null,
      2,
    ),
  );
} finally {
  api.kill('SIGTERM');
}
