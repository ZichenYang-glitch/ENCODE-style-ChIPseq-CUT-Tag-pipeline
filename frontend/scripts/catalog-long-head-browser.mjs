import { spawn } from 'node:child_process';
import { createServer } from 'node:net';
import { fileURLToPath } from 'node:url';
import { chromium } from '@playwright/test';

const repoRoot = fileURLToPath(new URL('../..', import.meta.url));
const frontendRoot = fileURLToPath(new URL('..', import.meta.url));
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

function freePort() {
  return new Promise((resolve, reject) => {
    const server = createServer();
    server.once('error', reject);
    server.listen(0, '127.0.0.1', () => {
      const address = server.address();
      const port = typeof address === 'object' && address ? address.port : 0;
      server.close(() => resolve(port));
    });
  });
}

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
      finish(
        new Error(`catalog API exited before it was ready (${code})\n${buffer}`),
      );
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

async function waitForHttp(url) {
  const deadline = Date.now() + 30_000;
  let lastError = `timed out waiting for ${url}`;
  while (Date.now() < deadline) {
    try {
      const response = await fetch(url);
      if (response.ok) return;
      lastError = `${url} returned ${response.status}`;
    } catch (error) {
      lastError = error instanceof Error ? error.message : String(error);
    }
    await new Promise((resolve) => setTimeout(resolve, 200));
  }
  throw new Error(lastError);
}

const api = spawn(python, [serverPath, '--head-revision', '401'], {
  cwd: repoRoot,
  env: {
    ...process.env,
    PYTHONPATH: fileURLToPath(new URL('../../src', import.meta.url)),
  },
  stdio: ['ignore', 'pipe', 'pipe'],
});
const vitePort = await freePort();
let vite;

try {
  const manifest = await waitForReady(api);
  const origin = `http://127.0.0.1:${vitePort}`;
  vite = spawn(
    process.execPath,
    [
      fileURLToPath(new URL('../node_modules/vite/bin/vite.js', import.meta.url)),
      '--host',
      '127.0.0.1',
      '--port',
      String(vitePort),
      '--strictPort',
    ],
    {
      cwd: frontendRoot,
      env: { ...process.env, VITE_API_PROXY_TARGET: manifest.api_base },
      stdio: ['ignore', 'pipe', 'pipe'],
    },
  );
  vite.stderr.on('data', (chunk) => process.stderr.write(chunk));
  await waitForHttp(`${origin}/login`);

  const browser = await chromium.launch({
    headless: true,
    args: ['--no-proxy-server'],
  });
  const page = await browser.newPage({ baseURL: origin });
  page.setDefaultTimeout(20_000);
  const historyPages = [];
  page.on('response', (response) => {
    const url = response.url();
    if (
      response.request().method() !== 'GET' ||
      !url.includes('/catalog/experiments/') ||
      !url.includes('/revisions')
    ) {
      return;
    }
    historyPages.push(response.json());
  });

  await page.goto('/login');
  await page.locator('#login-username').fill(manifest.admin_username);
  await page.locator('#login-password').fill(manifest.admin_password);
  await page.getByRole('button', { name: 'Sign in' }).click();
  await page.waitForURL(/\/workflows$/);
  await page.goto('/catalog');
  const search = page.getByRole('form', { name: 'Find experiments' });
  await search.getByLabel('Project').selectOption(manifest.project_id);
  await search.getByLabel('Assay').selectOption([manifest.assay_term_id]);
  await search.getByLabel('Cell type').selectOption([manifest.cell_type_term_id]);
  await search.getByRole('button', { name: 'Find experiments' }).click();
  await page
    .getByRole('list', { name: 'Experiments' })
    .getByRole('button', { name: 'History 401' })
    .click();

  const detail = page.getByRole('region', { name: 'Experiment revision' });
  await detail.getByRole('heading', { name: 'History 401', exact: true }).waitFor();
  await detail.locator('p', { hasText: '当前 401' }).waitFor();
  const revise = page.getByRole('form', { name: 'Revise experiment classification' });
  if ((await revise.getByLabel('Display name').inputValue()) !== 'History 401') {
    throw new Error(`form opened at ${await revise.getByLabel('Display name').inputValue()}`);
  }

  const revision401Page = page.waitForResponse(async (response) => {
    const url = response.url();
    if (
      response.request().method() !== 'GET' ||
      !url.includes('/catalog/experiments/') ||
      !url.includes('/revisions')
    ) {
      return false;
    }
    const body = await response.json();
    return (body.experiment_revisions ?? []).some(
      (item) => item.revision_number === 401,
    );
  });
  for (let step = 0; step < 20; step += 1) {
    await page.getByRole('button', { name: 'Next revision page' }).click();
  }
  await revision401Page;
  await page.getByText('Revision page 21').waitFor();
  const revisionLabels = await detail.getByLabel('Revision').locator('option').allTextContents();
  if (!revisionLabels.some((label) => label.includes('当前 401'))) {
    throw new Error(`revision 401 was not selectable: ${revisionLabels.join(' | ')}`);
  }
  const pages = await Promise.all(historyPages);
  const numbers = pages.map((body) =>
    (body.experiment_revisions ?? []).map((item) => item.revision_number),
  );
  const first = numbers[0] ?? [];
  if (Math.max(...first) > 20) {
    throw new Error(`first history page was not a single batch: ${first.join(',')}`);
  }
  if (!numbers.some((pageNumbers) => pageNumbers.includes(401))) {
    throw new Error(`revision 401 was not on a later history page: ${JSON.stringify(numbers)}`);
  }

  await revise.getByLabel('Display name').fill('History 402');
  const saved = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response.url().includes(`/experiments/${manifest.experiment_id}/revisions`),
  );
  await revise.getByRole('button', { name: 'Revise classification' }).click();
  const savedResponse = await saved;
  const savedBody = savedResponse.request().postDataJSON();
  if (savedResponse.status() !== 200 || savedBody.expected_revision_number !== 401) {
    throw new Error(
      `expected to save 402 from 401, got ${savedResponse.status()} ${JSON.stringify(savedBody)}`,
    );
  }
  await page.getByText('Revised classification to History 402.').waitFor();
  await detail.locator('p', { hasText: '当前 402' }).waitFor();
  await browser.close();
  console.log(
    JSON.stringify(
      {
        reopened_at_401: true,
        reached_revision_page_21: true,
        saved_402_from_401: true,
        history_batches: numbers.length,
      },
      null,
      2,
    ),
  );
} finally {
  if (vite) vite.kill('SIGTERM');
  api.kill('SIGTERM');
}
