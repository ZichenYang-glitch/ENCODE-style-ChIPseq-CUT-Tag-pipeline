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

async function signIn(page, username, password) {
  await page.goto('/login');
  await page.locator('#login-username').fill(username);
  await page.locator('#login-password').fill(password);
  await page.getByRole('button', { name: 'Sign in' }).click();
  await page.waitForURL(/\/workflows$/);
}

async function signOut(page, username) {
  await page.getByRole('button', { name: username }).click();
  await page.getByRole('button', { name: 'Sign out' }).click();
  await page.waitForURL(/\/login$/);
}

async function openHistory(page, manifest) {
  await page.goto('/catalog');
  const search = page.getByRole('form', { name: 'Find experiments' });
  await search.getByLabel('Project').selectOption(manifest.project_id);
  await search.getByLabel('Assay').selectOption([manifest.assay_term_id]);
  await search.getByLabel('Cell type').selectOption([manifest.cell_type_term_id]);
  await search.getByRole('button', { name: 'Find experiments' }).click();
  await page
    .getByRole('list', { name: 'Experiments' })
    .getByRole('button', { name: manifest.current_name })
    .click();
}

function revisionSelect(page) {
  return page.getByRole('combobox', { name: 'Revision', exact: true });
}

async function revisionNumbers(select) {
  const labels = await select.locator('option').allTextContents();
  return labels.map((label) => {
    const match = label.match(/(?:当前|旧版)\s+(\d+)/);
    if (!match) throw new Error(`unparsed revision option: ${label}`);
    return Number(match[1]);
  });
}

function sameNumbers(actual, expected) {
  if (actual.length !== expected.length) return false;
  const seen = new Set(actual);
  return expected.every((number) => seen.has(number));
}

const api = spawn(python, [serverPath, '--long-history'], {
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
  const report = {
    member_opened_cross_page_revision: false,
    member_reachable_revisions: 0,
    admin_saw_revision_24: false,
    admin_created_revision_25: false,
    concurrent_conflict: false,
  };

  await signIn(page, manifest.member_username, manifest.member_password);
  await openHistory(page, manifest);
  const detail = page.getByRole('region', { name: 'Experiment revision' });
  await detail.getByRole('heading', { name: 'History 24', exact: true }).waitFor();
  await detail.locator('p', { hasText: '当前 24' }).waitFor();
  await revisionSelect(page)
    .locator('option', { hasText: /^旧版 1:/ })
    .waitFor({ state: 'attached' });
  const firstPage = await revisionNumbers(revisionSelect(page));
  const firstExpected = [...Array(20).keys()].map((index) => index + 1).concat(24);
  if (!sameNumbers(firstPage, firstExpected)) {
    throw new Error(`first revision page was ${firstPage.join(',')}`);
  }
  await page.getByRole('button', { name: 'Next revision page' }).click();
  await page.getByText('Revision page 2').waitFor();
  await revisionSelect(page)
    .locator('option', { hasText: /^旧版 21:/ })
    .waitFor({ state: 'attached' });
  const secondPage = await revisionNumbers(revisionSelect(page));
  if (!sameNumbers(secondPage, [21, 22, 23, 24])) {
    throw new Error(`second revision page was ${secondPage.join(',')}`);
  }
  report.member_reachable_revisions = new Set([...firstPage, ...secondPage]).size;
  if (report.member_reachable_revisions !== 24) {
    throw new Error(`expected revisions 1-24 across pages, saw ${report.member_reachable_revisions}`);
  }
  await page.getByRole('button', { name: 'Previous revision page' }).click();
  await page.getByText('Revision page 1').waitFor();
  await revisionSelect(page).selectOption(manifest.older_revision_id);
  await detail.getByRole('heading', { name: 'History 1', exact: true }).waitFor();
  await detail.locator('p', { hasText: '旧版 1' }).waitFor();
  report.member_opened_cross_page_revision = true;

  await signOut(page, manifest.member_username);
  await signIn(page, manifest.admin_username, manifest.admin_password);
  await openHistory(page, manifest);
  await detail.locator('p', { hasText: '当前 24' }).waitFor();
  const revise = page.getByRole('form', { name: 'Revise experiment classification' });
  const displayName = revise.getByLabel('Display name');
  if ((await displayName.inputValue()) !== 'History 24') {
    throw new Error(`revise form started at ${await displayName.inputValue()}`);
  }
  report.admin_saw_revision_24 = true;

  await displayName.fill('History 25');
  const created = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response.url().includes(`/experiments/${manifest.experiment_id}/revisions`),
  );
  await revise.getByRole('button', { name: 'Revise classification' }).click();
  const createdResponse = await created;
  const createdBody = createdResponse.request().postDataJSON();
  if (createdBody.expected_revision_number !== 24 || createdResponse.status() !== 200) {
    throw new Error(
      `expected revision 24 to create 25, got ${createdResponse.status()} ${JSON.stringify(createdBody)}`,
    );
  }
  await page.getByText('Revised classification to History 25.').waitFor();
  await detail.locator('p', { hasText: '当前 25' }).waitFor();
  report.admin_created_revision_25 = true;

  const cookies = await page.context().cookies();
  const csrf = cookies.find((cookie) => cookie.name === 'helixweave_csrf')?.value;
  if (!csrf) throw new Error('missing CSRF cookie');
  const winner = await page.request.post(
    `${origin}/api/v1/catalog/experiments/${manifest.experiment_id}/revisions`,
    {
      data: {
        expected_revision_number: 25,
        display_name: 'History 26',
        assay: manifest.write.assay,
        target: manifest.write.target,
        batch: manifest.write.batch,
        samples: [
          {
            sample_revision_id: manifest.write.sample_revision_id,
            replicate_role: manifest.write.replicate_role,
          },
        ],
      },
      headers: {
        'Content-Type': 'application/json',
        'Idempotency-Key': 'history-conflict-26',
        'X-CSRF-Token': csrf,
      },
    },
  );
  if (winner.status() !== 200) {
    throw new Error(`concurrent revise failed ${winner.status()} ${await winner.text()}`);
  }

  await displayName.fill('History loser');
  const losing = page.waitForResponse(
    (response) =>
      response.request().method() === 'POST' &&
      response.url().includes(`/experiments/${manifest.experiment_id}/revisions`) &&
      (response.request().postData() ?? '').includes('History loser'),
  );
  await revise.getByRole('button', { name: 'Revise classification' }).click();
  const losingResponse = await losing;
  const losingPayload = await losingResponse.json();
  if (
    losingResponse.status() !== 409 ||
    losingPayload.issues?.[0]?.code !== 'catalog_revision_conflict'
  ) {
    throw new Error(
      `expected catalog_revision_conflict, got ${losingResponse.status()} ${JSON.stringify(losingPayload)}`,
    );
  }
  await page.getByRole('alert').filter({ hasText: '冲突。这次提交没有写入。' }).waitFor();
  await detail.getByRole('heading', { name: 'History 26', exact: true }).waitFor();
  await detail.locator('p', { hasText: '当前 26' }).waitFor();
  if ((await page.getByText('Revised classification to History loser.').count()) !== 0) {
    throw new Error('losing revision was presented as saved');
  }
  report.concurrent_conflict = true;

  await browser.close();
  console.log(JSON.stringify(report, null, 2));
} finally {
  if (vite) vite.kill('SIGTERM');
  api.kill('SIGTERM');
}
