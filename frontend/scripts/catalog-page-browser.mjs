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

async function selectOptionContaining(locator, text) {
  const labels = await locator.locator('option').allTextContents();
  const match = labels.find((label) => label.includes(text));
  if (match === undefined) {
    throw new Error(`no option containing ${text}: ${labels.join(' | ')}`);
  }
  await locator.selectOption({ label: match });
}

async function findExperiments(page, manifest) {
  await page.goto('/catalog');
  await page.getByRole('heading', { name: 'Experiment catalog' }).waitFor();
  const search = page.getByRole('form', { name: 'Find experiments' });
  await search.getByLabel('Project').selectOption(manifest.project_id);
  await search.getByLabel('Assay').selectOption([manifest.assay_term_id]);
  await search.getByLabel('Cell type').selectOption([manifest.cell_type_term_id]);
  await search.getByRole('button', { name: 'Find experiments' }).click();
  await page
    .getByRole('list', { name: 'Experiments' })
    .getByRole('button', { name: manifest.current_name })
    .waitFor();
}

const api = spawn(python, [serverPath], {
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
      env: {
        ...process.env,
        VITE_API_PROXY_TARGET: manifest.api_base,
      },
      stdio: ['ignore', 'pipe', 'pipe'],
    },
  );
  vite.stderr.on('data', (chunk) => process.stderr.write(chunk));
  await waitForHttp(`http://127.0.0.1:${vitePort}/login`);

  const browser = await chromium.launch({
    headless: true,
    args: ['--no-proxy-server'],
  });
  const page = await browser.newPage({
    locale: 'zh-CN',
    baseURL: `http://127.0.0.1:${vitePort}`,
  });
  page.setDefaultTimeout(20_000);
  const report = {
    api_base: manifest.api_base,
    member_maintenance_visible: null,
    admin_maintenance_visible: null,
    member_unbound: false,
    member_old_revision: false,
    member_second_page: false,
    excluded_hidden: false,
    admin_revised: false,
    admin_recorded_experiment: false,
  };

  await signIn(page, manifest.member_username, manifest.member_password);
  await findExperiments(page, manifest);
  const experiments = page.getByRole('list', { name: 'Experiments' });
  if ((await experiments.getByRole('button').count()) !== 20) {
    throw new Error(`expected 20 experiments on the first page`);
  }
  report.member_maintenance_visible = (await page
    .getByRole('form', { name: 'Record an experiment' })
    .count()) > 0;
  if (report.member_maintenance_visible) {
    throw new Error('member can see the experiment classification form');
  }

  await experiments
    .getByRole('button', { name: manifest.current_name })
    .click();
  const trace = page.getByRole('region', { name: 'Experiment trace' });
  await trace.getByRole('heading', { name: '未绑定' }).waitFor();
  const traceText = await trace.innerText();
  if (!traceText.includes('experiment_snapshot_binding_absent')) {
    throw new Error(traceText);
  }
  if (!traceText.includes('也不把运行样本写成文件贡献者')) {
    throw new Error(traceText);
  }
  report.member_unbound = true;
  const scope = page.getByRole('region', { name: 'Run input scope' });
  const scopeText = await scope.innerText();
  if (!scopeText.includes('运行输入范围') || !scopeText.includes('不是文件贡献者')) {
    throw new Error(scopeText);
  }
  if (!scopeText.includes('这条样本修订没有运行结果')) {
    throw new Error(scopeText);
  }

  await page
    .getByRole('combobox', { name: 'Revision', exact: true })
    .selectOption(manifest.older_revision_id);
  await page.getByRole('heading', { name: manifest.older_name }).waitFor();
  report.member_old_revision = true;

  await page.getByRole('button', { name: 'Next page' }).click();
  await experiments.getByRole('button', { name: manifest.page_two_name }).waitFor();
  await page.getByRole('heading', { name: manifest.older_name }).waitFor();
  report.opened_experiment_kept_after_paging = true;
  report.excluded_hidden =
    (await experiments.getByRole('button', { name: manifest.excluded_name }).count()) ===
    0;
  report.member_second_page =
    report.excluded_hidden &&
    (await experiments.getByRole('button').count()) === 1;
  if (!report.member_second_page) {
    throw new Error(await experiments.innerText());
  }

  await signOut(page, manifest.member_username);
  await signIn(page, manifest.admin_username, manifest.admin_password);
  await findExperiments(page, manifest);
  report.admin_maintenance_visible =
    (await page.getByRole('form', { name: 'Record an experiment' }).count()) > 0;
  if (!report.admin_maintenance_visible) {
    throw new Error('administrator cannot see the classification form');
  }

  const termForm = page.getByRole('form', { name: 'Record a term' });
  await termForm.getByLabel('Term type').selectOption('tag');
  await termForm.getByLabel('Display name').fill('browser-tag');
  await termForm.getByRole('button', { name: 'Record term' }).click();
  await page.getByText('Recorded term browser-tag.').waitFor();

  await experiments.getByRole('button', { name: manifest.current_name }).click();
  const revise = page.getByRole('form', { name: 'Revise experiment classification' });
  const revisedName = revise.getByLabel('Display name');
  await revisedName.fill('Treg CUT&Tag v3');
  if ((await revisedName.inputValue()) !== 'Treg CUT&Tag v3') {
    throw new Error(`revise display name is ${await revisedName.inputValue()}`);
  }
  await revise.getByRole('button', { name: 'Revise classification' }).click();
  try {
    await page.getByText('Revised classification to Treg CUT&Tag v3.').waitFor();
  } catch (error) {
    const alerts = await page.getByRole('alert').allTextContents();
    const detail = alerts.filter((item) => item.trim() !== '').join('\n');
    throw new Error(detail || (error instanceof Error ? error.message : String(error)));
  }
  report.admin_revised = true;

  const record = page.getByRole('form', { name: 'Record an experiment' });
  await record.getByLabel('Stable key').fill('browser-entered');
  await record.getByLabel('Display name').fill('Browser entered experiment');
  await record.getByLabel('Assay term').selectOption({ label: 'CUT&Tag' });
  await record.getByLabel('Target term').selectOption({ label: 'H3K27me3' });
  await record.getByLabel('Batch term').selectOption({ label: 'batch-1' });
  await selectOptionContaining(record.getByLabel('Sample revision'), 'treg');
  await record.getByLabel('Replicate role').selectOption({ label: 'treatment' });
  await record.getByRole('button', { name: 'Record experiment' }).click();
  await page.getByText('Recorded experiment Browser entered experiment.').waitFor();
  report.admin_recorded_experiment = true;

  await page
    .getByRole('form', { name: 'Find experiments' })
    .getByRole('button', { name: 'Find experiments' })
    .click();
  await page.getByRole('button', { name: 'Next page' }).click();
  await experiments.getByRole('button', { name: 'Browser entered experiment' }).waitFor();
  await browser.close();
  console.log(JSON.stringify(report, null, 2));
} finally {
  if (vite) vite.kill('SIGTERM');
  api.kill('SIGTERM');
}
