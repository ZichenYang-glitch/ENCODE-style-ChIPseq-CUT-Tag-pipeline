import { spawn } from 'node:child_process';
import { mkdtempSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { chromium } from '@playwright/test';

const root = new URL('../..', import.meta.url).pathname;
const python =
  process.env.HELIXWEAVE_PYTHON ??
  join(root, '.local/envs/ci-fast/bin/python');
const work = mkdtempSync(join(tmpdir(), 'washu-browser-'));
const apiPort = 8765;
const webPort = 4174;
const children = [];

function start(command, args, env = {}) {
  const child = spawn(command, args, {
    cwd: root,
    env: { ...process.env, ...env },
    stdio: ['ignore', 'pipe', 'pipe'],
    detached: true,
  });
  children.push(child);
  return child;
}

function waitFor(stream, pattern, timeoutMs) {
  return new Promise((resolve, reject) => {
    let text = '';
    const timer = setTimeout(() => {
      reject(new Error(`timed out waiting for ${pattern}\n${text}`));
    }, timeoutMs);
    const onData = (chunk) => {
      text += chunk.toString();
      const match = text.match(pattern);
      if (match) {
        clearTimeout(timer);
        resolve({ match, text });
      }
    };
    stream.on('data', onData);
  });
}

async function stopAll() {
  for (const child of children) {
    try {
      process.kill(-child.pid, 'SIGTERM');
    } catch {
      child.kill('SIGTERM');
    }
  }
}

try {
  const api = start(python, [
    'test/browser/washu_fixture_server.py',
    '--port',
    String(apiPort),
    '--database',
    join(work, 'platform.db'),
    '--workspace',
    join(work, 'workspaces'),
  ]);
  api.stderr.on('data', (chunk) => process.stderr.write(chunk));
  const ready = await waitFor(api.stdout, /WASHU_FIXTURE run_id=(\S+)/, 60_000);
  const runId = ready.match[1];
  const web = start(
    'npm',
    ['--prefix', 'frontend', 'run', 'dev', '--', '--host', '127.0.0.1', '--port', String(webPort), '--strictPort'],
    { VITE_API_PROXY_TARGET: `http://127.0.0.1:${apiPort}` },
  );
  web.stdout.on('data', (chunk) => process.stdout.write(chunk));
  web.stderr.on('data', (chunk) => process.stderr.write(chunk));
  await waitFor(web.stdout, /Local:/, 60_000);

  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext();
  await context.addCookies([
    {
      name: 'helixweave_session',
      value: 'a'.repeat(43),
      domain: '127.0.0.1',
      path: '/',
      httpOnly: true,
      sameSite: 'Lax',
    },
    {
      name: 'helixweave_csrf',
      value: 'b'.repeat(43),
      domain: '127.0.0.1',
      path: '/',
      httpOnly: false,
      sameSite: 'Lax',
    },
  ]);
  const page = await context.newPage();
  const external = [];
  const ranged = [];
  const apiResponses = [];
  const consoleMessages = [];
  page.on('console', (message) => consoleMessages.push(`${message.type()}: ${message.text()}`));
  page.on('pageerror', (error) => consoleMessages.push(`pageerror: ${error.message}`));
  page.on('response', (response) => {
    const url = response.url();
    if (url.includes('/api/')) apiResponses.push(`${response.status()} ${url}`);
    if (url.includes('/browser-tracks/') && response.status() === 206) {
      ranged.push(url);
    }
  });
  page.on('requestfailed', (request) => {
    apiResponses.push(`FAIL ${request.failure()?.errorText ?? ''} ${request.url()}`);
  });
  await page.route('**/*', (route) => {
    const host = new URL(route.request().url()).hostname;
    if (host !== '127.0.0.1' && host !== 'localhost') {
      external.push(route.request().url());
      return route.abort();
    }
    return route.continue();
  });
  await page.goto(`http://127.0.0.1:${webPort}/runs/${runId}?view=artifacts`);
  await page.getByRole('button', { name: 'View tracks in genome browser' }).click();
  const viewer = page.getByTestId('washu-genome-viewer');
  await viewer.waitFor({ timeout: 30_000 });
  try {
    await page.getByText(/bp in \d+px/).waitFor({ timeout: 20_000 });
  } catch (error) {
    await page.screenshot({ path: '/tmp/washu-browser-proof.png', fullPage: true });
    throw new Error(
      `${error instanceof Error ? error.message : 'draw timeout'}\n` +
        `${apiResponses.filter((line) => line.includes('browser-tracks') || line.startsWith('FAIL')).join('\n')}\n` +
        `${consoleMessages.filter((line) => line.startsWith('error') || line.startsWith('pageerror')).slice(-15).join('\n')}`,
    );
  }
  await page.screenshot({ path: '/tmp/washu-browser-proof.png', fullPage: true });
  const text = await page.locator('body').innerText();
  if (!text.includes('PET records') || !text.includes('Signal')) {
    throw new Error('Interaction and signal labels were not shown together.');
  }
  if (!ranged.some((url) => url.includes('artifact-pet')) || !ranged.some((url) => url.includes('artifact-signal'))) {
    throw new Error(`Expected ranged reads for both tracks, saw ${ranged.join('\n')}`);
  }
  if (external.length) {
    throw new Error(`External hosts were requested: ${external.join('\n')}`);
  }
  await page.goto('about:blank');
  await context.close();
  await browser.close();
  console.log('WASHU_PROOF ok');
} finally {
  await stopAll();
}
