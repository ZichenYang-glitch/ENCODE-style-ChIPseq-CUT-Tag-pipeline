import { spawn } from 'node:child_process';
import { existsSync, mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import net from 'node:net';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { chromium } from '@playwright/test';

const root = new URL('../..', import.meta.url).pathname;
const python =
  process.env.HELIXWEAVE_PYTHON ?? join(root, '.local/envs/ci-fast/bin/python');
const browsers =
  process.env.PLAYWRIGHT_BROWSERS_PATH ??
  '/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/ms-playwright';
if (existsSync(browsers)) process.env.PLAYWRIGHT_BROWSERS_PATH = browsers;

const work = process.env.WASHU_SELECTION_DIR
  ? (mkdirSync(process.env.WASHU_SELECTION_DIR, { recursive: true }),
    process.env.WASHU_SELECTION_DIR)
  : mkdtempSync(join(tmpdir(), 'helixweave-washu-selection-'));
const evidencePath = join(work, 'network.json');
const children = [];
const byteRequests = [];
let browser;
let page;
let port = 0;
let armPath = '';
let runId = '';

function start(command, args) {
  const child = spawn(command, args, {
    cwd: root,
    env: process.env,
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
      reject(new Error(`timed out waiting for ${pattern}\n${text.slice(-2000)}`));
    }, timeoutMs);
    stream.on('data', (chunk) => {
      text += chunk.toString();
      const match = text.match(pattern);
      if (match) {
        clearTimeout(timer);
        resolve(match);
      }
    });
  });
}

function reserve(startPort) {
  return new Promise((resolve, reject) => {
    const attempt = (candidate) => {
      if (candidate > startPort + 30) {
        reject(new Error(`no free port from ${startPort}`));
        return;
      }
      const server = net.createServer();
      server.once('error', () => attempt(candidate + 1));
      server.listen(candidate, '127.0.0.1', () => {
        server.close(() => resolve(candidate));
      });
    };
    attempt(startPort);
  });
}

async function stopAll() {
  if (browser) await browser.close().catch(() => {});
  for (const child of children) {
    try {
      process.kill(-child.pid, 'SIGTERM');
    } catch {
      child.kill('SIGTERM');
    }
  }
  await new Promise((resolve) => setTimeout(resolve, 300));
}

function urlsFor(artifactId, since = 0) {
  return byteRequests.slice(since).filter((entry) => entry.url.includes(`/${artifactId}/`));
}

async function settle(quietMs, timeoutMs) {
  const started = Date.now();
  let lastChange = Date.now();
  let count = byteRequests.length;
  while (Date.now() - started < timeoutMs) {
    if (byteRequests.length !== count) {
      count = byteRequests.length;
      lastChange = Date.now();
    } else if (Date.now() - lastChange >= quietMs) {
      return;
    }
    await page.waitForTimeout(100);
  }
}

async function refresh(status) {
  const pending = page.waitForResponse(
    (response) =>
      response.url().includes('/browser-session') &&
      response.request().method() === 'GET' &&
      response.status() === status,
    { timeout: 20_000 },
  );
  await page.getByRole('button', { name: 'Refresh list' }).click();
  const response = await pending;
  const body = await response.json();
  return body;
}

async function main() {
  port = await reserve(4310);
  const database = join(work, 'selection.sqlite');
  const workspace = join(work, 'workspace');
  mkdirSync(workspace, { recursive: true });
  const server = start(python, [
    'test/browser/washu_selection_server.py',
    '--port',
    String(port),
    '--database',
    database,
    '--workspace',
    workspace,
  ]);
  server.stderr.on('data', (chunk) => {
    writeFileSync(join(work, 'server.log'), chunk.toString(), { flag: 'a' });
  });
  const ready = await waitFor(
    server.stdout,
    /WASHU_SELECTION run_id=(\S+) arm=(\S+)/,
    60_000,
  );
  runId = ready[1];
  armPath = ready[2];
  browser = await chromium.launch({ headless: true });
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
  page = await context.newPage();
  page.on('console', (message) => {
    writeFileSync(
      join(work, 'console.log'),
      `${message.type()}: ${message.text()}\n`,
      { flag: 'a' },
    );
  });
  page.on('pageerror', (error) => {
    writeFileSync(join(work, 'console.log'), `pageerror: ${error.message}\n`, { flag: 'a' });
  });
  page.on('request', (request) => {
    const url = request.url();
    if (url.includes('/bytes/')) {
      byteRequests.push({ url, method: request.method() });
    }
  });
  await page.goto(`http://127.0.0.1:${port}/viewer.html?run=${runId}`);
  const trackA = page.getByRole('checkbox', { name: /track-a\.narrowPeak/ });
  const trackB = page.getByRole('checkbox', { name: /track-b\.bedpe\.gz/ });
  const trackC = page.getByRole('checkbox', { name: /track-c\.narrowPeak/ });
  await trackA.waitFor({ timeout: 20_000 });
  if (byteRequests.length !== 0) {
    throw new Error(`opened with track reads\n${JSON.stringify(byteRequests)}`);
  }
  await page.screenshot({ path: join(work, 'open.png'), fullPage: true });
  await trackA.check();
  await trackB.check();
  const started = Date.now();
  while (Date.now() - started < 20_000) {
    if (urlsFor('track-a').length > 0 && urlsFor('track-b').length > 0) break;
    await page.waitForTimeout(100);
  }
  if (urlsFor('track-a').length === 0 || urlsFor('track-b').length === 0) {
    throw new Error(`selected tracks were not read\n${JSON.stringify(byteRequests)}`);
  }
  if (urlsFor('track-c').length !== 0) {
    throw new Error(`unselected track was read\n${JSON.stringify(urlsFor('track-c'))}`);
  }
  await page.screenshot({ path: join(work, 'selected.png'), fullPage: true });
  await trackA.uncheck();
  await settle(800, 15_000);
  const removedAt = byteRequests.length;
  const sameGeneration = await refresh(200);
  if (!(await trackB.isChecked()) || (await trackA.isChecked()) || (await trackC.isChecked())) {
    throw new Error('same-generation refresh did not keep only track B');
  }
  if (urlsFor('track-a', removedAt).length !== 0 || urlsFor('track-c', removedAt).length !== 0) {
    throw new Error('refresh read a track that was not selected');
  }
  await trackB.uncheck();
  await settle(800, 15_000);
  const emptyAt = byteRequests.length;
  await page.waitForTimeout(1200);
  if (byteRequests.length !== emptyAt) {
    throw new Error('tracks were read after the last selection was removed');
  }
  writeFileSync(armPath, 'change\n');
  const conflictPending = page.waitForResponse(
    (response) =>
      response.url().includes('/browser-session') && response.status() === 409,
    { timeout: 20_000 },
  );
  await page.getByRole('button', { name: 'Refresh list' }).click();
  const conflict = await conflictPending;
  const conflictBody = await conflict.json();
  const reason = conflictBody.issues?.[0]?.context?.reason_code;
  if (reason !== 'generation_changed') {
    throw new Error(`expected generation_changed, received ${JSON.stringify(conflictBody)}`);
  }
  await page.getByText('These files changed.').waitFor({ timeout: 10_000 });
  const duringConflict = byteRequests.length;
  await page.waitForTimeout(1000);
  if (byteRequests.length !== duringConflict) {
    throw new Error('new files were read while the generation conflict was showing');
  }
  await page.screenshot({ path: join(work, 'conflict.png'), fullPage: true });
  const recovered = await refresh(200);
  await trackA.waitFor({ timeout: 10_000 });
  if (
    (await trackA.isChecked()) ||
    (await trackB.isChecked()) ||
    (await trackC.isChecked())
  ) {
    throw new Error('the new generation kept the previous selection');
  }
  const beforeReselect = byteRequests.length;
  await page.waitForTimeout(1000);
  if (byteRequests.length !== beforeReselect) {
    throw new Error('the new generation was read before it was chosen again');
  }
  if (recovered.artifact_generation === sameGeneration.artifact_generation) {
    throw new Error('generation did not change across the conflict');
  }
  await trackA.check();
  const reread = Date.now();
  while (Date.now() - reread < 20_000) {
    if (urlsFor('track-a', beforeReselect).length > 0) break;
    await page.waitForTimeout(100);
  }
  const rereadUrls = urlsFor('track-a', beforeReselect).map((entry) => entry.url);
  if (
    rereadUrls.length === 0 ||
    !rereadUrls.some((url) => url.includes(`/${recovered.artifact_generation}/`))
  ) {
    throw new Error(`reselected track did not read the new generation\n${rereadUrls.join('\n')}`);
  }
  if (urlsFor('track-c').length !== 0) {
    throw new Error('track C was read during the selection proof');
  }
  const evidence = {
    runId,
    openedWithoutByteReads: true,
    selectedOnly: ['track-a', 'track-b'],
    sameGenerationKept: 'track-b',
    emptySelectionStoppedReads: true,
    conflictStatus: 409,
    conflictReason: reason,
    previousGeneration: sameGeneration.artifact_generation,
    nextGeneration: recovered.artifact_generation,
    reread: rereadUrls,
    byteRequests,
  };
  writeFileSync(evidencePath, `${JSON.stringify(evidence, null, 2)}\n`);
  console.log(`WASHU_SELECTION ok ${evidencePath}`);
}

main()
  .then(async () => {
    await stopAll();
  })
  .catch(async (error) => {
    const message = error instanceof Error ? error.stack ?? error.message : String(error);
    writeFileSync(join(work, 'error.log'), `${message}\n`);
    if (page) {
      const text = await page.locator('body').innerText().catch(() => '');
      const html = await page.content().catch(() => '');
      writeFileSync(join(work, 'page.txt'), text);
      writeFileSync(join(work, 'page.html'), html);
      await page.screenshot({ path: join(work, 'failure.png'), fullPage: true }).catch(() => {});
    }
    console.error(message);
    console.error(`selection work ${work}`);
    process.exitCode = 1;
    await stopAll();
  });
