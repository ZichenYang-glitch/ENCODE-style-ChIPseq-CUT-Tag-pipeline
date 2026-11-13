import { spawn } from 'node:child_process';
import { existsSync, mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import net from 'node:net';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import zlib from 'node:zlib';
import { chromium } from '@playwright/test';

const root = new URL('../..', import.meta.url).pathname;
const python =
  process.env.HELIXWEAVE_PYTHON ?? join(root, '.local/envs/ci-fast/bin/python');
const browsers =
  process.env.PLAYWRIGHT_BROWSERS_PATH ??
  '/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/ms-playwright';
if (existsSync(browsers)) process.env.PLAYWRIGHT_BROWSERS_PATH = browsers;

const negative = process.argv.includes('--negative');
const work = process.env.WASHU_PROOF_DIR
  ? (mkdirSync(process.env.WASHU_PROOF_DIR, { recursive: true }), process.env.WASHU_PROOF_DIR)
  : mkdtempSync(join(tmpdir(), 'helixweave-washu-proof-'));
const screenshotPath = join(work, negative ? 'negative.png' : 'positive.png');
const logPath = join(work, 'proof.log');
const children = [];
const reads = [];
const apiUrls = [];
const consoleMessages = [];
const external = [];
let browser;
let preflightOk = false;
let readFailed = false;
let apiPort = 0;
let webPort = 0;

function log(line) {
  writeFileSync(logPath, `${line}\n`, { flag: 'a' });
}

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
      reject(new Error(`timed out waiting for ${pattern}\n${text.slice(-2000)}`));
    }, timeoutMs);
    stream.on('data', (chunk) => {
      text += chunk.toString();
      const match = text.match(pattern);
      if (match) {
        clearTimeout(timer);
        resolve({ match, text });
      }
    });
  });
}

function reserve(start) {
  const skipped = new Set([4174, 5173, 8765]);
  return new Promise((resolve, reject) => {
    const attempt = (port) => {
      if (port > start + 30) {
        reject(new Error(`no free port from ${start}`));
        return;
      }
      if (skipped.has(port)) {
        attempt(port + 1);
        return;
      }
      const server = net.createServer();
      server.once('error', () => attempt(port + 1));
      server.listen(port, '127.0.0.1', () => {
        server.close(() => resolve(port));
      });
    };
    attempt(start);
  });
}

function portFree(port) {
  return new Promise((resolve) => {
    const server = net.createServer();
    server.once('error', () => resolve(false));
    server.listen(port, '127.0.0.1', () => {
      server.close(() => resolve(true));
    });
  });
}

async function stopAll(ports) {
  if (browser) {
    await browser.close().catch(() => {});
    browser = undefined;
  }
  for (const child of children) {
    try {
      process.kill(-child.pid, 'SIGTERM');
    } catch {
      child.kill('SIGTERM');
    }
  }
  await new Promise((resolve) => setTimeout(resolve, 400));
  for (const child of children) {
    try {
      process.kill(-child.pid, 'SIGKILL');
    } catch {
      // The process group has already exited.
    }
  }
  await new Promise((resolve) => setTimeout(resolve, 200));
  const busy = [];
  for (const port of ports) {
    if (!(await portFree(port))) busy.push(port);
  }
  if (busy.length) throw new Error(`proof ports still listening: ${busy.join(',')}`);
}

function trackOf(url) {
  if (url.includes('artifact-pet')) return 'pet';
  if (url.includes('artifact-peak')) return 'peak';
  if (url.includes('artifact-signal')) return 'signal';
  return 'other';
}

function inflateBgzf(buffer) {
  const parts = [];
  let offset = 0;
  while (offset + 18 <= buffer.length) {
    if (buffer[offset] !== 0x1f || buffer[offset + 1] !== 0x8b) break;
    const extraLength = buffer.readUInt16LE(offset + 10);
    let blockSize = null;
    let cursor = offset + 12;
    const extraEnd = cursor + extraLength;
    while (cursor + 4 <= extraEnd && cursor + 4 <= buffer.length) {
      const subfieldLength = buffer.readUInt16LE(cursor + 2);
      if (buffer[cursor] === 66 && buffer[cursor + 1] === 67 && subfieldLength === 2) {
        blockSize = buffer.readUInt16LE(cursor + 4);
      }
      cursor += 4 + subfieldLength;
    }
    if (blockSize == null) break;
    const blockLength = blockSize + 1;
    if (offset + blockLength > buffer.length) break;
    parts.push(zlib.gunzipSync(buffer.subarray(offset, offset + blockLength)));
    offset += blockLength;
  }
  return Buffer.concat(parts);
}

function signalHasValue(bodies) {
  const needle = Buffer.from('00009040', 'hex');
  for (const body of bodies) {
    if (body.includes(Buffer.from('chr1')) && body.includes(needle)) return true;
    for (let start = 0; start < body.length - 2; start += 1) {
      if (body[start] !== 0x78 || body[start + 1] !== 0x9c) continue;
      for (let end = start + 8; end <= body.length; end += 1) {
        try {
          const inflated = zlib.inflateSync(body.subarray(start, end));
          if (inflated.includes(needle)) return true;
        } catch {
          // This slice is not a complete zlib member.
        }
      }
    }
  }
  return false;
}

function assertProof(viewerText, pageText, geometry, viewHits) {
  const preflight = reads.filter((read) => read.preflight);
  const data = reads.filter((read) => !read.preflight && !read.index);
  const indexes = reads.filter((read) => !read.preflight && read.index);
  preflightOk = ['pet', 'peak', 'signal'].every((track) =>
    preflight.some((read) => read.track === track && read.status === 206),
  );
  readFailed = reads.some((read) => !read.preflight && read.status >= 400);
  const loadErrors = consoleMessages.filter(
    (line) => line.startsWith('error:') || line.startsWith('pageerror:'),
  );
  if (!preflightOk) {
    throw new Error(`missing one-byte preflight 206\n${JSON.stringify(preflight)}`);
  }
  for (const track of ['pet', 'peak', 'signal']) {
    const hit = data.find(
      (read) =>
        read.track === track &&
        (read.status === 200 || read.status === 206) &&
        read.body.length > 1 &&
        !read.contentRange.startsWith('bytes 0-0/'),
    );
    if (!hit) throw new Error(`missing ${track} data read beyond preflight`);
  }
  for (const track of ['pet', 'peak']) {
    const hit = indexes.find(
      (read) =>
        read.track === track &&
        (read.status === 200 || read.status === 206) &&
        read.body.length > 1,
    );
    if (!hit) throw new Error(`missing ${track} index read`);
  }
  const petText = Buffer.concat(
    data.filter((read) => read.track === 'pet').map((read) => inflateBgzf(read.body)),
  ).toString('utf8');
  const peakText = Buffer.concat(
    data.filter((read) => read.track === 'peak').map((read) => inflateBgzf(read.body)),
  ).toString('utf8');
  if (!petText.includes('chr1:1600-1800,1') || petText.includes('\t211\t')) {
    throw new Error(`PET bytes were not the display anchors\n${petText.slice(0, 300)}`);
  }
  if (!peakText.includes('peak-a\t80')) {
    throw new Error(`peak bytes did not include peak-a score 80\n${peakText.slice(0, 300)}`);
  }
  const signalBodies = data.filter((read) => read.track === 'signal').map((read) => read.body);
  if (!signalBodies.some((body) => body.includes(Buffer.from('chr1')))) {
    throw new Error('bigWig read did not include the chr1 header');
  }
  if (!signalHasValue(signalBodies)) {
    throw new Error('bigWig read did not include the interval value 4.5');
  }
  if (!pageText.includes('PET records') || !pageText.includes('not a significant loop')) {
    throw new Error('PET label was missing or no longer says it is not a significant loop');
  }
  if (!pageText.includes('Signal') || !/2\d{3} bp in \d+px/.test(viewerText)) {
    throw new Error(`signal label or the chr1:500-2700 window was not drawn\n${viewerText.slice(0, 400)}`);
  }
  const onWindow = viewHits.some(
    (hit) => hit.start >= 0 && hit.start <= 1000 && hit.end >= 2000 && hit.end <= 4000,
  );
  if (!onWindow) {
    throw new Error(`viewer is not on chr1:500-2700 ${JSON.stringify(viewHits)}`);
  }
  const knownWindow = viewHits.some(
    (hit) => hit.start >= 400 && hit.start <= 600 && hit.end >= 2600 && hit.end <= 2800,
  );
  if (!knownWindow || !/\b1000\b/.test(viewerText) || !/\b2000\b/.test(viewerText)) {
    throw new Error(
      `known window was not displayed\n${viewerText.slice(0, 400)}\n${JSON.stringify(viewHits)}`,
    );
  }
  const drawn = geometry.filter((box) => box.width >= 200 && box.height >= 24);
  if (drawn.length < 1) {
    throw new Error(`viewer canvases do not show drawn tracks\n${JSON.stringify(geometry)}`);
  }
  if (loadErrors.length || readFailed) {
    throw new Error(`load errors\n${loadErrors.join('\n')}\nfailed=${readFailed}`);
  }
  if (external.length) throw new Error(`external hosts\n${external.join('\n')}`);
  return {
    preflight: preflight.length,
    data: data.map((read) => `${read.track}:${read.status}:${read.body.length}`),
    indexes: indexes.map((read) => `${read.track}:${read.status}:${read.body.length}`),
    drawn,
    work,
  };
}

try {
  apiPort = await reserve(18765);
  webPort = await reserve(14174);
  log(`ports api=${apiPort} web=${webPort} negative=${negative} work=${work}`);
  const api = start(
    python,
    [
      'test/browser/washu_fixture_server.py',
      '--port',
      String(apiPort),
      '--database',
      join(work, 'platform.db'),
      '--workspace',
      join(work, 'workspaces'),
    ],
    negative ? { WASHU_PROOF_REPLACE_INDEX: '1' } : {},
  );
  api.stderr.on('data', (chunk) => process.stderr.write(chunk));
  const ready = await waitFor(api.stdout, /WASHU_FIXTURE run_id=(\S+)/, 60_000);
  const runId = ready.match[1];
  const web = start(
    'npm',
    [
      '--prefix',
      'frontend',
      'run',
      'dev',
      '--',
      '--host',
      '127.0.0.1',
      '--port',
      String(webPort),
      '--strictPort',
    ],
    { VITE_API_PROXY_TARGET: `http://127.0.0.1:${apiPort}` },
  );
  web.stdout.on('data', (chunk) => process.stdout.write(chunk));
  web.stderr.on('data', (chunk) => process.stderr.write(chunk));
  await waitFor(web.stdout, /Local:/, 60_000);

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
  const page = await context.newPage();
  page.on('console', (message) => {
    consoleMessages.push(`${message.type()}: ${message.text()}`);
  });
  page.on('pageerror', (error) => {
    consoleMessages.push(`pageerror: ${error.message}`);
  });
  page.on('request', (request) => {
    const url = request.url();
    if (url.includes('/api/')) apiUrls.push(`${request.method()} ${request.headers().range ?? '-'} ${url}`);
  });
  page.on('response', (response) => {
    const url = response.url();
    if (!url.includes('/bytes/')) return;
    const range = response.request().headers().range ?? '';
    const entry = {
      url,
      range,
      status: response.status(),
      contentRange: response.headers()['content-range'] ?? '',
      body: Buffer.alloc(0),
      track: trackOf(url),
      index: url.includes('.tbi'),
      preflight: range === 'bytes=0-0',
    };
    reads.push(entry);
    // Do not await here. Waiting for the body inside this handler deadlocks
    // the route and drops the read before it is recorded.
    void response.body().then((buf) => {
      entry.body = Buffer.from(buf);
    }).catch(() => {});
  });
  page.on('requestfailed', (request) => {
    if (!request.url().includes('/bytes/')) return;
    reads.push({
      url: request.url(),
      range: request.headers().range ?? '',
      status: 0,
      contentRange: '',
      body: Buffer.alloc(0),
      track: trackOf(request.url()),
      index: request.url().includes('.tbi'),
      preflight: (request.headers().range ?? '') === 'bytes=0-0',
    });
  });
  await page.route('**/*', async (route) => {
    const url = route.request().url();
    const host = new URL(url).hostname;
    if (host !== '127.0.0.1' && host !== 'localhost') {
      external.push(url);
      return route.abort();
    }
    return route.continue();
  });
  await page.goto(`http://127.0.0.1:${webPort}/runs/${runId}?view=artifacts`);
  await page.getByRole('button', { name: 'View tracks in genome browser' }).click();
  const viewer = page.getByTestId('washu-genome-viewer');
  await viewer.waitFor({ timeout: 30_000 });
  if (negative) {
    const started = Date.now();
    while (Date.now() - started < 15_000) {
      if (reads.some((read) => !read.preflight && read.status >= 400)) break;
      await page.waitForTimeout(200);
    }
  } else {
    await viewer.getByText(/bp in \d+px/).waitFor({ timeout: 25_000 });
    const started = Date.now();
    while (Date.now() - started < 20_000) {
      const ready = ['pet', 'peak', 'signal'].every((track) =>
        reads.some(
          (read) =>
            read.track === track &&
            !read.preflight &&
            !read.index &&
            (read.status === 200 || read.status === 206) &&
            read.body.length > 1,
        ),
      );
      const indexed = ['pet', 'peak'].every((track) =>
        reads.some(
          (read) =>
            read.track === track &&
            read.index &&
            !read.preflight &&
            (read.status === 200 || read.status === 206) &&
            read.body.length > 1,
        ),
      );
      if (ready && indexed) break;
      await page.waitForTimeout(200);
    }
  }
  if (!negative) {
    await page.waitForFunction(() => {
      const input = document.querySelector('[aria-label="Genomic region"]');
      return input instanceof HTMLInputElement && input.value === 'chr1:500-2700';
    }, undefined, { timeout: 5_000 });
    await page
      .waitForFunction(() => {
        const viewer = document.querySelector('[data-testid="washu-genome-viewer"]');
        return viewer instanceof HTMLElement && !viewer.innerText.includes('Loading view');
      }, undefined, { timeout: 8_000 })
      .catch(() => {});
    await page.waitForTimeout(2_500);
  }
  await page.screenshot({ path: screenshotPath, fullPage: true });
  const viewerText = await viewer.innerText();
  const pageText = await page.locator('body').innerText();
  log(viewerText);
  const sessionState = await page.evaluate(() => {
    const root = document.querySelector('[data-testid="washu-genome-viewer"]');
    const key = Object.getOwnPropertyNames(root ?? {}).find((name) =>
      name.startsWith('__reactFiber'),
    );
    const fiber = key ? root[key] : null;
    const hits = [];
    const stack = fiber ? [fiber] : [];
    let guard = 0;
    while (stack.length && guard < 8000 && hits.length < 6) {
      guard += 1;
      const node = stack.pop();
      const props = node?.memoizedProps;
      const region = props?.viewRegion;
      if (region && typeof region.getWidth === 'function') {
        let text = '';
        try {
          text = region.toString();
        } catch {
          text = '';
        }
        let start = null;
        let end = null;
        try {
          const coords = region.getContextCoordinates();
          start = coords.start;
          end = coords.end;
        } catch {
          start = null;
        }
        hits.push({ width: region.getWidth(), start, end });
      }
      if (node?.child) stack.push(node.child);
      if (node?.sibling) stack.push(node.sibling);
    }
    return { hits, guard };
  });
  log(`SESSION ${JSON.stringify(sessionState)}`);
  log(
    reads
      .map((read) => `${read.track} ${read.index ? 'index' : 'data'} ${read.status} ${read.range || '-'} ${read.body.length}`)
      .join('\n'),
  );
  const geometry = await viewer.locator('canvas').evaluateAll((nodes) =>
    nodes.map((node) => {
      const box = node.getBoundingClientRect();
      return { width: box.width, height: box.height };
    }),
  );
  const evidence = assertProof(viewerText, pageText, geometry, sessionState.hits);
  if (negative) {
    throw new Error(`negative case unexpectedly satisfied the proof ${JSON.stringify(evidence)}`);
  }
  await stopAll([apiPort, webPort]);
  console.log(`WASHU_PROOF ok ${JSON.stringify({ ...evidence, screenshot: screenshotPath })}`);
} catch (error) {
  const message = error instanceof Error ? error.message : String(error);
  log(apiUrls.join('\n'));
  log(message);
  log(consoleMessages.filter((line) => line.startsWith('error') || line.startsWith('pageerror')).join('\n'));
  const rejection = reads
    .filter((read) => !read.preflight && read.status === 409)
    .map((read) => read.body.toString('utf8'))
    .join('\n');
  const replaced = rejection.includes('derivative_changed') && !rejection.includes(work);
  if (negative && preflightOk && replaced) {
    console.error(`WASHU_PROOF rejected ${message.split('\n')[0]}`);
    console.error(`proof work ${work}`);
    process.exitCode = 1;
  } else {
    console.error(message);
    console.error(`proof work ${work}`);
    process.exitCode = 1;
  }
  try {
    await stopAll([apiPort, webPort].filter((port) => port > 0));
  } catch (stopError) {
    console.error(stopError instanceof Error ? stopError.message : stopError);
  }
}
