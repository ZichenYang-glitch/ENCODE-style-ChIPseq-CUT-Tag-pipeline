/**
 * Chromium proof for local tracks overlaid on an authorized run.
 *
 *   PYTHONPATH=src node frontend/scripts/washu-local-overlay-proof.mjs
 *   PYTHONPATH=src node frontend/scripts/washu-local-overlay-proof.mjs --fixtures-only
 *
 * The packaged page must expose the selectors in
 * test/fixtures/browser/local-overlay-consensus.json. Fixture measurement
 * does not need that page.
 */
import { spawn } from 'node:child_process';
import { closeSync, existsSync, mkdirSync, mkdtempSync, openSync, readFileSync, readSync, writeFileSync } from 'node:fs';
import net from 'node:net';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import {
  bufferHasFloat,
  bytesForSlices,
  inflateBgzf,
  isExternalUrl,
  isLocalUpload,
  isWashuProxy,
  pairLines,
  parseRegion,
  regionCovers,
  sumSliceBytes,
  textHasGene,
} from './washu-local-overlay-observe.mjs';

const root = new URL('../..', import.meta.url).pathname;
const consensus = JSON.parse(
  readFileSync(join(root, 'test/fixtures/browser/local-overlay-consensus.json'), 'utf8'),
);
const selectors = consensus.selectors;
const names = consensus.filenames;
const fixturesOnly = process.argv.includes('--fixtures-only');
const browsers =
  process.env.PLAYWRIGHT_BROWSERS_PATH ??
  '/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/ms-playwright';
if (existsSync(browsers)) process.env.PLAYWRIGHT_BROWSERS_PATH = browsers;

const work = process.env.WASHU_LOCAL_OVERLAY_DIR
  ? (mkdirSync(process.env.WASHU_LOCAL_OVERLAY_DIR, { recursive: true }), process.env.WASHU_LOCAL_OVERLAY_DIR)
  : mkdtempSync(join(tmpdir(), 'helixweave-local-overlay-'));
const fixtureDir = join(work, 'fixtures');
const logPath = join(work, 'proof.log');
const children = [];
const requests = [];
const bodies = [];
let browser;
let page;

function log(line) {
  writeFileSync(logPath, `${line}\n`, { flag: 'a' });
}

function pythonExecutable() {
  const candidates = [
    process.env.HELIXWEAVE_PYTHON,
    join(root, '.local/envs/ci-fast/bin/python'),
    '/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/envs/ci-fast/bin/python',
  ].filter(Boolean);
  const found = candidates.find((candidate) => existsSync(candidate));
  if (!found) throw new Error('HelixWeave Python was not found');
  return found;
}

function runProcess(command, args, env = {}) {
  const pythonPath = [join(root, 'src'), process.env.PYTHONPATH ?? ''].filter(Boolean).join(':');
  const child = spawn(command, args, {
    cwd: root,
    env: { ...process.env, PYTHONPATH: pythonPath, PYTHONUNBUFFERED: '1', ...env },
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
  return new Promise((resolve, reject) => {
    const attempt = (port) => {
      if (port > start + 30) {
        reject(new Error(`no free port from ${start}`));
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

function readAt(file, start, length) {
  const handle = openSync(file, 'r');
  try {
    const buffer = Buffer.alloc(length);
    const got = readSync(handle, buffer, 0, length, start);
    return buffer.subarray(0, got);
  } finally {
    closeSync(handle);
  }
}

function localMarkers() {
  const first = pairLines('chr7', consensus.longrange_pairs.egfr[0]);
  return [first[0], 'egfr-peak', 'tp53-peak'];
}

function assertCleanNetwork() {
  const markers = localMarkers();
  const external = [];
  const uploads = [];
  const genes = [];
  const proxies = [];
  for (const request of requests) {
    const blob = `${request.method} ${request.url}\n${request.body}`;
    if (isExternalUrl(request.url)) external.push(request.url);
    if (isWashuProxy(request.url) || isWashuProxy(request.body)) proxies.push(request.url);
    if (isLocalUpload(request.method, request.body, markers)) uploads.push(request.url);
    if (textHasGene(request.url) || textHasGene(request.body)) genes.push(request.url);
  }
  if (external.length || uploads.length || genes.length || proxies.length) {
    throw new Error(
      `network left the browser ${JSON.stringify({ external, uploads, genes, proxies })}`,
    );
  }
}

async function collectReads() {
  const reads = [];
  for (const frame of page.frames()) {
    const chunk = await frame.evaluate(() => window.__localReads ?? []).catch(() => []);
    reads.push(...chunk);
  }
  return reads;
}

async function clearReads() {
  for (const frame of page.frames()) {
    await frame.evaluate(() => {
      if (Array.isArray(window.__localReads)) window.__localReads.length = 0;
    }).catch(() => {});
  }
}

function longrangeText(reads) {
  const slices = reads.filter((record) => record.name === names.longrange);
  return inflateBgzf(bytesForSlices(join(fixtureDir, names.longrange), slices, readAt)).toString('utf8');
}

async function frameWithLabel(label) {
  const deadline = Date.now() + 20_000;
  while (Date.now() < deadline) {
    for (const frame of page.frames()) {
      const count = await frame.getByLabel(label).count().catch(() => 0);
      if (count > 0) return frame;
    }
    await page.waitForTimeout(100);
  }
  throw new Error(`missing control ${label}`);
}

async function frameWithText(text) {
  const deadline = Date.now() + 20_000;
  while (Date.now() < deadline) {
    for (const frame of page.frames()) {
      const count = await frame.getByText(text, { exact: false }).count().catch(() => 0);
      if (count > 0) return frame;
    }
    await page.waitForTimeout(100);
  }
  throw new Error(`missing text ${text}`);
}

async function regionValue(scope) {
  return scope.getByLabel(selectors.region_input).inputValue();
}

async function submitRegion(scope, value) {
  const field = scope.getByLabel(selectors.region_input);
  await field.fill(value);
  await field.press('Enter');
}

async function waitForRegion(scope, predicate) {
  const deadline = Date.now() + 15_000;
  let last = '';
  while (Date.now() < deadline) {
    last = await regionValue(scope);
    const parsed = parseRegion(last);
    if (predicate(parsed, last)) return parsed;
    await page.waitForTimeout(100);
  }
  throw new Error(`region stayed ${last}`);
}

async function waitForLongrange(includes, excludes) {
  const deadline = Date.now() + 30_000;
  let last = '';
  while (Date.now() < deadline) {
    const reads = await collectReads();
    const text = longrangeText(reads);
    const dataBytes = sumSliceBytes(reads, names.longrange);
    const ready = includes.every((line) => text.includes(line))
      && excludes.every((line) => !text.includes(line))
      && dataBytes > 0;
    if (ready) return { reads, text, dataBytes };
    last = text.slice(0, 240);
    await page.waitForTimeout(200);
  }
  throw new Error(`longrange window was not read\n${last}`);
}

function coversPair(region, chrom, pair) {
  return regionCovers(region, chrom, pair.a[0], pair.a[1])
    && regionCovers(region, chrom, pair.b[0], pair.b[1]);
}

async function openLocal(port) {
  await page.goto(`http://127.0.0.1:${port}/`);
  const opener = page
    .getByRole('button', { name: selectors.open_local_tracks })
    .or(page.getByRole('link', { name: selectors.open_local_tracks }));
  await opener.click();
  await frameWithText(selectors.user_assembly_text);
  return frameWithLabel(selectors.region_input);
}

async function proveStandalone(scope, observation) {
  const tight = `${consensus.windows.egfr.chrom}:${consensus.windows.egfr.start}-${consensus.windows.egfr.end}`;
  await submitRegion(scope, tight);
  await waitForRegion(scope, (_parsed, raw) => raw.trim() === tight);
  const input = scope.getByLabel(selectors.local_file_input);
  await input.setInputFiles([
    join(fixtureDir, names.longrange),
    join(fixtureDir, names.longrange_index),
    join(fixtureDir, names.signal),
    join(fixtureDir, names.peaks),
  ]);
  const viewer = scope.getByTestId(selectors.viewer);
  await viewer.waitFor({ timeout: 30_000 });
  const first = pairLines('chr7', consensus.longrange_pairs.egfr[0]);
  const mid = pairLines('chr7', consensus.longrange_pairs.egfr[1])[0];
  const loaded = await waitForLongrange(first, [mid]);
  const dataLimit = observation.data_file_bytes / 1000;
  if (loaded.dataBytes >= dataLimit) {
    throw new Error(`coordinate view read ${loaded.dataBytes} of ${observation.data_file_bytes} data bytes`);
  }
  const signalBytes = bytesForSlices(
    join(fixtureDir, names.signal),
    loaded.reads.filter((record) => record.name === names.signal),
    readAt,
  );
  if (!bufferHasFloat(signalBytes, 4.5)) throw new Error('bigWig window did not include signal 4.5');
  const peakBytes = bytesForSlices(
    join(fixtureDir, names.peaks),
    loaded.reads.filter((record) => record.name === names.peaks),
    readAt,
  ).toString('utf8');
  if (!peakBytes.includes('chr7\t55019050\t55019250\tegfr-peak\t80')) {
    throw new Error('BED window did not include egfr-peak 55019050-55019250');
  }
  const canvases = await viewer.locator('canvas').evaluateAll((nodes) =>
    nodes.map((node) => {
      const box = node.getBoundingClientRect();
      return { width: box.width, height: box.height };
    }),
  );
  if (!canvases.some((box) => box.width >= 200 && box.height >= 24)) {
    throw new Error(`GenomeHub did not draw the window ${JSON.stringify(canvases)}`);
  }
  await page.screenshot({ path: join(work, 'coordinate.png'), fullPage: true });

  await clearReads();
  await submitRegion(scope, 'EGFR');
  const egfr = await waitForRegion(scope, (parsed) => {
    const endPair = consensus.longrange_pairs.egfr[2];
    return coversPair(parsed, 'chr7', consensus.longrange_pairs.egfr[0])
      && coversPair(parsed, 'chr7', endPair);
  });
  const egfrLines = [
    ...pairLines('chr7', consensus.longrange_pairs.egfr[0]),
    ...pairLines('chr7', consensus.longrange_pairs.egfr[2]),
  ];
  const egfrRead = await waitForLongrange(egfrLines, []);
  if (egfrRead.dataBytes >= dataLimit) {
    throw new Error(`EGFR view read ${egfrRead.dataBytes} data bytes`);
  }

  await clearReads();
  await submitRegion(scope, 'TP53');
  await waitForRegion(scope, (parsed) =>
    coversPair(parsed, 'chr17', consensus.longrange_pairs.tp53[0])
    && coversPair(parsed, 'chr17', consensus.longrange_pairs.tp53[2]),
  );
  const tp53Lines = [
    ...pairLines('chr17', consensus.longrange_pairs.tp53[0]),
    ...pairLines('chr17', consensus.longrange_pairs.tp53[2]),
  ];
  const tp53Read = await waitForLongrange(tp53Lines, []);
  const tp53Signal = bytesForSlices(
    join(fixtureDir, names.signal),
    tp53Read.reads.filter((record) => record.name === names.signal),
    readAt,
  );
  if (!bufferHasFloat(tp53Signal, 2.25)) throw new Error('TP53 window did not include signal 2.25');

  const beforeAmbiguous = await regionValue(scope);
  await submitRegion(scope, consensus.genes.ambiguous.symbol);
  const matches = scope.getByRole('listbox', { name: selectors.gene_matches });
  await matches.waitFor({ timeout: 10_000 });
  await matches.getByRole('option', { name: /chr1/ }).waitFor();
  await matches.getByRole('option', { name: /chr15/ }).waitFor();
  if ((await regionValue(scope)) !== beforeAmbiguous) {
    throw new Error('ambiguous gene navigated before a locus was chosen');
  }
  await matches.getByRole('option', { name: /chr1/ }).click();
  await waitForRegion(scope, (parsed) =>
    Boolean(parsed)
    && parsed.chrom === 'chr1'
    && parsed.start < 24894
    && parsed.end > 11120,
  );

  const beforeMissing = await regionValue(scope);
  await submitRegion(scope, consensus.genes.missing_symbol);
  await scope.getByText(selectors.no_gene_status).waitFor({ timeout: 10_000 });
  if ((await regionValue(scope)) !== beforeMissing) {
    throw new Error('missing gene changed the current region');
  }

  await scope.getByRole('button', { name: `${selectors.remove_prefix}${names.peaks}` }).click();
  await scope.getByRole('button', { name: `${selectors.remove_prefix}${names.signal}` }).waitFor();
  if (await scope.getByRole('button', { name: `${selectors.remove_prefix}${names.peaks}` }).count()) {
    throw new Error('removed peak track is still listed');
  }

  await scope.getByLabel(selectors.local_file_input).setInputFiles([
    join(fixtureDir, names.longrange_missing_index),
  ]);
  await scope.getByText(selectors.missing_index_status).waitFor({ timeout: 10_000 });
  await scope.getByRole('button', { name: `${selectors.remove_prefix}${names.signal}` }).waitFor();

  await scope.getByLabel(selectors.local_file_input).setInputFiles([
    join(fixtureDir, names.longrange_wrong),
    join(fixtureDir, names.longrange_wrong_index),
  ]);
  await scope.getByText(selectors.wrong_index_status).waitFor({ timeout: 10_000 });

  await scope.getByLabel(selectors.assembly).selectOption(selectors.assembly_mm10);
  await scope.getByText(selectors.assembly_mismatch_status).waitFor({ timeout: 10_000 });
  await scope.getByText(selectors.user_assembly_text).waitFor();
  await scope.getByLabel(selectors.assembly).selectOption(selectors.assembly_hg38);
  return { egfr, coordinateDataBytes: loaded.dataBytes, geneDataBytes: egfrRead.dataBytes };
}

async function proveRuns(port, ids) {
  await page.goto(`http://127.0.0.1:${port}/runs/${ids.run_a}?view=artifacts`);
  const openBrowser = page.getByRole('button', { name: selectors.open_run_browser });
  if (await openBrowser.count()) await openBrowser.click();
  const scope = await frameWithLabel(selectors.local_file_input);
  const signal = scope.getByRole('checkbox', { name: new RegExp(ids.run_a) }).filter({ hasText: /Signal/ });
  const signalBox = (await signal.count())
    ? signal.first()
    : scope.getByRole('checkbox', { name: /Signal/ }).first();
  await signalBox.check();
  const signalDeadline = Date.now() + 25_000;
  while (Date.now() < signalDeadline) {
    if (bodies.some((entry) => entry.url.includes(ids.signal_a) && bufferHasFloat(entry.body, 7.25))) break;
    await page.waitForTimeout(200);
  }
  if (!bodies.some((entry) => entry.url.includes(ids.signal_a) && bufferHasFloat(entry.body, 7.25))) {
    throw new Error('run A bigWig 7.25 was not read from the authorized API');
  }
  await scope.getByLabel(selectors.local_file_input).setInputFiles([
    join(fixtureDir, names.longrange),
    join(fixtureDir, names.longrange_index),
    join(fixtureDir, names.signal),
    join(fixtureDir, names.peaks),
  ]);
  await waitForLongrange(pairLines('chr7', consensus.longrange_pairs.egfr[0]), []);
  await scope.getByRole('button', { name: selectors.add_run }).click();
  const runs = scope.getByRole('listbox', { name: selectors.authorized_runs });
  await runs.getByRole('option', { name: new RegExp(ids.run_b) }).click();
  const peak = scope.getByRole('checkbox', { name: new RegExp(ids.run_b) });
  await peak.first().check();
  const peakDeadline = Date.now() + 25_000;
  while (Date.now() < peakDeadline) {
    if (bodies.some((entry) => entry.body.includes('run-b-peak'))) break;
    await page.waitForTimeout(200);
  }
  if (!bodies.some((entry) => entry.body.includes('run-b-peak'))) {
    throw new Error('run B peak was not read');
  }
  const revised = await page.request.post(`http://127.0.0.1:${port}/proof-control/revise-run-b`);
  if (!revised.ok()) throw new Error(`revision control returned ${revised.status()}`);
  await scope.getByRole('button', { name: selectors.refresh_runs }).click();
  await scope.getByText(selectors.stale_status).waitFor({ timeout: 15_000 });
  await page.waitForTimeout(1000);
  if (bodies.some((entry) => entry.body.includes('run-b-revised'))) {
    throw new Error('the revised run B peak was read without being chosen again');
  }
  await scope.getByRole('button', { name: `${selectors.remove_prefix}${names.signal}` }).waitFor();
  await scope.getByRole('button', { name: selectors.add_run }).click();
  await runs.getByRole('option', { name: new RegExp(ids.run_mm10) }).click();
  await scope.getByText(selectors.assembly_mismatch_status).waitFor({ timeout: 10_000 });
  if (requests.some((request) => request.url.includes(ids.signal_mm10))) {
    throw new Error('mm10 run bytes were requested on the hg38 view');
  }
}

let apiPort;
try {
  const python = pythonExecutable();
  const fixture = runProcess(python, [
    'test/browser/local_overlay_fixtures.py',
    '--dest',
    fixtureDir,
  ]);
  fixture.stderr.on('data', (chunk) => process.stderr.write(chunk));
  const fixtureReady = await waitFor(fixture.stdout, /LOCAL_OVERLAY_FIXTURE (\{.*\})/, 180_000);
  const summary = JSON.parse(fixtureReady.match[1]);
  log(`FIXTURE ${fixtureReady.match[1]}`);
  const observation = JSON.parse(readFileSync(join(fixtureDir, 'observation.json'), 'utf8'));
  if (fixturesOnly) {
    console.log(`LOCAL_OVERLAY_PROOF fixtures ${JSON.stringify(summary)}`);
    await stopAll([]);
    process.exit(0);
  }

  apiPort = await reserve(19110);
  log(`port ${apiPort} work ${work}`);
  const api = runProcess(python, [
    'test/browser/washu_local_overlay_server.py',
    '--port',
    String(apiPort),
    '--database',
    join(work, 'platform.db'),
    '--workspace',
    join(work, 'workspaces'),
  ]);
  api.stderr.on('data', (chunk) => process.stderr.write(chunk));
  const ready = await waitFor(
    api.stdout,
    /LOCAL_OVERLAY_SERVER run_a=(\S+) run_b=(\S+) run_mm10=(\S+) signal_a=(\S+) peak_b=(\S+) signal_mm10=(\S+)/,
    120_000,
  );
  const ids = {
    run_a: ready.match[1],
    run_b: ready.match[2],
    run_mm10: ready.match[3],
    signal_a: ready.match[4],
    peak_b: ready.match[5],
    signal_mm10: ready.match[6],
  };
  const { chromium } = await import('@playwright/test');
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
  await context.addInitScript(() => {
    window.__localReads = [];
    window.__helixweaveCsp = [];
    document.addEventListener('securitypolicyviolation', (event) => {
      window.__helixweaveCsp.push(event.blockedURI);
    });
    const originalSlice = Blob.prototype.slice;
    Blob.prototype.slice = function slice(start, end, contentType) {
      const size = this.size || 0;
      const from = start == null ? 0 : start < 0 ? Math.max(0, size + start) : start;
      const to = end == null ? size : end < 0 ? Math.max(0, size + end) : end;
      window.__localReads.push({
        name: typeof this.name === 'string' ? this.name : '',
        start: from,
        end: to,
        bytes: Math.max(0, to - from),
      });
      return originalSlice.call(this, start, end, contentType);
    };
    const originalBuffer = Blob.prototype.arrayBuffer;
    Blob.prototype.arrayBuffer = function arrayBuffer() {
      window.__localReads.push({
        name: typeof this.name === 'string' ? this.name : '',
        start: 0,
        end: this.size || 0,
        bytes: this.size || 0,
      });
      return originalBuffer.call(this);
    };
  });
  page = await context.newPage();
  page.on('request', (request) => {
    requests.push({
      method: request.method(),
      url: request.url(),
      body: request.postData() ?? '',
    });
  });
  page.on('response', (response) => {
    const entry = { url: response.url(), status: response.status(), body: Buffer.alloc(0) };
    bodies.push(entry);
    void response.body().then((buffer) => {
      entry.body = Buffer.from(buffer);
    }).catch(() => {});
  });
  const standalone = await proveStandalone(await openLocal(apiPort), observation);
  await proveRuns(apiPort, ids);
  assertCleanNetwork();
  const violations = [];
  for (const frame of page.frames()) {
    const chunk = await frame.evaluate(() => window.__helixweaveCsp ?? []).catch(() => []);
    violations.push(...chunk);
  }
  if (violations.some((uri) => isWashuProxy(uri) || (uri && isExternalUrl(uri)))) {
    throw new Error(`proxy or external CSP violation ${JSON.stringify(violations)}`);
  }
  const result = {
    work,
    fixture: summary,
    coordinateDataBytes: standalone.coordinateDataBytes,
    geneDataBytes: standalone.geneDataBytes,
    egfr: standalone.egfr,
    requests: requests.length,
  };
  writeFileSync(join(work, 'browser.json'), `${JSON.stringify(result, null, 2)}\n`);
  await stopAll([apiPort]);
  console.log(`LOCAL_OVERLAY_PROOF ok ${JSON.stringify(result)}`);
} catch (error) {
  const message = error instanceof Error ? error.message : String(error);
  log(message);
  if (page) await page.screenshot({ path: join(work, 'failure.png'), fullPage: true }).catch(() => {});
  console.error(message);
  console.error(`proof work ${work}`);
  try {
    await stopAll(apiPort == null ? [] : [apiPort]);
  } catch (stopError) {
    console.error(stopError instanceof Error ? stopError.message : stopError);
  }
  process.exitCode = 1;
}
