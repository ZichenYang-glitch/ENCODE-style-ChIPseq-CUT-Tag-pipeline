// Run against a packaged HelixWeave frontend using HELIXWEAVE_PYTHON.
import assert from 'node:assert/strict';
import { spawn } from 'node:child_process';
import { once } from 'node:events';
import { existsSync } from 'node:fs';
import { mkdir, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { createInterface } from 'node:readline';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { chromium, expect } from '@playwright/test';

const root = fileURLToPath(new URL('../..', import.meta.url));
const work = process.env.WASHU_OVERLAY_DIR
  ? resolve(process.env.WASHU_OVERLAY_DIR)
  : await mkdtemp(join(tmpdir(), 'helixweave-overlay-'));
await mkdir(work, { recursive: true });
const fixture = await mkdtemp(join(work, 'fixture-'));
const workspace = join(fixture, 'workspace');
// The control waits for the two surviving tracks to paint before editing the
// coordinates. The default exercises the same enabled input immediately after removal.
const waitForRemovalDraw = process.argv.includes('--wait-for-removal-draw');
const evidence = {
  source: 'synthetic',
  navigationAfterRemoval: waitForRemovalDraw ? 'after redraw' : 'immediate',
};
const expectedPeaks = {
  'peak-primary': [5500, 5700, 61],
  'peak-local': [5800, 6000, null],
  'peak-second': [6000, 6200, 73],
};
const remaining = ['peak-primary', 'peak-local'];
const allPeaks = [...remaining, 'peak-second'];
const timeout = 20_000;
const checkoutPython = join(root, '.local/envs/ci-fast/bin/python');
const python = process.env.HELIXWEAVE_PYTHON ?? (existsSync(checkoutPython) ? checkoutPython : 'python3');
const server = spawn(python, [
  '-I', 'test/browser/washu_overlay_server.py',
  '--workspace', workspace, '--database', join(fixture, 'overlay.sqlite'),
], { cwd: root, stdio: ['ignore', 'pipe', 'inherit'] });
const serverExit = once(server, 'exit');
const lines = createInterface({ input: server.stdout });
let browser;
let page;
let frame;
let network;

function decodedBody(bytes) {
  return bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes;
}

function monitorBrowser(context, origin, localFile) {
  const requests = [];
  const reads = [];
  const forbidden = [];
  const policyReports = [];
  const pageErrors = [];
  context.on('page', tab => {
    tab.on('pageerror', error => pageErrors.push(error.message));
    tab.on('console', message => {
      if (message.type() === 'error' && /minified React error|maximum update depth exceeded|Genome View Error/i.test(message.text())) {
        pageErrors.push(message.text());
      }
    });
  });
  context.on('request', request => {
    const url = new URL(request.url());
    if (url.protocol === 'data:' || url.protocol === 'blob:') return;
    requests.push(url);
    const body = request.postDataBuffer();
    const payload = [decodeURIComponent(url.href), body?.toString('utf8'), ...Object.values(request.headers())].join('\n');
    const localContent = payload.includes(localFile) || payload.includes('peak-local') || payload.includes('5800\t6000');
    const report = url.origin === origin && url.pathname === '/csp-report' && request.method() === 'POST';
    if (report) policyReports.push(JSON.parse(body.toString('utf8'))['csp-report']?.['blocked-uri']);
    if (url.origin !== origin || localContent || url.pathname.startsWith('/__washu_overlay__/') ||
        (!report && (body !== null || !['GET', 'HEAD'].includes(request.method())))) {
      forbidden.push({ method: request.method(), path: url.origin + url.pathname, localContent });
    }
  });
  context.on('response', response => {
    const url = new URL(response.url());
    if (!url.pathname.includes('/bytes/')) return;
    const entry = {
      path: url.pathname, origin: url.origin, status: response.status(),
      range: response.request().headers().range ?? null, body: null,
    };
    reads.push(entry);
    // A remount can abort an in-flight read. Only a completed response can prove a read.
    void response.body().then(bytes => { entry.body = bytes; }).catch(() => {});
  });
  return { origin, requests, reads, forbidden, policyReports, pageErrors };
}

function networkReceipt(log) {
  return {
    httpRequests: log.requests.length, byteRequests: log.requests.filter(url => url.pathname.includes('/bytes/')).length,
    externalRequests: log.requests.filter(url => url.origin !== log.origin).length,
    localUploads: log.forbidden.filter(item => item.localContent).length,
    forbidden: log.forbidden, cspReports: log.policyReports, pageErrors: log.pageErrors,
  };
}

async function platformRead(network, origin, track, peak, since = 0) {
  async function completed(path) {
    let hit;
    await expect.poll(() => {
      hit = network.reads.slice(since).find(read => read.origin === origin && read.path === path &&
        [200, 206].includes(read.status) && read.body?.length > 2);
      return Boolean(hit);
    }, { timeout, message: 'Expected a completed browser read of ' + path }).toBe(true);
    return hit;
  }
  const data = await completed(track.content_path);
  const index = await completed(track.index_path);
  const [start, end, score] = expectedPeaks[peak];
  const row = ['chr1', start, end, peak, score, '+'].join('\t');
  assert(decodedBody(data.body).toString('utf8').split('\n').includes(row));
  assert.equal(decodedBody(index.body).subarray(0, 4).toString('latin1'), 'TBI\x01');
  const receipt = ({ path, status, range, body }) => ({ path, status, range, bytes: body.length });
  return { origin, data: receipt(data), index: receipt(index), row };
}

// GenomeHub 62 exposes its visible genomic interval and BED placements in its
// committed React tree. Pair each placed feature with the rectangle actually in the SVG.
function renderedHub(host) {
  function* descendants(node) {
    if (!node) return;
    yield node;
    for (let child = node.child; child; child = child.sibling) yield* descendants(child);
  }
  const key = Object.keys(host).find(name => name.startsWith('__reactFiber'));
  if (!key) return null;
  let rootFiber = host[key];
  while (rootFiber.return) rootFiber = rootFiber.return;
  let current;
  for (const node of descendants(rootFiber.stateNode.current)) {
    if (node.stateNode === host) { current = node; break; }
  }
  const rects = [...host.querySelectorAll('svg rect')];
  let models;
  let window;
  let ruler;
  const peaks = {};
  for (const node of descendants(current)) {
    const props = node.memoizedProps;
    if (!props) continue;
    if (props.genomeName && props.storeConfig && Array.isArray(props.tracks)) {
      models = props.tracks.map(track => ({
        name: track.name, type: track.type, url: track.url ?? '', file: track.fileObj?.name ?? null,
      }));
    }
    const intervals = props.userViewRegion?.getGenomeIntervals?.();
    if (intervals?.length === 1) {
      const { chr, start, end } = intervals[0];
      window = { chr, start, end };
    }
    if (props.trackModel?.type === 'ruler' && props.windowWidth > 0) {
      const selected = props.selectedRegion?.getGenomeIntervals?.();
      if (selected?.length === 1) {
        const { chr, start, end } = selected[0];
        ruler = { chr, start, end, pixels: props.windowWidth };
      }
    }
    if (!Array.isArray(props.placements)) continue;
    for (const { feature, xSpan, placedFeatures } of props.placements) {
      if (!feature || !placedFeatures?.length || !(xSpan?.end > xSpan.start)) continue;
      const mark = rects.find(rect => Math.abs(rect.x.baseVal.value - xSpan.start) < 0.01 &&
        Math.abs(rect.width.baseVal.value - (xSpan.end - xSpan.start)) < 0.01);
      if (!mark) continue;
      const box = mark.getBoundingClientRect();
      if (box.width <= 0 || box.height <= 0 || !mark.checkVisibility()) continue;
      peaks[feature[3]] = {
        chr: feature.chr, start: feature.start, end: feature.end,
        score: feature[4] == null ? null : Number(feature[4]),
        svgWidth: box.width, svgHeight: box.height,
      };
    }
  }
  return { window, ruler, models, peaks };
}

async function expectDrawing(region, names, urls, localFile) {
  const viewer = frame.getByTestId('washu-genome-viewer');
  await expect(viewer).toHaveCount(1, { timeout });
  await viewer.scrollIntoViewIfNeeded();
  const [first, last] = region.split(':')[1].split('-').map(Number);
  let state;
  try {
    await expect.poll(async () => {
      state = await viewer.evaluate(renderedHub);
      if (state?.window?.chr !== 'chr1' || state.window.start !== first - 1 ||
          state.window.end !== last || state.models?.length !== names.length + 1) return false;
      if (state.ruler?.chr !== 'chr1' || state.ruler.start !== first - 1 || state.ruler.end !== last) return false;
      const local = state.models.some(track => track.type === 'bed' && track.file === localFile && track.url === '');
      if (!local || !state.models.some(track => track.type === 'ruler') ||
          !urls.every(url => state.models.some(track => track.type === 'bed' && track.url === url))) return false;
      if (Object.keys(state.peaks).length !== names.length) return false;
      return names.every(name => {
        const peak = state.peaks[name];
        const [start, end, score] = expectedPeaks[name];
        const width = (end - start) * state.ruler.pixels / (state.ruler.end - state.ruler.start);
        return peak?.chr === 'chr1' && peak.start === start && peak.end === end &&
          peak.score === score && Math.abs(peak.svgWidth - width) < 1;
      });
    }, { timeout }).toBe(true);
  } catch (error) {
    evidence.drawingFailure = { requested: region, expectedPeaks: names, observed: state };
    await viewer.screenshot({ path: join(work, 'overlay-window-mismatch.png'), timeout: 2000 }).catch(() => {});
    throw new Error('Expected ' + names.join(', ') + ' at ' + region + '. Observed ' + JSON.stringify(state), { cause: error });
  }
  return { viewerCount: 1, ...state };
}

async function main(manifest) {
  const origin = 'http://127.0.0.1:' + manifest.port;
  browser = await chromium.launch();
  const context = await browser.newContext({ viewport: { width: 1500, height: 1000 } });
  const login = await context.request.post(origin + '/api/v1/auth/login', {
    headers: { Origin: origin }, data: { username: 'overlay-reviewer', password: manifest.password },
  });
  assert.equal(login.status(), 200);
  evidence.build = {
    chromium: browser.version(), package: manifest.package_version,
    wheelSha256: manifest.wheel_sha256, frontendIdentity: manifest.frontend_identity,
  };
  const api = context.request;
  const primaryId = manifest.primary_run_id;
  const secondId = manifest.second_run_id;
  const foreignId = manifest.foreign_run_id;
  const ids = [primaryId, secondId, foreignId];
  const history = await api.get(origin + '/api/v1/runs?limit=100');
  assert.equal(history.status(), 200);
  assert(!(await history.text()).includes(workspace));
  const listed = (await history.json()).runs.map(run => run.run_id);
  const tracks = [];
  for (const [index, id] of ids.entries()) {
    assert(listed.includes(id));
    const session = await api.get(origin + '/api/v1/runs/' + id + '/browser-session');
    assert.equal(session.status(), 200);
    assert(!(await session.text()).includes(workspace));
    const data = await session.json();
    const scaffold = index === 2 ? 'mm10' : 'hg38';
    assert.equal(data.reference.viewer_scaffold, scaffold);
    assert.equal(data.reference.assembly_label, scaffold);
    const track = data.tracks.find(item => item.artifact_id === manifest.artifact_id);
    assert.equal(track?.preparation, 'ready');
    const base = '/api/v1/runs/' + id + '/browser-tracks/' + manifest.artifact_id + '/bytes/';
    assert.equal(track.content_path, base + data.artifact_generation + '/' + track.revision);
    assert.equal(track.index_path, track.content_path + '.tbi');
    tracks.push({
      ...track, run: id, generation: data.artifact_generation,
      referenceId: data.reference.revision_id, identity: data.reference.identity_sha256,
    });
  }
  const [primary, second] = tracks;
  assert(primary.identity && second.identity && primary.identity !== second.identity);
  assert.equal(primary.generation, manifest.primary_generation);
  assert.equal(primary.revision, manifest.primary_revision);
  assert.equal(second.generation, manifest.second_generation);
  assert.equal(second.revision, manifest.second_revision);
  const anonymousRanges = [];
  for (const track of [primary, second]) {
    for (const path of [track.content_path, track.index_path]) {
      const denial = await fetch(origin + path, { headers: { Range: 'bytes=0-3' } });
      assert.equal(denial.status, 401);
      const text = await denial.text();
      assert(!text.includes(workspace) && !/peak-primary|peak-second/.test(text));
      anonymousRanges.push({ path, status: denial.status });
    }
  }
  async function rejected(path, status, reason) {
    const res = await api.get(origin + path, { headers: { Range: 'bytes=0-3' } });
    assert.equal(res.status(), status);
    const text = await res.text();
    assert(text.includes(reason) && !text.includes(workspace));
    return { status, reason };
  }
  evidence.protection = {
    login: login.status(), anonymousRanges,
    missing: await rejected(primary.content_path.replace('/' + manifest.artifact_id + '/', '/missing-artifact/'), 404, 'RUN_ARTIFACT_NOT_FOUND'),
    staleRevision: await rejected(primary.content_path.replace(primary.revision, 'artifactrev-' + 'd'.repeat(64)), 409, 'revision_changed'),
    staleGeneration: await rejected(primary.content_path.replace(primary.generation, 'artifactgen-' + 'e'.repeat(64)), 409, 'generation_changed'),
  };

  network = monitorBrowser(context, origin, manifest.local_bed);
  // Record and stop any attempted off-origin transfer before contacting another server.
  await context.route('**/*', route => {
    const url = new URL(route.request().url());
    return url.origin === origin || ['blob:', 'data:'].includes(url.protocol) ? route.continue() : route.abort();
  });
  page = await context.newPage();
  const nav = await page.goto(origin + '/runs/' + primaryId + '?view=artifacts');
  const frameResponse = page.waitForResponse(r => new URL(r.url()).pathname === '/viewer.html' && r.request().isNavigationRequest());
  await page.getByRole('button', { name: 'View tracks in genome browser' }).click();
  const viewerDoc = await frameResponse;
  frame = viewerDoc.frame();
  const documentPolicy = nav.headers()['content-security-policy'];
  const viewerPolicy = viewerDoc.headers()['content-security-policy'];
  const frameOptions = viewerDoc.headers()['x-frame-options'];
  for (const policy of [documentPolicy, viewerPolicy]) {
    assert(policy.startsWith("connect-src 'self'; report-uri /csp-report"));
  }
  assert(viewerPolicy.includes("frame-ancestors 'self'"));
  assert.equal(frameOptions, 'SAMEORIGIN');
  Object.assign(evidence.protection, { documentPolicy, viewerPolicy, frameOptions });

  const article = id => frame.getByTestId('browser-run').and(frame.getByRole('article', { name: 'Run ' + id, exact: true }));
  const primaryBox = article(primaryId).getByRole('checkbox');
  const region = frame.getByRole('textbox', { name: 'Genomic region' });
  const picker = frame.getByRole('combobox', { name: 'Authorized runs' });
  async function addRun(id) {
    await picker.selectOption(id);
    await frame.getByRole('button', { name: 'Add run', exact: true }).click();
  }
  await expect(primaryBox).toBeEnabled({ timeout });
  // Include the initial delayed region settling in the zero-read observation.
  await page.waitForTimeout(1250);
  await expect(primaryBox).not.toBeChecked();
  await expect(frame.getByTestId('browser-run')).toHaveCount(1);
  await expect(frame.getByTestId('washu-genome-viewer')).toHaveCount(0);
  evidence.initial = {
    checked: await frame.getByRole('checkbox', { checked: true }).count(),
    byteRequests: network.requests.filter(url => url.pathname.includes('/bytes/')).length,
    summaryRequests: network.requests.filter(url => url.pathname.includes('/region-summary/')).length,
  };
  assert.deepEqual(evidence.initial, { checked: 0, byteRequests: 0, summaryRequests: 0 });

  await primaryBox.check();
  await frame.getByTestId('local-track-file-input').setInputFiles(join(root, 'test/fixtures/browser', manifest.local_bed));
  await expect(frame.getByTestId('local-browser-track')).toBeVisible();
  await addRun(secondId);
  const secondArticle = article(secondId);
  const secondBox = secondArticle.getByRole('checkbox');
  for (const track of [primary, second]) {
    const run = article(track.run);
    await expect(run).toHaveAttribute('data-reference-identity', track.identity);
    await expect(run).toHaveAttribute('data-generation', track.generation);
    await expect(run.getByTestId('browser-track')).toHaveAttribute('data-revision', track.revision);
  }
  await expect(secondArticle.getByTestId('reference-conflict')).toHaveCount(0);
  await secondBox.check();
  await expect(region).toHaveValue(manifest.window);
  const primaryUrl = origin + primary.content_path;
  const bothUrls = [primaryUrl, origin + second.content_path];
  const parsed = await expectDrawing(manifest.window, allPeaks, bothUrls, manifest.local_bed);
  evidence.overlay = { ...parsed, platform: [] };
  for (const [track, peak] of [[primary, 'peak-primary'], [second, 'peak-second']]) {
    evidence.overlay.platform.push({
      run: track.run, artifact: track.artifact_id, generation: track.generation, revision: track.revision,
      referenceId: track.referenceId, identity: track.identity,
      read: await platformRead(network, origin, track, peak),
    });
  }
  const hub = frame.getByTestId('washu-genome-viewer');
  await hub.screenshot({ path: join(work, 'overlay-three.png') });

  await secondArticle.getByRole('button', { name: 'Remove run ' + secondId }).click();
  await expect(secondArticle).toHaveCount(0);
  if (waitForRemovalDraw) await expectDrawing(manifest.window, remaining, [primaryUrl], manifest.local_bed);
  const nextRegion = 'chr1:5200-6400';
  await region.fill(nextRegion);
  evidence.removal = await expectDrawing(nextRegion, remaining, [primaryUrl], manifest.local_bed);
  await hub.screenshot({ path: join(work, 'overlay-after-remove.png') });

  await addRun(foreignId);
  const foreign = article(foreignId);
  const conflict = foreign.getByTestId('reference-conflict');
  await expect(conflict).toContainText('mm10');
  await expect(foreign.getByRole('checkbox')).toBeDisabled();
  await foreign.screenshot({ path: join(work, 'overlay-foreign.png') });
  const foreignReads = network.requests.filter(url => url.pathname.includes('/runs/' + foreignId + '/browser-tracks/'));
  assert.equal(foreignReads.length, 0);
  evidence.foreign = { run: foreignId, assembly: 'mm10', reason: await conflict.innerText(), trackRequests: foreignReads.length };
  await expectDrawing(nextRegion, remaining, [primaryUrl], manifest.local_bed);
  await foreign.getByRole('button', { name: 'Remove run ' + foreignId }).click();
  await expect(foreign).toHaveCount(0);

  await addRun(secondId);
  await secondBox.check();
  await expectDrawing(nextRegion, allPeaks, bothUrls, manifest.local_bed);

  const changed = await api.post(origin + '/__washu_overlay__/invalidate', {
    headers: { 'x-washu-overlay-token': manifest.token }, data: { run_id: secondId },
  });
  assert.equal(changed.status(), 200);
  const replacement = await changed.json();
  assert.equal(replacement.run_id, secondId);
  assert.notEqual(replacement.generation, second.generation);
  assert.notEqual(replacement.revision, second.revision);
  const readsBeforeExpiry = network.reads.length;
  // Navigate while the stale track is still selected. The page must discover
  // the conflict on its own versioned byte route and keep the other sources usable.
  const expiredWindow = 'chr1:5300-6350';
  const invalidRead = page.waitForResponse(r => r.url() === bothUrls[1] && r.status() === 409);
  await region.fill(expiredWindow);
  const stale = await invalidRead;
  const reason = await stale.text();
  assert(reason.includes('generation_changed') && !reason.includes(workspace));
  const alert = frame.getByRole('alert').filter({ hasText: secondId });
  await expect(alert).toContainText('generation_changed');
  const failureMessage = await alert.innerText();
  const kept = await api.get(primaryUrl, { headers: { Range: 'bytes=0-3' } });
  assert.equal(kept.status(), 206);
  assert.equal((await kept.body()).length, 4);
  const afterFailure = await expectDrawing(expiredWindow, remaining, [primaryUrl], manifest.local_bed);
  await hub.screenshot({ path: join(work, 'overlay-after-expiry.png') });
  evidence.expiry = {
    run: secondId, oldGeneration: second.generation, oldRevision: second.revision,
    newGeneration: replacement.generation, newRevision: replacement.revision,
    browserByteStatus: stale.status(), reason: 'generation_changed', failureMessage,
    primaryRangeStatus: kept.status(), primaryContentRange: kept.headers()['content-range'], afterFailure,
    primaryBrowserRead: await platformRead(network, origin, primary, 'peak-primary', readsBeforeExpiry),
  };
  await secondArticle.getByRole('button', { name: 'Refresh run ' + secondId }).click();
  await expect(secondArticle).toHaveAttribute('data-generation', replacement.generation);
  await expect(secondArticle.getByTestId('browser-track')).toHaveAttribute('data-revision', replacement.revision);
  await expect(secondBox).toBeEnabled();
  await expect(secondBox).not.toBeChecked();
  evidence.expiry.afterRefresh = await expectDrawing(expiredWindow, remaining, [primaryUrl], manifest.local_bed);
}

try {
  let manifest;
  for await (const line of lines) {
    if (line.startsWith('WASHU_OVERLAY ')) { manifest = JSON.parse(line.slice('WASHU_OVERLAY '.length)); break; }
  }
  assert(manifest, 'overlay fixture exited before startup');
  await main(manifest);
} catch (error) {
  evidence.error = error.message;
  if (frame) evidence.viewerText = await frame.getByTestId('washu-genome-viewer').innerText({ timeout: 1000 }).catch(() => '');
  if (page) await page.screenshot({ path: join(work, 'failure.png') }).catch(() => {});
  console.error('Overlay evidence', work);
  throw error;
} finally {
  try {
    if (browser) await browser.close();
  } finally {
    lines.close();
    if (server.exitCode === null && server.signalCode === null) server.kill();
    await serverExit;
    if (network) {
      evidence.network = networkReceipt(network);
      if (evidence.expiry) {
        evidence.expiry.replacementByteRequests = network.requests.filter(url =>
          url.pathname.includes('/bytes/' + evidence.expiry.newGeneration + '/')).length;
      }
      if (evidence.foreign) {
        evidence.foreign.trackRequests = network.requests.filter(url =>
          url.pathname.includes('/runs/' + evidence.foreign.run + '/browser-tracks/')).length;
      }
    }
    await writeFile(join(work, 'evidence.json'), JSON.stringify(evidence, null, 2) + '\n');
  }
}
assert.equal(evidence.expiry.replacementByteRequests, 0);
assert.equal(evidence.foreign.trackRequests, 0);
assert.deepEqual(evidence.network.forbidden, []);
assert.deepEqual(evidence.network.pageErrors, []);
// At base-level zoom the pinned ruler can ask for the public hg38 sequence.
// Production CSP blocks it and its fallback. A report about any track is a failure.
const referenceOnly = new Set([
  'https://vizhub.wustl.edu/public/hg38/hg38.2bit',
  'https://epigenome.wustl.edu/cors/vizhub.wustl.edu/public/hg38/hg38.2bit',
]);
for (const blocked of evidence.network.cspReports) assert(referenceOnly.has(blocked), 'Unexpected blocked URL: ' + blocked);
console.log('WASHU_OVERLAY ok', join(work, 'evidence.json'));
