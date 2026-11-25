// Run from the repo with Node 22 and the locked Playwright Chromium.
// HELIXWEAVE_PYTHON selects a wheel with the integrated archive feature and
// freshly packaged frontend. WASHU_FIXTURE_PYTHON uses workflow/envs/deeptools.lock.
// WASHU_ARCHIVE_DIR optionally selects the parent for a fresh receipt directory.
// Run `node frontend/scripts/catalog-archive-browser.mjs`.
import assert from 'node:assert/strict';
import { execFile, spawn } from 'node:child_process';
import { mkdir, mkdtemp, rename, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { createInterface } from 'node:readline';
import { fileURLToPath } from 'node:url';
import { promisify } from 'node:util';
import { chromium, expect, request as apiRequest } from '@playwright/test';

const root = fileURLToPath(new URL('../../', import.meta.url));
const python = process.env.HELIXWEAVE_PYTHON ?? 'python3';
const fixturePython = process.env.WASHU_FIXTURE_PYTHON ?? python;
for (const key of Object.keys(process.env)) {
  if (/proxy$/i.test(key)) delete process.env[key];
}
const env = Object.fromEntries(
  Object.entries(process.env).filter(([key]) => !/^(ENCODE_PIPELINE_|HELIXWEAVE_)/.test(key)),
);
const exec = promisify(execFile);
const size = 17_995_618;
const rangeBudget = 4 * 1024 * 1024;
const targetWindow = {
  region: 'chr7:55100000-55110000',
  chr: 'chr7',
  left: 55_104_400,
  right: 55_105_900,
};
const values = { 1: 4.5, 2: 9 };
const trackName = 'signal.bw';
const archiveContent = /\/archived-tracks\/[^/]+\/content(?:\.tbi)?$/;
const loginPath = '/api/v1/auth/login';

async function catalogManifest(server) {
  const lines = createInterface({ input: server.stdout, signal: AbortSignal.timeout(120_000) });
  for await (const line of lines) {
    if (line.startsWith('CATALOG_PAGE_READY ')) {
      return JSON.parse(line.slice('CATALOG_PAGE_READY '.length));
    }
  }
  throw new Error('The packaged archive fixture did not become ready. Check the API and CLI output.');
}

// Read live drawing props from the pinned WashU GenomeHub. The ruler must
// agree with the requested viewport and the measured plateau with the planted
// genomic interval. The caller compares the observed value to its version.
function renderedSignal(target) {
  const host = document.querySelector('[data-testid="washu-genome-viewer"]');
  if (!host?.querySelector('canvas')) return null;
  const key = Object.keys(host).find((candidate) => candidate.startsWith('__reactFiber'));
  let root = host[key];
  if (!root) return null;
  while (root.return) root = root.return;

  function* subtree(node) {
    if (!node) return;
    yield node;
    for (let child = node.child; child; child = child.sibling) yield* subtree(child);
  }

  // The DOM fiber may belong to the previous render after a viewport change.
  let viewer;
  for (const node of subtree(root.stateNode.current)) {
    if (node.stateNode === host) {
      viewer = node;
      break;
    }
  }
  let viewport;
  let signal;
  for (const node of subtree(viewer)) {
    const { trackModel, selectedRegion, windowWidth, data, viewRegion } = node.memoizedProps ?? {};
    if (trackModel?.type === 'ruler' && windowWidth > 0) {
      const span = selectedRegion?.getGenomeIntervals?.()[0];
      if (span && target.region === `${span.chr}:${span.start + 1}-${span.end}`) {
        viewport = { chr: span.chr, start: span.start, end: span.end };
      }
    }
    if (!Array.isArray(data) || !data.length) continue;
    const interval = viewRegion?.getGenomeIntervals?.().find((span) => span.chr === target.chr);
    if (!interval) continue;
    const bpPerPixel = (interval.end - interval.start) / data.length;
    const middle = Math.floor(((target.left + target.right) / 2 - interval.start) / bpPerPixel);
    const value = data[middle];
    if (!Number.isFinite(value) || value <= 0) continue;
    let first = middle;
    let last = middle + 1;
    while (first > 0 && data[first - 1] === value) first--;
    while (last < data.length && data[last] === value) last++;
    const left = interval.start + first * bpPerPixel;
    const right = interval.start + last * bpPerPixel;
    if (Math.abs(left - target.left) <= 3 * bpPerPixel &&
        Math.abs(right - target.right) <= 3 * bpPerPixel) {
      signal = { value, left, right, bpPerPixel };
    }
  }
  return viewport && signal ? { viewport, ...signal } : null;
}

function viewerIdentity(url, origin, fixture) {
  return url.origin === origin &&
    url.pathname === '/viewer.html' &&
    url.searchParams.size === 3 &&
    url.searchParams.get('project') === fixture.project_id &&
    url.searchParams.get('experiment') === fixture.experiment_id &&
    url.searchParams.get('archive') === fixture.archive_version_id;
}

async function checkTrackFacts(scope, fixture) {
  await expect(scope.getByTestId('archived-track-label')).toHaveText(trackName);
  await expect(scope.getByTestId('archived-track-version-label')).toHaveText(
    `第 ${fixture.version_number} 版`,
  );
  await expect(scope.getByTestId('archived-track-source')).toHaveText(
    `登记 ${fixture.source.result_id} · 修订 ${fixture.source.result_revision_id}`,
  );
  await expect(scope.getByTestId('archived-track-assembly')).toHaveText('hg38');
  await expect(scope.getByTestId('archived-track-experiment-revision')).toHaveText(
    fixture.experiment_revision_id,
  );
}

async function readTranscript(requests, contentPath) {
  assert(requests.length > 0, 'WashU made no archive content requests');
  const reads = await Promise.all(requests.map(async ({ request, phase }) => {
    const { pathname: path, search } = new URL(request.url());
    assert.equal(path + search, contentPath, 'archive read used a different version or member');
    const response = await request.response();
    assert(response, `archive read failed for ${path}`);
    const method = request.method();
    const range = await request.headerValue('range');
    const hasSession = /\bhelixweave_session=/.test(await request.headerValue('cookie') ?? '');
    const headers = response.headers();
    const { responseBodySize: bytes } = await request.sizes();
    assert(hasSession, 'the WashU worker did not send the member session');
    assert.equal(headers['accept-ranges'], 'bytes');
    assert.match(headers['cache-control'] ?? '', /\bprivate\b/);
    assert.match(headers['cache-control'] ?? '', /\bno-store\b/);
    if (method === 'HEAD' && range === null) {
      assert.equal(response.status(), 200);
      assert.equal(Number(headers['content-length']), size);
      assert.equal(bytes, 0);
    } else {
      const match = /^bytes=(\d+)-(\d+)$/.exec(range ?? '');
      assert(match, `expected a single bounded Range, got ${range}`);
      const start = Number(match[1]);
      const end = Math.min(Number(match[2]), size - 1);
      assert(start <= end && end - start + 1 <= rangeBudget, `range is too large: ${range}`);
      assert.equal(response.status(), 206);
      assert.equal(headers['content-range'], `bytes ${start}-${end}/${size}`);
      assert.equal(Number(headers['content-length']), end - start + 1);
      assert.equal(bytes, method === 'GET' ? end - start + 1 : 0);
    }
    return {
      phase, method, path, range, status: response.status(),
      content_range: headers['content-range'] ?? null,
      bytes,
    };
  }));
  const bytes = reads.reduce((sum, read) => sum + read.bytes, 0);
  assert(bytes > 0 && bytes <= rangeBudget, `WashU transferred ${bytes} bytes`);
  return { reads, bytes };
}

async function revokeMemberSessions(manifest, userId) {
  const admin = await apiRequest.newContext({ baseURL: manifest.api_base });
  try {
    const login = await admin.post(loginPath, {
      data: { username: manifest.admin_username, password: manifest.admin_password },
    });
    assert.equal(login.status(), 200);
    const { cookies } = await admin.storageState();
    const csrf = cookies.find((cookie) => cookie.name === 'helixweave_csrf');
    assert(csrf, 'the administrator login did not issue a CSRF cookie');
    const revoked = await admin.post(`/api/v1/auth/accounts/${userId}/sessions/revoke`, {
      headers: { 'X-CSRF-Token': csrf.value },
    });
    assert.equal(revoked.status(), 200);
    const { revoked_count: count } = await revoked.json();
    assert(count > 0, 'no member session was revoked');
    return count;
  } finally {
    await admin.dispose();
  }
}

async function readDenial(page, contentPath, status, viewerCode) {
  const observed = await page.evaluate(async (path) => {
    const response = await fetch(path, {
      headers: { Range: 'bytes=0-0' },
      credentials: 'same-origin',
      cache: 'no-store',
    });
    const payload = await response.json().catch(() => null);
    return { status: response.status, issue: payload?.issues?.[0]?.code ?? null };
  }, contentPath);
  assert.equal(observed.status, status);
  await page.reload();
  await expect(page.getByTestId('archived-track-error')).toHaveAttribute(
    'data-error-code', viewerCode,
  );
  await expect(page.getByTestId('washu-genome-viewer')).toHaveCount(0);
  return { ...observed, viewer_code: viewerCode };
}

async function openVersion(browser, manifest, fixture, out, work) {
  const origin = manifest.api_base;
  const experimentPath = `/api/v1/catalog/projects/${manifest.project_id}/experiments/${manifest.experiment_id}`;
  const descriptorPath = `${experimentPath}/archived-tracks/${fixture.archive_version_id}`;
  const context = await browser.newContext({
    baseURL: origin,
    viewport: { width: 1440, height: 1000 },
    serviceWorkers: 'block',
  });
  context.setDefaultTimeout(20_000);
  const external = [];
  const writes = [];
  const requests = [];
  let chooserEvents = 0;
  let logins = 0;
  let phase = 'open';
  try {
    context.on('page', (page) => page.on('filechooser', () => chooserEvents++));
    context.on('request', (request) => {
      if (archiveContent.test(new URL(request.url()).pathname)) requests.push({ request, phase });
    });
    await context.route('**/*', async (route) => {
      const request = route.request();
      const url = new URL(request.url());
      const method = request.method();
      if (['blob:', 'data:'].includes(url.protocol)) return route.continue();
      if (url.origin !== origin) {
        external.push(url.href);
        return route.abort();
      }
      if (method === 'POST' && url.pathname === loginPath && logins === 0) {
        logins++;
      } else if (!['GET', 'HEAD'].includes(method) || request.postDataBuffer() !== null) {
        writes.push(`${method} ${url.pathname}`);
        return route.abort();
      }
      return route.continue();
    });

    const catalog = await context.newPage();
    const home = await catalog.goto('/');
    assert.equal(home.status(), 200);
    assert.match(home.headers()['content-type'] ?? '', /^text\/html\b/);
    await catalog.locator('#login-username').fill(manifest.member_username);
    await catalog.locator('#login-password').fill(manifest.member_password);
    const [login] = await Promise.all([
      catalog.waitForResponse((response) => new URL(response.url()).pathname === loginPath),
      catalog.getByRole('button', { name: 'Sign in', exact: true }).click(),
    ]);
    assert.equal(login.status(), 200);
    const { principal } = await login.json();
    assert.equal(principal.username, manifest.member_username);
    assert.equal(principal.role, 'member');
    await catalog.getByRole('navigation', { name: 'Primary' })
      .getByRole('link', { name: 'Catalog', exact: true }).click();
    const search = catalog.getByRole('form', { name: 'Find experiments' });
    await search.getByLabel('Project').selectOption(manifest.project_id);
    await search.getByRole('button', { name: 'Find experiments' }).click();
    await catalog.getByRole('list', { name: 'Experiments' })
      .getByRole('button', { name: manifest.current_name }).click();
    await catalog.getByRole('region', { name: 'Experiment revision' })
      .getByRole('heading', { name: manifest.current_name, exact: true }).waitFor();

    const row = catalog.getByTestId('archived-track-version').and(
      catalog.locator(`[data-archive-version-id="${fixture.archive_version_id}"]`),
    );
    await checkTrackFacts(row, fixture);
    const action = row.getByRole('link', {
      name: `在 WashU 查看 ${trackName} 第 ${fixture.version_number} 版`,
      exact: true,
    });
    const href = await action.getAttribute('href');
    assert(viewerIdentity(new URL(href, origin), origin, fixture), 'wrong catalog viewer link');
    assert.equal(requests.length, 0, 'the catalog loaded bytes before a version was opened');
    const openedAt = performance.now();
    const [document, descriptor] = await Promise.all([
      context.waitForEvent('response', {
        predicate: (response) => response.request().isNavigationRequest() &&
          viewerIdentity(new URL(response.url()), origin, fixture),
      }),
      context.waitForEvent('response', {
        predicate: (response) => response.request().method() === 'GET' &&
          response.url() === origin + descriptorPath,
      }),
      action.click(),
    ]);
    assert.equal(document.status(), 200);
    assert.match(document.headers()['content-type'] ?? '', /^text\/html\b/);
    assert.equal(descriptor.status(), 200);
    const selected = await descriptor.json();
    for (const field of [
      'archive_version_id', 'project_id', 'experiment_id', 'experiment_revision_id', 'version_number',
    ]) {
      assert.equal(selected[field], fixture[field], field);
    }
    assert.deepEqual(selected.source, fixture.source);
    assert.equal(selected.assembly, 'hg38');
    assert.equal(selected.display_name, trackName);
    assert.equal(selected.track_type, 'bigwig');
    assert.equal(selected.has_index, false);
    const viewer = catalog;
    await viewer.waitForURL((url) => viewerIdentity(url, origin, fixture));
    await checkTrackFacts(viewer.getByTestId('archived-track-viewer'), fixture);
    const region = viewer.getByRole('textbox', { name: 'Genomic region' });
    await region.fill(targetWindow.region);
    phase = 'target-window';
    await region.press('Enter');
    const rendered = await viewer.waitForFunction(renderedSignal, targetWindow, { timeout: 30_000 });
    const elapsedMs = Math.round(performance.now() - openedAt);
    const geometry = await rendered.jsonValue();
    const expectedValue = values[fixture.version_number];
    assert.equal(geometry.value, expectedValue);
    await viewer.getByTestId('washu-genome-viewer').screenshot({
      path: join(out, `v${fixture.version_number}.png`),
    });

    const contentPath = `${descriptorPath}/content`;
    const transcript = await readTranscript(requests, contentPath);
    assert.equal(transcript.reads.length, requests.length, 'archive reads continued after the capture');
    const evidence = {
      version: fixture.version_number,
      archive_version_id: fixture.archive_version_id,
      experiment_revision_id: fixture.experiment_revision_id,
      source: selected.source,
      size_bytes: fixture.size_bytes,
      sha256: fixture.content_sha256,
      expected_value: expectedValue,
      open_to_render_ms: elapsedMs,
      geometry,
      transferred_bytes: transcript.bytes,
      ranges: transcript.reads,
    };
    if (fixture.version_number === 1) {
      const archiveRoot = join(work, 'archive');
      const offlineRoot = join(work, 'offline-archive');
      await rename(archiveRoot, offlineRoot);
      try {
        evidence.offline = await readDenial(viewer, contentPath, 409, 'archive_offline');
        assert.equal(evidence.offline.issue, 'archive_offline');
      } finally {
        await rename(offlineRoot, archiveRoot);
      }
      const revoked = await revokeMemberSessions(manifest, principal.user_id);
      evidence.revoked = {
        ...(await readDenial(viewer, contentPath, 401, 'archive_unauthorized')),
        sessions_revoked: revoked,
      };
      assert(requests.every(({ request }) => request.url() === origin + contentPath),
        'a denied read caused a request for another archived version');
    }
    assert.deepEqual(external, [], 'external requests');
    assert.deepEqual(writes, [], 'unexpected writes or uploads');
    assert.equal(chooserEvents, 0, 'file chooser events');
    assert.equal(logins, 1);
    return {
      ...evidence,
      uploads: writes.length,
      chooser_events: chooserEvents,
      external_requests: external.length,
    };
  } finally {
    await context.close();
  }
}

const outputParent = resolve(process.env.WASHU_ARCHIVE_DIR ?? tmpdir());
await mkdir(outputParent, { recursive: true });
const out = await mkdtemp(join(outputParent, 'washu-archive-'));
const work = await mkdtemp(join(tmpdir(), 'helixweave-archive-fixture-'));
let server;
let serverClosed;
let browser;
try {
  await exec(python, [
    '-I', '-S', join(root, 'scripts/checkout_bootstrap.py'),
    '--repository-root', root, 'installed-artifact',
  ], { cwd: work, env, timeout: 20_000 });
  await exec(fixturePython, [
    '-I', join(root, 'test/browser/archive_track_fixtures.py'), '--dest', join(work, 'pool'),
  ], { env, timeout: 45_000 });
  server = spawn(python, [
    '-I', join(root, 'test/browser/catalog_page_server.py'), '--archive-work', work,
  ], { cwd: work, env, stdio: ['ignore', 'pipe', 'inherit'] });
  server.on('error', (error) => process.stderr.write(`${error.message}\n`));
  serverClosed = new Promise((done) => server.once('close', done));
  const manifest = await catalogManifest(server);
  const fixtures = manifest.archive_versions;
  assert.deepEqual(fixtures.map((fixture) => fixture.version_number), [1, 2]);
  assert(fixtures.every((fixture) => fixture.size_bytes === size));
  assert.notEqual(fixtures[0].content_sha256, fixtures[1].content_sha256);
  assert.notEqual(fixtures[0].archive_version_id, fixtures[1].archive_version_id);
  browser = await chromium.launch({ args: ['--no-proxy-server'] });
  const evidence = [];
  for (const fixture of fixtures.toReversed()) {
    evidence.push(await openVersion(browser, manifest, fixture, out, work));
  }
  const receipt = {
    chromium: browser.version(),
    frontend_identity: manifest.frontend_identity,
    region: targetWindow.region,
    range_budget: rangeBudget,
    account_role: 'member',
    retired_sources: true,
    evidence,
  };
  await writeFile(join(out, 'receipt.json'), JSON.stringify(receipt, null, 2) + '\n');
  console.log(JSON.stringify({ receipt: join(out, 'receipt.json'), evidence }, null, 2));
} finally {
  try {
    await browser?.close();
  } finally {
    if (server && server.exitCode === null && server.signalCode === null) server.kill('SIGTERM');
    await serverClosed;
    await rm(work, { recursive: true, force: true });
  }
}
