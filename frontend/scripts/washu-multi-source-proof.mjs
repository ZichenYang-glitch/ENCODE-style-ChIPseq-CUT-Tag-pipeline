// Default proof: two hg38 runs, the undeclared local peak, and an mm10 run that is refused.
// Mouse: WASHU_OVERLAY_ASSEMBLY=mm10|GRCm38 and optional WASHU_OVERLAY_FOREIGN=hg38|mm39.
// The local peak then declares mm10. Version, permission, and removal checks stay.
import { spawn } from 'node:child_process';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import net from 'node:net';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import zlib from 'node:zlib';
import { chromium } from '@playwright/test';

const root = fileURLToPath(new URL('../..', import.meta.url));
const python =
  process.env.HELIXWEAVE_PYTHON ?? join(root, '.local/envs/ci-fast/bin/python');
const work = process.env.WASHU_OVERLAY_DIR
  ? (mkdirSync(process.env.WASHU_OVERLAY_DIR, { recursive: true }), process.env.WASHU_OVERLAY_DIR)
  : mkdtempSync(join(tmpdir(), 'helixweave-washu-overlay-'));
const requestedAssembly = process.env.WASHU_OVERLAY_ASSEMBLY ?? '';
const requestedForeign = process.env.WASHU_OVERLAY_FOREIGN ?? '';
const evidencePath = join(work, 'evidence.json');
const enforcedPolicy = "connect-src 'self'; report-uri /csp-report";
const navigatedRegion = 'chr1:5200-6400';
const sessionCookie = 'a'.repeat(43);
const csrfCookie = 'b'.repeat(43);
const children = [];
const reads = [];
const posts = [];
const external = [];
const consoleMessages = [];
let browser;
let page;
let port = 0;
let manifest;
let workspace = '';
let stage = 'open';

const evidence = {
  source: null,
  window: null,
  api: null,
  protection: null,
  tracks: null,
  removal: null,
  expiry: null,
  foreign: null,
  network: null,
};

function assemblyKey(value) {
  const key = String(value ?? '').trim().toLowerCase();
  if (/^(?:hg38|grch38)(?:\.p[0-9]+)?$/.test(key)) return 'hg38';
  if (key === 'grcm38') return 'mm10';
  if (key === 'grcm39') return 'mm39';
  return key;
}

function localBedPath(record) {
  if (record.local_bed.startsWith('/')) return record.local_bed;
  return join(root, 'test/fixtures/browser', record.local_bed);
}

function bedDeclares(text, declaration) {
  const escaped = declaration.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  return (
    new RegExp(`^browser genome ${escaped}$`, 'm').test(text) ||
    new RegExp(`(?:^|\\s)(?:db|genome)=["']?${escaped}(?=$|\\s)`, 'im').test(text)
  );
}

function log(line) {
  writeFileSync(join(work, 'proof.log'), `${line}\n`, { flag: 'a' });
}

function writeEvidence() {
  writeFileSync(evidencePath, `${JSON.stringify(evidence, null, 2)}\n`);
}

function start(command, args) {
  const child = spawn(command, args, {
    cwd: root,
    env: { ...process.env, PYTHONPATH: join(root, 'src') },
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

function portFree(target) {
  return new Promise((resolve) => {
    const server = net.createServer();
    server.once('error', () => resolve(false));
    server.listen(target, '127.0.0.1', () => {
      server.close(() => resolve(true));
    });
  });
}

async function waitUntilListening() {
  const started = Date.now();
  let lastError = 'fixture server did not accept connections';
  while (Date.now() - started < 20_000) {
    try {
      await fetch(`http://127.0.0.1:${port}/viewer.html`);
      return;
    } catch (error) {
      lastError = error instanceof Error ? error.message : String(error);
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
  }
  throw new Error(lastError);
}

async function stopAll() {
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
  if (port && !(await portFree(port))) {
    throw new Error(`proof port still listening: ${port}`);
  }
}

function cookieHeader() {
  return `helixweave_session=${sessionCookie}; helixweave_csrf=${csrfCookie}`;
}

async function apiGet(path, { cookie = true } = {}) {
  const headers = {};
  if (cookie) headers.cookie = cookieHeader();
  const response = await fetch(`http://127.0.0.1:${port}${path}`, { headers });
  const body = Buffer.from(await response.arrayBuffer());
  return { status: response.status, text: body.toString('utf8'), body };
}

function bytePath(runId, generation, revision, artifactId = manifest.artifact_id) {
  return `/api/v1/runs/${runId}/browser-tracks/${artifactId}/bytes/${generation}/${revision}`;
}

function assertClosed(response, label) {
  if (response.text.includes(workspace)) {
    throw new Error(`${label} leaked the workspace path`);
  }
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

function parseTrackUrl(url) {
  const match = url.match(
    /\/runs\/([^/]+)\/browser-tracks\/([^/]+)\/bytes\/([^/]+)\/([^/?]+)/,
  );
  if (!match) return null;
  return {
    runId: match[1],
    artifactId: match[2],
    generation: match[3],
    revision: decodeURIComponent(match[4].replace(/\.tbi$/, '')),
    index: match[4].includes('.tbi') || url.includes('.tbi'),
  };
}

function isExternal(url) {
  const parsed = new URL(url);
  if (parsed.protocol === 'blob:' || parsed.protocol === 'data:') return false;
  return parsed.hostname !== '127.0.0.1' && parsed.hostname !== 'localhost';
}

function runArticle(frame, runId) {
  return frame.locator(`[data-testid="browser-run"][data-run-id="${runId}"]`);
}

async function readHub(frame) {
  return frame.evaluate(() => {
    const root = document.querySelector('[data-testid="washu-genome-viewer"]');
    if (!(root instanceof HTMLElement)) return null;
    const key = Object.getOwnPropertyNames(root).find((name) => name.startsWith('__reactFiber'));
    const fiber = key ? root[key] : null;
    let models = [];
    const regions = [];
    const drawn = [];
    const stack = fiber ? [fiber] : [];
    let guard = 0;
    while (stack.length && guard < 12000) {
      guard += 1;
      const node = stack.pop();
      const props = node?.memoizedProps;
      if (props && Array.isArray(props.tracks)) {
        const next = props.tracks.map((track) => ({
          type: typeof track?.type === 'string' ? track.type : '',
          name: typeof track?.name === 'string' ? track.name : '',
          url: typeof track?.url === 'string' ? track.url : '',
          file: Boolean(track?.fileObj || (Array.isArray(track?.files) && track.files.length)),
        }));
        if (next.some((track) => track.type === 'ruler') && next.length >= models.length) {
          models = next;
        }
      }
      const region = props?.viewRegion;
      if (region && typeof region.getWidth === 'function') {
        let start = null;
        let end = null;
        try {
          const coords = region.getContextCoordinates();
          start = coords.start;
          end = coords.end;
        } catch {
          start = null;
        }
        regions.push({ start, end });
      }
      if (props) {
        for (const bucket of [props.data, props.items, props.features, props.placements, props.arrangedData]) {
          if (!Array.isArray(bucket)) continue;
          for (const item of bucket) {
            const feature = item?.feature;
            if (feature?.chr === 'chr1' && item.placedFeatures?.length && item.xSpan?.end > item.xSpan?.start) {
              drawn.push(`${feature.chr}:${feature.start}-${feature.end} ${feature[3] ?? ''}`);
            }
          }
        }
      }
      if (node?.child) stack.push(node.child);
      if (node?.sibling) stack.push(node.sibling);
    }
    return {
      viewers: document.querySelectorAll('[data-testid="washu-genome-viewer"]').length,
      models,
      regions,
      drawn,
    };
  });
}

function showsRegion(regions, displayRegion) {
  const [, start, end] = displayRegion.match(/^chr1:(\d+)-(\d+)$/);
  return regions.some(
    (region) =>
      region.start === Number(start) - 1 && region.end === Number(end),
  );
}

async function waitForHub(frame, ready, timeoutMs) {
  const started = Date.now();
  let last = null;
  while (Date.now() - started < timeoutMs) {
    last = await readHub(frame);
    if (last && ready(last)) return last;
    await page.waitForTimeout(200);
  }
  throw new Error(`genome hub did not reach the expected tracks\n${JSON.stringify(last)}`);
}

function dataReads(runId, since = 0) {
  return reads.slice(since).filter((read) => {
    const parsed = parseTrackUrl(read.url);
    return (
      parsed &&
      parsed.runId === runId &&
      parsed.artifactId === manifest.artifact_id &&
      !parsed.index &&
      !read.preflight &&
      (read.status === 200 || read.status === 206) &&
      read.body.length > 1
    );
  });
}

async function protectionChecks() {
  const primarySession = await apiGet(`/api/v1/runs/${manifest.primary_run_id}/browser-session`);
  const secondSession = await apiGet(`/api/v1/runs/${manifest.second_run_id}/browser-session`);
  const foreignSession = await apiGet(`/api/v1/runs/${manifest.foreign_run_id}/browser-session`);
  for (const [label, response] of [
    ['primary session', primarySession],
    ['second session', secondSession],
    ['foreign session', foreignSession],
  ]) {
    if (response.status !== 200) throw new Error(`${label} returned ${response.status}`);
    assertClosed(response, label);
  }
  const primary = JSON.parse(primarySession.text);
  const second = JSON.parse(secondSession.text);
  const foreign = JSON.parse(foreignSession.text);
  if (
    primary.reference.viewer_scaffold !== manifest.scaffold ||
    second.reference.viewer_scaffold !== manifest.scaffold
  ) {
    throw new Error(`the two authorized runs are not both on the ${manifest.scaffold} scaffold`);
  }
  if (
    primary.reference.assembly_label !== manifest.assembly ||
    second.reference.assembly_label !== manifest.assembly
  ) {
    throw new Error(`the two authorized runs do not both declare ${manifest.assembly}`);
  }
  if (
    !primary.reference.identity_sha256 ||
    primary.reference.identity_sha256 === second.reference.identity_sha256
  ) {
    throw new Error('the second run does not have its own reference identity');
  }
  if (
    foreign.reference.viewer_scaffold !== manifest.foreign_scaffold ||
    foreign.reference.assembly_label !== manifest.foreign_assembly
  ) {
    throw new Error(
      `foreign run was ${foreign.reference.assembly_label}/${foreign.reference.viewer_scaffold}, not ${manifest.foreign_assembly}/${manifest.foreign_scaffold}`,
    );
  }
  if (manifest.foreign_scaffold === manifest.scaffold) {
    throw new Error(`${manifest.foreign_assembly} shares the ${manifest.scaffold} window`);
  }
  const history = await apiGet('/api/v1/runs?limit=100');
  assertClosed(history, 'run history');
  if (history.status !== 200) throw new Error(`run history returned ${history.status}`);
  const listed = new Set(JSON.parse(history.text).runs.map((run) => run.run_id));
  for (const runId of [manifest.primary_run_id, manifest.second_run_id, manifest.foreign_run_id]) {
    if (!listed.has(runId)) throw new Error(`authorized run ${runId} was not listed`);
  }
  for (const [label, body, generation, revision] of [
    ['primary', primary, manifest.primary_generation, manifest.primary_revision],
    ['second', second, manifest.second_generation, manifest.second_revision],
  ]) {
    const track = body.tracks.find((item) => item.artifact_id === manifest.artifact_id);
    if (!track || track.preparation !== 'ready' || track.revision !== revision) {
      throw new Error(`${label} peak was not ready at the seeded revision`);
    }
    if (body.artifact_generation !== generation) {
      throw new Error(`${label} generation did not match the seed`);
    }
  }
  const foreignTrack = foreign.tracks.find((item) => item.artifact_id === manifest.artifact_id);
  if (!foreignTrack) throw new Error('foreign run did not publish the shared artifact id');

  const anonymous = await apiGet(
    bytePath(manifest.primary_run_id, manifest.primary_generation, manifest.primary_revision),
    { cookie: false },
  );
  assertClosed(anonymous, 'anonymous byte read');
  if (anonymous.status !== 401 || anonymous.text.includes('peak-primary')) {
    throw new Error(`anonymous byte read returned ${anonymous.status}`);
  }
  const missing = await apiGet(
    bytePath(
      manifest.primary_run_id,
      manifest.primary_generation,
      manifest.primary_revision,
      'missing-artifact',
    ),
  );
  assertClosed(missing, 'missing artifact');
  if (missing.status !== 404 || !missing.text.includes('RUN_ARTIFACT_NOT_FOUND')) {
    throw new Error(`missing artifact returned ${missing.status} ${missing.text.slice(0, 180)}`);
  }
  const staleRevision = await apiGet(
    bytePath(
      manifest.primary_run_id,
      manifest.primary_generation,
      `artifactrev-${'d'.repeat(64)}`,
    ),
  );
  assertClosed(staleRevision, 'stale revision');
  if (staleRevision.status !== 409 || !staleRevision.text.includes('revision_changed')) {
    throw new Error(`stale revision returned ${staleRevision.status} ${staleRevision.text.slice(0, 180)}`);
  }
  const staleGeneration = await apiGet(
    bytePath(
      manifest.primary_run_id,
      `artifactgen-${'e'.repeat(64)}`,
      manifest.primary_revision,
    ),
  );
  assertClosed(staleGeneration, 'stale generation');
  if (staleGeneration.status !== 409 || !staleGeneration.text.includes('generation')) {
    throw new Error(`stale generation returned ${staleGeneration.status}`);
  }
  evidence.api = {
    primary_identity: primary.reference.identity_sha256,
    second_identity: second.reference.identity_sha256,
    foreign_identity: foreign.reference.identity_sha256,
    primary_revision_id: primary.reference.revision_id,
    second_revision_id: second.reference.revision_id,
    scaffolds: {
      primary: primary.reference.viewer_scaffold,
      second: second.reference.viewer_scaffold,
      foreign: foreign.reference.viewer_scaffold,
    },
    labels: {
      primary: primary.reference.assembly_label,
      second: second.reference.assembly_label,
      foreign: foreign.reference.assembly_label,
    },
    artifact_id: manifest.artifact_id,
    primary_generation: manifest.primary_generation,
    second_generation: manifest.second_generation,
  };
  evidence.protection = {
    anonymous: anonymous.status,
    missing_artifact: missing.status,
    revision_changed: staleRevision.status,
    generation_rejected: staleGeneration.status,
    workspace_leaked: false,
  };
}

function inflatedText(runId) {
  return Buffer.concat(dataReads(runId).map((read) => inflateBgzf(read.body))).toString('utf8');
}

function sameOriginRead(runId, generation, revision) {
  const hits = reads.filter((read) => {
    const parsed = parseTrackUrl(read.url);
    if (!parsed || parsed.runId !== runId) return false;
    const url = new URL(read.url);
    return (
      url.hostname === '127.0.0.1' &&
      url.pathname.startsWith('/api/') &&
      parsed.generation === generation &&
      parsed.revision === revision &&
      (read.status === 200 || read.status === 206)
    );
  });
  const data = hits.find((read) => !parseTrackUrl(read.url).index && !read.preflight && read.body.length > 1);
  const index = hits.find((read) => parseTrackUrl(read.url).index && !read.preflight && read.body.length > 1);
  if (!data || !index) return null;
  return { data: data.url, index: index.url };
}

async function main() {
  port = await reserve(18940);
  workspace = join(work, 'workspace');
  mkdirSync(workspace, { recursive: true });
  const serverArgs = [
    'test/browser/washu_overlay_server.py',
    '--port',
    String(port),
    '--database',
    join(work, 'overlay.sqlite'),
    '--workspace',
    workspace,
  ];
  if (requestedAssembly) serverArgs.push('--assembly', requestedAssembly);
  if (requestedForeign) serverArgs.push('--foreign', requestedForeign);
  const server = start(python, serverArgs);
  server.stderr.on('data', (chunk) => {
    writeFileSync(join(work, 'server.log'), chunk.toString(), { flag: 'a' });
  });
  const ready = await waitFor(server.stdout, /WASHU_OVERLAY (\{.*\})/, 90_000);
  manifest = JSON.parse(ready[1]);
  if (requestedAssembly && assemblyKey(manifest.assembly) !== assemblyKey(requestedAssembly)) {
    throw new Error(`server assembly ${manifest.assembly} does not match ${requestedAssembly}`);
  }
  if (requestedForeign && assemblyKey(manifest.foreign_assembly) !== assemblyKey(requestedForeign)) {
    throw new Error(`server foreign assembly ${manifest.foreign_assembly} does not match ${requestedForeign}`);
  }
  const localBed = localBedPath(manifest);
  if (!existsSync(localBed)) throw new Error(`local BED fixture is missing: ${localBed}`);
  const localText = readFileSync(localBed, 'utf8');
  if (!localText.includes('peak-local') || !localText.includes('5800\t6000')) {
    throw new Error('local BED fixture is not the synthetic peak-local interval');
  }
  if (manifest.local_declaration && !bedDeclares(localText, manifest.local_declaration)) {
    throw new Error(`local peak does not declare ${manifest.local_declaration}`);
  }
  if (!manifest.local_declaration && manifest.scaffold !== 'hg38') {
    throw new Error('a mouse window requires a local peak that declares mm10');
  }
  evidence.source = {
    assembly: manifest.assembly,
    scaffold: manifest.scaffold,
    organism: manifest.organism,
    foreign_assembly: manifest.foreign_assembly,
    foreign_scaffold: manifest.foreign_scaffold,
    local_declaration: manifest.local_declaration,
    local_bed: manifest.local_bed,
    window: manifest.window,
  };
  evidence.window = manifest.window;
  await waitUntilListening();
  await protectionChecks();
  writeEvidence();

  browser = await chromium.launch({ headless: true });
  const context = await browser.newContext();
  await context.addCookies([
    {
      name: 'helixweave_session',
      value: sessionCookie,
      domain: '127.0.0.1',
      path: '/',
      httpOnly: true,
      sameSite: 'Lax',
    },
    {
      name: 'helixweave_csrf',
      value: csrfCookie,
      domain: '127.0.0.1',
      path: '/',
      httpOnly: false,
      sameSite: 'Lax',
    },
  ]);
  page = await context.newPage();
  page.on('console', (message) => {
    consoleMessages.push(`${message.type()}: [${stage}] ${message.text()}`);
  });
  page.on('pageerror', (error) => {
    consoleMessages.push(`pageerror: ${error.message}`);
  });
  page.on('request', (request) => {
    const url = request.url();
    if (isExternal(url)) external.push(`${request.method()} ${url}`);
    if (request.method() !== 'GET' && request.method() !== 'HEAD') {
      posts.push({ method: request.method(), url, body: request.postData() ?? '' });
    }
  });
  page.on('response', (response) => {
    const url = response.url();
    if (!url.includes('/bytes/')) return;
    const range = response.request().headers().range ?? '';
    const entry = {
      url,
      range,
      status: response.status(),
      body: Buffer.alloc(0),
      preflight: range === 'bytes=0-0' || range === 'bytes=0-1',
    };
    reads.push(entry);
    void response.body().then((buf) => {
      entry.body = Buffer.from(buf);
    }).catch(() => {});
  });

  const navigation = await page.goto(
    `http://127.0.0.1:${port}/runs/${manifest.primary_run_id}?view=artifacts`,
  );
  const documentPolicy = navigation?.headers()['content-security-policy'] ?? '';
  await page.getByRole('button', { name: 'View tracks in genome browser' }).click();
  await page.getByTitle('Genome browser').waitFor();
  let viewerFrame;
  for (let attempt = 0; attempt < 50; attempt += 1) {
    viewerFrame = page.frames().find((frame) => frame.url().includes('/viewer.html'));
    if (viewerFrame) break;
    await page.waitForTimeout(100);
  }
  if (!viewerFrame) throw new Error('viewer frame did not load');
  const viewerResponse = await page.request.get(viewerFrame.url());
  const viewerPolicy = viewerResponse.headers()['content-security-policy'] ?? '';
  const viewerFrameOptions = viewerResponse.headers()['x-frame-options'] ?? '';
  evidence.protection = {
    ...evidence.protection,
    document_policy: documentPolicy,
    viewer_policy: viewerPolicy,
    frame_options: viewerFrameOptions,
  };
  if (
    !documentPolicy.startsWith(enforcedPolicy) ||
    !viewerPolicy.startsWith(enforcedPolicy) ||
    !viewerPolicy.includes("frame-ancestors 'self'") ||
    viewerFrameOptions !== 'SAMEORIGIN'
  ) {
    throw new Error(
      `viewer protection differed from the packaged contract ${JSON.stringify(evidence.protection)}`,
    );
  }
  const overlayReady = await viewerFrame
    .locator('[data-testid="browser-run"]')
    .waitFor({ timeout: 20_000 })
    .then(() => true)
    .catch(() => false);
  if (!overlayReady) {
    throw new Error(
      'WASHU_OVERLAY ui-not-ready: the packaged page has no browser-run selector yet',
    );
  }
  const primaryBox = runArticle(viewerFrame, manifest.primary_run_id).getByRole('checkbox');
  await primaryBox.waitFor({ timeout: 20_000 });
  if ((await primaryBox.isChecked()) || reads.length !== 0) {
    throw new Error(`opened with a selected track or ${reads.length} byte reads`);
  }
  const secondBoxEarly = runArticle(viewerFrame, manifest.second_run_id);
  if (await secondBoxEarly.count()) {
    throw new Error('the second run was open before it was added');
  }

  stage = 'select primary';
  await primaryBox.check();
  await viewerFrame.locator('[data-testid="local-track-file-input"]').setInputFiles(localBed);
  const localRow = viewerFrame.locator('[data-testid="local-browser-track"][data-source="local"]');
  await localRow.waitFor({ timeout: 20_000 });
  const localAssembly = await localRow.getAttribute('data-assembly');
  if (localAssembly !== manifest.scaffold) {
    throw new Error(`local peak was placed on ${localAssembly}, not ${manifest.scaffold}`);
  }
  if (manifest.local_declaration) {
    const signature = await viewerFrame.locator('[data-testid="local-track-assembly"]').innerText();
    if (!signature.includes(manifest.local_declaration)) {
      throw new Error(`local assembly signature did not record ${manifest.local_declaration}: ${signature}`);
    }
  }
  await viewerFrame.getByRole('combobox', { name: 'Authorized runs' }).selectOption(manifest.second_run_id);
  await viewerFrame.getByRole('button', { name: 'Add run' }).click();
  const secondArticle = runArticle(viewerFrame, manifest.second_run_id);
  await secondArticle.waitFor({ timeout: 20_000 });
  await secondArticle.getByRole('checkbox').waitFor({ timeout: 20_000 });
  const secondIdentity = await secondArticle.getAttribute('data-reference-identity');
  const primaryIdentity = await runArticle(viewerFrame, manifest.primary_run_id).getAttribute(
    'data-reference-identity',
  );
  const conflict = (await secondArticle.locator('[data-testid="reference-conflict"]').allInnerTexts()).join('\n');
  if (conflict) {
    throw new Error(`same-assembly run was rejected\n${conflict}`);
  }
  if (primaryIdentity !== evidence.api.primary_identity || secondIdentity !== evidence.api.second_identity || primaryIdentity === secondIdentity) {
    throw new Error('the open runs do not show two reference identities');
  }
  const primaryScaffold = await runArticle(viewerFrame, manifest.primary_run_id).getAttribute('data-scaffold');
  const secondScaffold = await secondArticle.getAttribute('data-scaffold');
  if (primaryScaffold !== manifest.scaffold || secondScaffold !== manifest.scaffold) {
    throw new Error(`open runs are not on ${manifest.scaffold} (${primaryScaffold}, ${secondScaffold})`);
  }
  const secondBox = secondArticle.getByRole('checkbox');
  if (await secondBox.isDisabled()) throw new Error('second run track could not be chosen');
  stage = 'select second';
  await secondBox.check();

  const started = Date.now();
  while (Date.now() - started < 25_000) {
    if (dataReads(manifest.primary_run_id).length && dataReads(manifest.second_run_id).length) break;
    await page.waitForTimeout(200);
  }
  const primaryBytes = inflatedText(manifest.primary_run_id);
  const secondBytes = inflatedText(manifest.second_run_id);
  if (!primaryBytes.includes('peak-primary\t61') || !primaryBytes.includes('chr1\t5500\t5700')) {
    throw new Error(`primary peak bytes were not read\n${primaryBytes.slice(0, 200)}`);
  }
  if (!secondBytes.includes('peak-second\t73') || !secondBytes.includes('chr1\t6000\t6200')) {
    throw new Error(`second peak bytes were not read\n${secondBytes.slice(0, 200)}`);
  }
  const primaryRead = sameOriginRead(
    manifest.primary_run_id,
    manifest.primary_generation,
    manifest.primary_revision,
  );
  const secondRead = sameOriginRead(
    manifest.second_run_id,
    manifest.second_generation,
    manifest.second_revision,
  );
  if (!primaryRead || !secondRead) {
    throw new Error('platform tracks were not read from their own generation and revision');
  }
  const hub = await waitForHub(
    viewerFrame,
    (state) => {
      const urls = state.models.map((track) => track.url);
      return (
        state.viewers === 1 &&
        urls.some((url) => url.includes(`/runs/${manifest.primary_run_id}/`)) &&
        urls.some((url) => url.includes(`/runs/${manifest.second_run_id}/`)) &&
        state.models.some((track) => track.file && track.url === '') &&
        showsRegion(state.regions, manifest.window) &&
        ['chr1:5500-5700', 'chr1:5800-6000', 'chr1:6000-6200'].every(locus => state.drawn.some(text => text.includes(locus)))
      );
    },
    25_000,
  );
  if (!hub.drawn.some((text) => text.includes('5800') || text.includes('peak-local'))) {
    throw new Error('local peak was not drawn in the shared window');
  }
  evidence.tracks = {
    viewers: hub.viewers,
    models: hub.models.map((track) => ({ type: track.type, name: track.name, file: track.file, same_origin: track.url.startsWith(`http://127.0.0.1:${port}/api/`) || track.url === '' })),
    primary: {
      run_id: manifest.primary_run_id,
      generation: manifest.primary_generation,
      revision: manifest.primary_revision,
      identity: primaryIdentity,
      read: primaryRead.data,
    },
    second: {
      run_id: manifest.second_run_id,
      generation: manifest.second_generation,
      revision: manifest.second_revision,
      identity: secondIdentity,
      read: secondRead.data,
    },
    plotted_features: [...new Set(hub.drawn)],
    local_assembly: localAssembly,
    local_declaration: manifest.local_declaration,
    local_drawn: hub.drawn.filter((text) => text.includes('peak-local') || text.includes('5800')).slice(0, 4),
    regions: hub.regions.slice(0, 4),
  };
  await page.screenshot({ path: join(work, 'overlay-three.png'), fullPage: true });

  await viewerFrame.getByRole('button', { name: `Remove run ${manifest.second_run_id}` }).click();
  await secondArticle.waitFor({ state: 'detached', timeout: 10_000 });
  const regionBox = viewerFrame.getByRole('textbox', { name: 'Genomic region' });
  stage = 'navigate after removal';
  await regionBox.fill(navigatedRegion);
  await viewerFrame.waitForFunction(
    (expected) => {
      const input = document.querySelector('[aria-label="Genomic region"]');
      return input instanceof HTMLInputElement && input.value === expected;
    },
    navigatedRegion,
    { timeout: 10_000 },
  );
  const afterRemoval = await waitForHub(
    viewerFrame,
    (state) => {
      const urls = state.models.map((track) => track.url);
      return (
        state.viewers === 1 &&
        urls.some((url) => url.includes(`/runs/${manifest.primary_run_id}/`)) &&
        !urls.some((url) => url.includes(`/runs/${manifest.second_run_id}/`)) &&
        state.models.some((track) => track.file && track.url === '') &&
        showsRegion(state.regions, navigatedRegion)
      );
    },
    20_000,
  );
  evidence.removal = {
    removed: manifest.second_run_id,
    region: navigatedRegion,
    remaining: afterRemoval.models.map((track) => track.name),
    regions: afterRemoval.regions.slice(0, 4),
  };
  await page.screenshot({ path: join(work, 'overlay-after-remove.png'), fullPage: true });

  await viewerFrame.getByRole('combobox', { name: 'Authorized runs' }).selectOption(manifest.second_run_id);
  await viewerFrame.getByRole('button', { name: 'Add run' }).click();
  const restored = runArticle(viewerFrame, manifest.second_run_id);
  await restored.waitFor({ timeout: 20_000 });
  stage = 'restore second';
  await restored.getByRole('checkbox').check();
  const beforeExpiry = reads.length;
  const restoredAt = Date.now();
  while (Date.now() - restoredAt < 20_000) {
    if (dataReads(manifest.second_run_id, beforeExpiry).length) break;
    await page.waitForTimeout(200);
  }
  const liveSecond = dataReads(manifest.second_run_id, beforeExpiry)[0];
  if (!liveSecond) throw new Error('the restored second run was not read');
  const oldUrl = new URL(liveSecond.url);
  const invalidateResponse = await fetch(`http://127.0.0.1:${port}/__washu_overlay__/invalidate`, {
    method: 'POST',
    headers: {
      'content-type': 'application/json',
      'x-washu-overlay-token': manifest.token,
    },
    body: JSON.stringify({ run_id: manifest.second_run_id }),
  });
  const invalidated = await invalidateResponse.json();
  if (!invalidateResponse.ok || invalidated.generation === manifest.second_generation) {
    throw new Error(`second run generation did not change ${JSON.stringify(invalidated)}`);
  }
  const rejected = await apiGet(`${oldUrl.pathname}${oldUrl.search}`);
  assertClosed(rejected, 'expired byte read');
  if (
    rejected.status !== 409 ||
    (!rejected.text.includes('generation_changed') && !rejected.text.includes('revision_changed'))
  ) {
    throw new Error(`old second-run bytes returned ${rejected.status} ${rejected.text.slice(0, 180)}`);
  }
  const primaryStill = await apiGet(
    bytePath(manifest.primary_run_id, manifest.primary_generation, manifest.primary_revision),
  );
  if (primaryStill.status !== 206 && primaryStill.status !== 200) {
    throw new Error(`primary bytes failed after the other run changed (${primaryStill.status})`);
  }
  const expiredWindow = 'chr1:5300-6350';
  stage = 'navigate after expiry';
  await regionBox.fill(expiredWindow);
  const expiredNotice = viewerFrame.getByRole('alert').filter({ hasText: manifest.second_run_id });
  await expiredNotice.getByText(/generation_changed|revision_changed/).waitFor({ timeout: 20_000 });
  const notice = await expiredNotice.innerText();
  const afterExpiry = await waitForHub(viewerFrame, state =>
    state.models.some(track => track.url.includes(`/runs/${manifest.primary_run_id}/`)) &&
    state.models.some(track => track.file && track.url === '') &&
    !state.models.some(track => track.url.includes(`/runs/${manifest.second_run_id}/`)) &&
    state.drawn.some(item => item.startsWith('chr1:5500-5700')) &&
    state.drawn.some(item => item.startsWith('chr1:5800-6000')), 20_000);
  const newGenerationRead = reads.some((read) => {
    const parsed = parseTrackUrl(read.url);
    return (
      parsed &&
      parsed.runId === manifest.second_run_id &&
      parsed.generation === invalidated.generation &&
      (read.status === 200 || read.status === 206)
    );
  });
  if (newGenerationRead) throw new Error('the replacement generation was read without a new choice');
  const expiredUrls = afterExpiry?.models.map((track) => track.url) ?? [];
  if (
    !expiredUrls.some((url) => url.includes(`/runs/${manifest.primary_run_id}/`)) ||
    !(afterExpiry?.models ?? []).some((track) => track.file && track.url === '')
  ) {
    throw new Error(`primary or local track left the hub\n${JSON.stringify(afterExpiry?.models)}`);
  }
  evidence.expiry = {
    run_id: manifest.second_run_id,
    old_generation: manifest.second_generation,
    old_revision: manifest.second_revision,
    new_generation: invalidated.generation,
    new_revision: invalidated.revision,
    rejected_status: rejected.status,
    rejected_reason: rejected.text.includes('generation_changed') ? 'generation_changed' : 'revision_changed',
    visible_failure: notice,
    remaining_features: [...new Set(afterExpiry.drawn)],
    primary_still_readable: primaryStill.status,
    new_generation_read: false,
  };

  const readsBeforeForeign = reads.length;
  await viewerFrame.getByRole('combobox', { name: 'Authorized runs' }).selectOption(manifest.foreign_run_id);
  await viewerFrame.getByRole('button', { name: 'Add run' }).click();
  const foreignArticle = runArticle(viewerFrame, manifest.foreign_run_id);
  stage = 'add incompatible source';
  await foreignArticle.getByTestId('reference-conflict').waitFor({ timeout: 20_000 });
  const foreignConflict = (
    await foreignArticle.locator('[data-testid="reference-conflict"]').allInnerTexts()
  ).join('\n');
  if (!foreignConflict.includes(manifest.foreign_scaffold)) {
    throw new Error(`${manifest.foreign_assembly} run was not rejected\n${foreignConflict}`);
  }
  const foreignScaffold = await foreignArticle.getAttribute('data-scaffold');
  if (foreignScaffold !== manifest.foreign_scaffold) {
    throw new Error(`foreign run scaffold on the page was ${foreignScaffold}`);
  }
  const foreignBox = foreignArticle.getByRole('checkbox');
  if (await foreignBox.isEnabled()) {
    throw new Error(`${manifest.foreign_assembly} track stayed selectable`);
  }
  await page.waitForTimeout(1000);
  const foreignReads = reads.slice(readsBeforeForeign).filter((read) => read.url.includes(manifest.foreign_run_id));
  const foreignBytes = foreignReads.reduce((sum, read) => sum + read.body.length, 0);
  if (foreignReads.length || foreignBytes !== 0) {
    throw new Error(
      `${manifest.foreign_assembly} bytes were read (${foreignBytes} B)\n${foreignReads.map((read) => read.url).join('\n')}`,
    );
  }
  const stillThere = await readHub(viewerFrame);
  const stillUrls = stillThere?.models.map((track) => track.url) ?? [];
  if (
    !stillUrls.some((url) => url.includes(`/runs/${manifest.primary_run_id}/`)) ||
    !(stillThere?.models ?? []).some((track) => track.file)
  ) {
    throw new Error(`rejecting ${manifest.foreign_assembly} removed the tracks already open`);
  }
  evidence.foreign = {
    run_id: manifest.foreign_run_id,
    assembly: manifest.foreign_assembly,
    scaffold: manifest.foreign_scaffold,
    reason: foreignConflict,
    bytes_read: foreignBytes,
  };
  await page.screenshot({ path: join(work, 'overlay-foreign.png'), fullPage: true });

  const uploaded = posts.filter(
    (post) => post.body.includes('peak-local') || post.body.includes('5800\t6000'),
  );
  const controlFromPage = posts.filter((post) => post.url.includes('/__washu_overlay__/'));
  if (uploaded.length || controlFromPage.length || external.length) {
    throw new Error(
      `local content or another origin left the browser\n${JSON.stringify({ uploaded, controlFromPage, external })}`,
    );
  }
  if (reads.some((read) => isExternal(read.url) || !new URL(read.url).pathname.startsWith('/api/'))) {
    throw new Error('a track read left the same-origin artifact API');
  }
  const renderingErrors = consoleMessages.filter(message =>
    /Minified React error|Maximum update depth exceeded|Genome View Error|^pageerror:/.test(message));
  if (renderingErrors.length) throw new Error(renderingErrors.join('\n'));
  const cspBlocked = posts.filter(post => new URL(post.url).pathname === '/csp-report').map(post => {
    const blocked = JSON.parse(post.body)['csp-report']?.['blocked-uri'];
    try { return new URL(blocked).origin; } catch { return blocked ?? null; }
  });
  evidence.network = {
    external: external.length,
    csp_blocked_origins: cspBlocked,
    rendering_errors: renderingErrors,
    local_uploads: uploaded.length,
    byte_reads: reads.length,
    posts: posts.map((post) => `${post.method} ${new URL(post.url).pathname}`),
  };
  writeEvidence();
  await stopAll();
  console.log(`WASHU_OVERLAY ok ${evidencePath}`);
}

try {
  await main();
} catch (error) {
  const message = error instanceof Error ? error.message : String(error);
  evidence.error = message;
  evidence.console = consoleMessages.filter(
    (line) => line.startsWith('error') || line.startsWith('pageerror'),
  );
  writeEvidence();
  log(message);
  console.error(message);
  console.error(`proof work ${work}`);
  console.error(`evidence ${evidencePath}`);
  process.exitCode = 1;
  try {
    await stopAll();
  } catch (stopError) {
    console.error(stopError instanceof Error ? stopError.message : stopError);
    process.exitCode = 1;
  }
}
