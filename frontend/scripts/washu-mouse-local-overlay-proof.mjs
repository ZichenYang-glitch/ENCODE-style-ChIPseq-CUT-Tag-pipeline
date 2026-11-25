// Positive mm10 check from the anonymous home page. Synthetic files only.
// The drawing check is washu-local-overlay-observe.mjs. hg38 stays on
// washu-local-overlay-proof.mjs unless WASHU_ASSEMBLY=mm10.
//
// HELIXWEAVE_PYTHON=/path/to/wheel/bin/python \
// WASHU_FIXTURE_PYTHON=/path/to/chipseq/bin/python \
//   node frontend/scripts/washu-mouse-local-overlay-proof.mjs
import assert from 'node:assert/strict';
import { spawn, spawnSync } from 'node:child_process';
import { once } from 'node:events';
import { mkdir, mkdtemp, readFile, stat, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { basename, join, resolve } from 'node:path';
import { createInterface } from 'node:readline';
import { fileURLToPath } from 'node:url';
import { chromium } from '@playwright/test';
import { installFileReadMeter } from './washu-file-reads.mjs';
import { renderedWindow } from './washu-local-overlay-observe.mjs';

const root = fileURLToPath(new URL('../../', import.meta.url));
const work = process.env.WASHU_MOUSE_DIR ?? await mkdtemp(join(tmpdir(), 'washu-mouse-'));
await mkdir(work, { recursive: true });
const fixtures = resolve(process.env.WASHU_MOUSE_FIXTURES ?? join(work, 'fixtures'));
if (!process.env.WASHU_MOUSE_FIXTURES) {
  const generate = spawnSync(process.env.WASHU_FIXTURE_PYTHON ?? 'python3', [
    'test/browser/local_overlay_fixtures.py', '--assembly', 'mm10', '--dest', fixtures,
  ], { cwd: root, stdio: 'inherit' });
  if (generate.error || generate.status !== 0) throw generate.error ?? new Error('Mouse fixture generation failed');
}
const planted = JSON.parse(await readFile(join(root, 'test/fixtures/browser/local-overlay-mm10.json'), 'utf8'));
assert.equal(planted.assembly, 'mm10');
const names = planted.filenames;
const data = join(fixtures, names.longrange);
const index = data + '.tbi';
const signal = join(fixtures, names.signal);
const peak = join(fixtures, names.peaks);
const observation = JSON.parse(await readFile(join(fixtures, 'observation.json'), 'utf8'));
assert.equal(observation.assembly, 'mm10');
if (!process.env.WASHU_MOUSE_FIXTURES) assert.equal(observation.record_count, 4);

function inside(gene, interval) {
  const [low, high] = gene.zero_based_half_open;
  return interval.chrom === gene.chrom && interval.start >= low && interval.end <= high;
}
function expectedDrawing(geneName) {
  const gene = planted.genes[geneName];
  const pair = planted.longrange_pairs[geneName.toLowerCase()][0];
  const peakRow = planted.peaks.find(item => inside(gene, item));
  const signalRow = planted.signals.find(item => inside(gene, item));
  assert.equal(gene.hub_region, `${gene.chrom}:${gene.gtf_start}-${gene.gtf_end}`);
  assert.equal(gene.zero_based_half_open[0], gene.gtf_start - 1);
  return {
    chrom: gene.chrom,
    pet: [...pair.a, ...pair.b],
    peak: [peakRow.start, peakRow.end],
    signal: signalRow,
    region: gene.hub_region,
    geneId: gene.gene_id,
  };
}
const egfr = expectedDrawing('Egfr');
const trp53 = expectedDrawing('Trp53');
assert.equal(egfr.chrom, 'chr11');
assert.equal(trp53.chrom, 'chr11');
assert.notDeepEqual(egfr.pet, trp53.pet);

async function drawn(page, expected) {
  const handle = await page.waitForFunction(renderedWindow, expected, { timeout: 30000 });
  return handle.jsonValue();
}
async function chooseGene(page, name, expected) {
  const field = page.getByLabel('Genomic region');
  await field.fill(name);
  const list = page.getByRole('listbox', { name: 'Gene matches' });
  const option = list.locator(
    `[data-symbol="${name}"][data-gene-id="${expected.geneId}"][data-region="${expected.region}"]`,
  );
  await option.waitFor();
  assert.equal(await list.getByRole('option').count(), 1, `${name} is not a single mouse locus`);
  await field.press('Enter');
  await page.waitForFunction(
    value => document.querySelector('[aria-label="Genomic region"]').value === value,
    expected.region,
  );
}
function bytesByFile(reads) {
  const totals = {};
  for (const read of reads) {
    const file = totals[read.name] ??= { document: 0, worker: 0 };
    file[read.scope] += read.bytes;
  }
  return totals;
}

const server = spawn(process.env.HELIXWEAVE_PYTHON ?? 'python3', [
  'test/browser/washu_local_overlay_server.py', '--work', resolve(work),
], { cwd: root, stdio: ['ignore', 'pipe', 'inherit'] });
const serverExit = once(server, 'exit');
const lines = createInterface({ input: server.stdout });
let browser;
let page;
try {
  let port;
  for await (const line of lines) {
    const match = /^WASHU_VIEWER_PORT=(\d+)$/.exec(line);
    if (match) { port = match[1]; break; }
  }
  assert(port, 'the installed HelixWeave server failed to start');
  const origin = `http://127.0.0.1:${port}`;
  browser = await chromium.launch();
  const context = await browser.newContext({ viewport: { width: 1500, height: 1100 } });
  const meter = await installFileReadMeter(context);
  const forbidden = [];
  const routes = new Set();
  await context.addInitScript(() => {
    window.__localViewerCsp = [];
    document.addEventListener('securitypolicyviolation', event => {
      if (event.effectiveDirective === 'connect-src') window.__localViewerCsp.push(event.blockedURI);
    });
  });
  await context.route('**/*', async route => {
    const request = route.request();
    const url = new URL(request.url());
    if (['blob:', 'data:'].includes(url.protocol)) return route.continue();
    const unsafe = url.origin !== origin || !['GET', 'HEAD'].includes(request.method()) ||
      request.postDataBuffer() !== null ||
      /Egfr|Trp53|ENSMUSG|MGI:95|pet\.longrange|signal\.bw|peak\.bed/i.test(url.href);
    if (unsafe) {
      forbidden.push([request.method(), url.href]);
      await route.abort();
    } else {
      routes.add(url.pathname + url.search);
      await route.continue();
    }
  });
  page = await context.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  page.on('console', message => {
    if (/Minified React error|Maximum update depth exceeded|Genome View Error/.test(message.text())) {
      errors.push(message.text());
    }
  });
  await page.goto(origin);
  await page.getByRole('link', { name: '打开本地轨道', exact: true }).click();
  await page.waitForURL(url => url.pathname === '/viewer.html' && !url.search);
  const field = page.getByLabel('Genomic region');
  await page.getByLabel('Assembly', { exact: true }).selectOption('mm10');
  await page.getByText('当前声明 GRCm38', { exact: true }).waitFor();
  const mark = meter.mark();
  const started = performance.now();
  await page.getByLabel('Local track files').setInputFiles([data, index, signal, peak]);
  const rows = page.getByTestId('local-track');
  await rows.nth(2).waitFor();
  const described = await rows.evaluateAll(nodes => nodes.map(node => ({
    assembly: node.dataset.assembly,
    type: node.dataset.trackType,
    read: node.dataset.readMethod,
    text: node.textContent,
  })));
  assert.deepEqual(described.map(row => row.assembly), ['mm10', 'mm10', 'mm10']);
  assert.deepEqual(described.map(row => row.type).sort(), ['bed', 'bigwig', 'longrange']);
  assert.deepEqual(described.map(row => row.read).sort(), [
    'wuepgg3-track-bigwig-file-slice',
    'wuepgg3-track-tabix-file-slice',
    'wuepgg3-track-text-file',
  ]);
  for (const row of described) {
    assert.match(row.text, /用户指定 (?:mm10|GRCm38)/);
    assert.equal(/GRCh38|\bhg38\b/.test(row.text), false);
  }
  assert.equal(await page.getByTestId('reference-conflict').count(), 0);
  assert.equal(await page.locator('[data-reason="assembly-conflict"]').count(), 0);

  await chooseGene(page, 'Egfr', egfr);
  const egfrDrawing = await drawn(page, egfr);
  const elapsedMs = Math.round(performance.now() - started);
  const firstWindow = {
    elapsedMs,
    bytesByFile: bytesByFile(meter.since(mark)),
    dataFileBytes: (await stat(data)).size,
    indexFileBytes: (await stat(index)).size,
  };
  const petReads = firstWindow.bytesByFile[basename(data)];
  assert.ok(petReads?.worker > 0 && firstWindow.bytesByFile[basename(index)]?.worker > 0);
  if (firstWindow.dataFileBytes > 32 * 1024 ** 2) {
    assert.ok(observation.record_count > 2_000_000);
    assert.ok(petReads.document + petReads.worker < firstWindow.dataFileBytes / 10);
  }
  assert.equal(await field.inputValue(), egfr.region);
  assert.equal(egfrDrawing.pet.a.start, egfr.pet[0]);
  assert.equal(egfrDrawing.peak.start, egfr.peak[0]);
  await page.getByTestId('washu-genome-viewer').screenshot({ path: join(work, 'Egfr.png') });

  await chooseGene(page, 'Trp53', trp53);
  const trp53Drawing = await drawn(page, trp53);
  assert.equal(await field.inputValue(), trp53.region);
  assert.equal(trp53Drawing.signal.value, trp53.signal.value);
  await page.getByTestId('washu-genome-viewer').screenshot({ path: join(work, 'Trp53.png') });

  const probe = planted.coordinate_probe.region;
  await field.fill(probe);
  assert.equal(await field.inputValue(), probe);
  const coordinateDrawing = await drawn(page, { ...egfr, region: probe });
  assert.equal(await page.evaluate(renderedWindow, egfr), null,
    'the ruler still matches the wider Egfr gene window');
  await page.getByTestId('washu-genome-viewer').screenshot({ path: join(work, 'coordinate.png') });

  const byteCounts = bytesByFile(meter.since(mark));
  const dataReads = byteCounts[basename(data)];
  assert.ok(dataReads?.worker > 0, 'longrange bytes were not read in the worker');
  assert.ok(byteCounts[basename(index)]?.worker > 0, 'TBI bytes were not read in the worker');
  assert.deepEqual(await page.evaluate(() => window.__localViewerCsp), []);
  assert.deepEqual(forbidden, []);
  assert.deepEqual([...routes].filter(path => path.includes('gene-index')).map(path =>
    path.includes('mm10-gene-index') ? 'mm10' : path), ['mm10']);
  assert.deepEqual(errors, []);
  const result = {
    source: 'synthetic',
    annotation: 'GENCODE M25',
    assembly: 'mm10',
    coordinateConvention: '1-based closed for search; BED half-open',
    records: observation.record_count,
    firstWindow,
    egfr: egfr.region,
    trp53: trp53.region,
    coordinate: probe,
    egfrDrawing,
    trp53Drawing,
    coordinateDrawing,
    byteCounts,
    tested: ['Egfr gene', 'Trp53 gene', '1-based coordinate', 'three local tracks', 'no upload'],
    networkPaths: [...routes],
    forbidden,
    errors,
  };
  await writeFile(join(work, 'browser.json'), JSON.stringify(result, null, 2) + '\n');
  console.log('WASHU_MOUSE_LOCAL_FILES ok', JSON.stringify({
    work,
    records: result.records,
    egfr: egfr.region,
    trp53: trp53.region,
    coordinate: probe,
    byteCounts,
    egfrPet: egfrDrawing.pet,
    trp53Signal: trp53Drawing.signal,
  }));
} catch (error) {
  if (page) {
    await page.screenshot({ path: join(work, 'failure.png') }).catch(() => {});
    await writeFile(join(work, 'failure.txt'), await page.locator('body').innerText().catch(() => ''));
  }
  console.error('WashU mouse browser evidence', work);
  throw error;
} finally {
  if (browser) await browser.close();
  lines.close();
  if (server.exitCode === null && server.signalCode === null) server.kill('SIGTERM');
  await serverExit;
}
