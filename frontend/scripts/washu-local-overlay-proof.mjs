// Start at the home page of an installed package. All chosen files are synthetic.
import assert from 'node:assert/strict';
import { spawn, spawnSync } from 'node:child_process';
import { once } from 'node:events';
import { copyFile, mkdir, mkdtemp, readFile, stat, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { basename, join, resolve } from 'node:path';
import { createInterface } from 'node:readline';
import { fileURLToPath } from 'node:url';
import { chromium } from '@playwright/test';
import { installFileReadMeter } from './washu-file-reads.mjs';
import { renderedChromosome, renderedWindow } from './washu-local-overlay-observe.mjs';

const root = fileURLToPath(new URL('../../', import.meta.url));
const work = process.env.WASHU_LOCAL_OVERLAY_DIR ?? await mkdtemp(join(tmpdir(), 'washu-local-'));
await mkdir(work, { recursive: true });
const fixtures = resolve(process.env.WASHU_LOCAL_FIXTURES ?? join(work, 'fixtures'));
if (!process.env.WASHU_LOCAL_FIXTURES) {
  const generate = spawnSync(process.env.WASHU_FIXTURE_PYTHON ?? 'python3', [
    'test/browser/local_overlay_fixtures.py', '--dest', fixtures,
  ], { cwd: root, stdio: 'inherit' });
  if (generate.error || generate.status !== 0) throw generate.error ?? new Error('Fixture generation failed');
}
const data = join(fixtures, 'pet.longrange.gz');
const index = data + '.tbi';
const signal = join(fixtures, 'signal.bw');
const peak = join(fixtures, 'peak.bed');
const observation = JSON.parse(await readFile(join(fixtures, 'observation.json'), 'utf8'));
const dataSize = (await stat(data)).size;
assert(dataSize > 32 * 1024 ** 2 && observation.record_count > 2_000_000);
const egfr = {
  chrom: 'chr7', pet: [55101000, 55102000, 55106000, 55107000],
  signal: { start: 55104400, end: 55105900, value: 4.5 }, peak: [55103200, 55104100],
};

async function drawn(page, expected) {
  const handle = await page.waitForFunction(renderedWindow, expected, { timeout: 20000 });
  return handle.jsonValue();
}
async function chooseGene(page, name, symbol, region) {
  const field = page.getByLabel('Genomic region');
  await field.fill(name);
  await page.getByRole('listbox', { name: 'Gene matches' }).locator(`[data-symbol="${symbol}"]`).first().waitFor();
  await field.press('Enter');
  await page.waitForFunction(value => document.querySelector('[aria-label="Genomic region"]').value === value, region);
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
      request.postDataBuffer() !== null || /EGFR|TP53|SHOX|ERBB1|PHOG|pet\.longrange|signal\.bw|peak\.bed/i.test(url.href);
    if (unsafe) {
      forbidden.push([request.method(), url.href]);
      await route.abort();
    } else {
      routes.add(url.pathname);
      await route.continue();
    }
  });
  page = await context.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(origin);
  await page.getByRole('link', { name: '打开本地轨道', exact: true }).click();
  await page.waitForURL(url => url.pathname === '/viewer.html' && !url.search);
  const field = page.getByLabel('Genomic region');
  const fileInput = page.getByLabel('Local track files');
  await page.getByLabel('Assembly', { exact: true }).selectOption('hg38');
  const mark = meter.mark();
  const start = performance.now();
  await fileInput.setInputFiles([data, index, signal, peak]);
  await chooseGene(page, 'EGFR', 'EGFR', 'chr7:55018820-55211628');
  const visible = await drawn(page, egfr);
  const elapsedMs = Math.round(performance.now() - start);
  const byteCounts = bytesByFile(meter.since(mark));
  const dataReads = byteCounts[basename(data)];
  assert(dataReads?.worker > 0 && dataReads.worker + dataReads.document < dataSize / 10);
  assert(byteCounts[basename(index)]?.worker > 0);
  await page.getByTestId('washu-genome-viewer').screenshot({ path: join(work, 'EGFR.png') });

  await chooseGene(page, 'TP53', 'TP53', 'chr17:7661779-7687546');
  await page.waitForFunction(renderedChromosome, 'chr17');
  await field.fill('chr7:55100000-55110000');
  await drawn(page, egfr);
  await page.getByTestId('washu-genome-viewer').screenshot({ path: join(work, 'coordinate.png') });
  await field.fill('SHOX');
  const matches = page.getByRole('listbox', { name: 'Gene matches' });
  await matches.locator('[data-region^="chrX:"]').waitFor();
  await matches.locator('[data-region^="chrY:"]').waitFor();
  await field.press('Enter');
  assert.equal(await field.inputValue(), 'SHOX', 'ambiguous input jumped without the user choosing');
  await matches.locator('[data-region^="chrY:"]').click();
  await page.waitForFunction(renderedChromosome, 'chrY');
  await field.fill('HELIXWEAVE_NONEXISTENT_GENE');
  await page.getByText(/no genes match|未找到|没有.*匹配/i).first().waitFor();
  await chooseGene(page, 'ERBB1', 'EGFR', 'chr7:55018820-55211628');
  await drawn(page, egfr);

  await fileInput.setInputFiles(data);
  await page.locator('[data-reason="missing-index"]').waitFor();
  assert.equal(await page.getByTestId('local-track').count(), 3);
  const wrongDir = join(work, 'wrong');
  await mkdir(wrongDir, { recursive: true });
  const wrong = join(wrongDir, basename(index));
  await copyFile(join(fixtures, 'wrong-index.tbi'), wrong);
  await fileInput.setInputFiles([data, wrong]);
  await page.locator('[role="alert"][data-reason*="index"]:not([data-reason="missing-index"])').waitFor();
  await drawn(page, egfr);
  await fileInput.setInputFiles({
    name: 'explicit-reference.bed', mimeType: 'text/plain',
    buffer: Buffer.from('track db=hg19\nchr7\t4\t10\n'),
  });
  await page.locator('[data-reason="assembly-conflict"]').waitFor();
  await page.getByLabel('Assembly', { exact: true }).selectOption('mm10');
  await page.getByTestId('reference-conflict').first().waitFor();
  await page.getByLabel('Assembly', { exact: true }).selectOption('hg38');
  await page.getByRole('button', { name: 'Remove local track pet.longrange.gz' }).click();
  assert.equal(await page.getByTestId('local-track').count(), 2);
  const renamed = join(work, 'index-original.tbi');
  await copyFile(index, renamed);
  await fileInput.setInputFiles([data, renamed]);
  const choice = page.getByTestId('local-track-choice');
  await choice.getByRole('button', { name: /index-original\.tbi/ }).click();
  await drawn(page, egfr);
  assert.equal(await page.getByTestId('local-track').count(), 3);
  assert.deepEqual(await page.evaluate(() => window.__localViewerCsp), []);
  assert.deepEqual(forbidden, []);
  assert.deepEqual(errors, []);
  const result = {
    source: 'synthetic', annotation: 'GENCODE v49', coordinateConvention: '1-based closed for search; BED half-open',
    records: observation.record_count, dataFileBytes: dataSize, indexFileBytes: (await stat(index)).size,
    elapsedMs, byteCounts, visible, secondGene: 'chr17:7661779-7687546', coordinate: 'chr7:55100000-55110000',
    tested: ['SHOX chrX/chrY', 'ERBB1 alias', 'missing gene', 'missing TBI', 'wrong valid TBI',
      'renamed TBI', 'explicit BED assembly conflict', 'switch assembly', 'remove and reopen'],
    networkPaths: [...routes], forbidden, errors,
  };
  await writeFile(join(work, 'browser.json'), JSON.stringify(result, null, 2) + '\n');
  console.log('WASHU_LOCAL_FILES ok', JSON.stringify({ work, records: result.records, elapsedMs, byteCounts, visible }));
} catch (error) {
  if (page) {
    await page.screenshot({ path: join(work, 'failure.png') }).catch(() => {});
    await writeFile(join(work, 'failure.txt'), await page.locator('body').innerText().catch(() => ''));
  }
  console.error('WashU browser evidence', work);
  throw error;
} finally {
  if (browser) await browser.close();
  lines.close();
  if (server.exitCode === null && server.signalCode === null) server.kill('SIGTERM');
  await serverExit;
}
