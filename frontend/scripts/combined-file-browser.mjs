import { createHash } from 'node:crypto';
import { spawn } from 'node:child_process';
import { existsSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import net from 'node:net';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from '@playwright/test';

const repoRoot = fileURLToPath(new URL('../..', import.meta.url));
const serverPath = fileURLToPath(
  new URL('../../test/browser/combined_file_server.py', import.meta.url),
);
const python =
  process.env.HELIXWEAVE_PYTHON ?? join(repoRoot, '.local/envs/ci-fast/bin/python');
const browsers =
  process.env.PLAYWRIGHT_BROWSERS_PATH ??
  '/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/ms-playwright';
if (existsSync(browsers)) process.env.PLAYWRIGHT_BROWSERS_PATH = browsers;

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

const work = mkdtempSync(join(tmpdir(), 'helixweave-combined-file-'));
const children = [];

function start(command, args) {
  const child = spawn(command, args, {
    cwd: repoRoot,
    env: {
      ...process.env,
      PYTHONPATH: fileURLToPath(new URL('../../src', import.meta.url)),
    },
    stdio: ['ignore', 'pipe', 'pipe'],
    detached: true,
  });
  children.push(child);
  return child;
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
      finish(new Error(`combined server exited before it was ready (${code})\n${buffer}`));
    });
    child.stderr.on('data', (chunk) => {
      process.stderr.write(chunk);
    });
    child.stdout.on('data', (chunk) => {
      buffer += chunk.toString();
      const line = buffer
        .split('\n')
        .find((item) => item.startsWith('COMBINED_FILE_READY '));
      if (!line) return;
      try {
        finish(null, JSON.parse(line.slice('COMBINED_FILE_READY '.length)));
      } catch (error) {
        finish(error instanceof Error ? error : new Error(String(error)));
      }
    });
  });
}

function sha256(bytes) {
  return createHash('sha256').update(bytes).digest('hex');
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

function step(message) {
  console.error(`STEP ${message}`);
}

async function stopAll(browser) {
  await Promise.race([
    browser ? browser.close().catch(() => {}) : Promise.resolve(),
    new Promise((resolve) => setTimeout(resolve, 3_000)),
  ]);
  for (const child of children) {
    try {
      process.kill(-child.pid, 'SIGTERM');
    } catch {
      child.kill('SIGTERM');
    }
  }
}

async function main() {
  const hardStop = setTimeout(() => {
    console.error('combined browser exceeded 90s');
    for (const child of children) {
      try {
        process.kill(-child.pid, 'SIGKILL');
      } catch {
        child.kill('SIGKILL');
      }
    }
    process.exit(1);
  }, 90_000);
  hardStop.unref?.();
  const port = await reserve(4367);
  const server = start(python, [
    serverPath,
    '--port',
    String(port),
    '--database',
    join(work, 'combined.sqlite'),
    '--workspace',
    join(work, 'workspace'),
    '--durable',
    join(work, 'durable'),
    '--request',
    join(work, 'relocate.request'),
    '--result',
    join(work, 'relocate.result.json'),
  ]);
  const manifest = await waitForReady(server);
  step('launch');
  const browser = await chromium.launch({
    headless: true,
    args: ['--no-proxy-server'],
  });
  const page = await browser.newPage({ baseURL: manifest.api_base });
  page.setDefaultTimeout(20_000);
  page.on('pageerror', (error) => {
    console.error(`PAGE ${error.message}`);
  });
  const responses = [];
  page.on('response', (response) => {
    const url = response.url();
    if (!url.includes('/api/')) return;
    responses.push(`${response.status()} ${url}`);
  });
  try {
    step('login');
    await page.goto('/');
    await page.locator('#login-username').fill(manifest.admin_username);
    await page.locator('#login-password').fill(manifest.admin_password);
    await page.getByRole('button', { name: 'Sign in' }).click();
    await page.waitForURL(/\/workflows$/);

    step('catalog');
    await page.goto('/catalog');
    await page.getByRole('heading', { name: 'Experiment catalog' }).waitFor();
    const search = page.getByRole('form', { name: 'Find experiments' });
    await search.getByLabel('Project').selectOption(manifest.project_id);
    await search.getByLabel('Assay').selectOption([manifest.assay_term_id]);
    await search.getByLabel('Cell type').selectOption([manifest.cell_type_term_id]);
    await search.getByRole('button', { name: 'Find experiments' }).click();
    await page.getByRole('button', { name: manifest.experiment_name }).click();
    await page.getByRole('heading', { name: manifest.experiment_name }).waitFor();
    const scope = page.getByRole('region', { name: 'Run input scope' });
    await scope.getByText(manifest.artifact_id).waitFor();
    const scopeText = await scope.innerText();
    if (scopeText.includes('这条样本修订没有运行结果')) {
      throw new Error(scopeText);
    }

    step(`scope ${scopeText.replace(/\s+/g, ' ').slice(0, 240)}`);
    step('viewer');
    const sessionResponse = page.waitForResponse(
      (response) =>
        response.url().includes('/browser-session') &&
        response.request().method() === 'GET' &&
        response.status() === 200,
    );
    await page.goto(`/viewer.html?run=${manifest.run_id}`);
    const session = await (await sessionResponse).json();
    const track = (session.tracks ?? []).find(
      (item) => item.artifact_id === manifest.artifact_id,
    );
    if (!track) {
      throw new Error(`track missing from session ${JSON.stringify(session)}`);
    }
    if (
      session.run_id !== manifest.run_id ||
      session.artifact_generation !== manifest.generation ||
      track.revision !== manifest.revision ||
      track.source_form !== 'display_derivative' ||
      track.preparation !== 'needed' ||
      track.provenance != null
    ) {
      throw new Error(
        `track was already converted before it was chosen ${JSON.stringify(track)}`,
      );
    }
    step('checkbox');
    const clicked = await Promise.race([
      page.evaluate(async (name) => {
        const label = [...document.querySelectorAll('label')].find((item) =>
          item.textContent.includes(name),
        );
        const input = label?.querySelector('input[type="checkbox"]');
        if (!input) return 'checkbox missing';
        input.click();
        const deadline = Date.now() + 5_000;
        let button = null;
        while (Date.now() < deadline) {
          button = [...document.querySelectorAll('button')].find((item) =>
            item.textContent.includes(name),
          );
          if (button) break;
          await new Promise((resolve) => setTimeout(resolve, 50));
        }
        if (!button) return 'prepare button missing';
        button.click();
        return 'clicked';
      }, manifest.peak_name),
      new Promise((resolve) => setTimeout(() => resolve('evaluate timed out'), 15_000)),
    ]);
    step(clicked);
    if (clicked !== 'clicked') {
      throw new Error(clicked);
    }
    const trackDeadline = Date.now() + 20_000;
    let trackRead = null;
    while (Date.now() < trackDeadline) {
      trackRead = responses.find(
        (line) =>
          line.includes(
            `/browser-tracks/${manifest.artifact_id}/bytes/${manifest.generation}/${manifest.revision}`,
          ) && !line.includes('.tbi'),
      );
      if (trackRead) break;
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    if (!trackRead) {
      throw new Error(`track bytes were not read ${JSON.stringify(responses)}`);
    }
    step(trackRead);
    const trackUrl = trackRead.slice(trackRead.indexOf('http'));
    const trackResponse = await page.request.get(trackUrl);
    if (trackResponse.status() !== 200 && trackResponse.status() !== 206) {
      throw new Error(`track read returned ${trackResponse.status()}`);
    }
    const trackBody = await trackResponse.body();
    const prepared = await page.request.get(
      `/api/v1/runs/${manifest.run_id}/browser-session`,
    );
    if (prepared.status() !== 200) {
      throw new Error(`browser session returned ${prepared.status()}`);
    }
    const preparedSession = await prepared.json();
    const preparedTrack = (preparedSession.tracks ?? []).find(
      (item) => item.artifact_id === manifest.artifact_id,
    );
    if (
      preparedTrack?.provenance?.source_artifact_id !== manifest.artifact_id ||
      preparedTrack?.provenance?.source_revision !== manifest.revision ||
      preparedTrack?.provenance?.source_generation !== manifest.generation
    ) {
      throw new Error(
        `prepared source identity drifted ${JSON.stringify(preparedTrack)}`,
      );
    }
    track.source_form = preparedTrack.source_form;
    step('download-before');
    const before = await page.request.get(manifest.download_path);
    if (before.status() !== 200) {
      throw new Error(`original download returned ${before.status()} ${await before.text()}`);
    }
    const beforeBody = await before.body();
    const beforeSha = sha256(beforeBody);
    if (beforeSha !== manifest.source_sha256 || beforeBody.length !== manifest.source_size) {
      throw new Error(
        `original download sha ${beforeSha} size ${beforeBody.length} does not match the source`,
      );
    }

    step('relocate');
    writeFileSync(manifest.relocate_request, 'relocate\n');
    const deadline = Date.now() + 30_000;
    let relocated = null;
    while (Date.now() < deadline) {
      if (existsSync(manifest.relocate_result)) {
        relocated = JSON.parse(readFileSync(manifest.relocate_result, 'utf8'));
        break;
      }
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    if (!relocated || relocated.ok !== 'true') {
      throw new Error(`relocation did not finish ${JSON.stringify(relocated)}`);
    }
    if (
      relocated.run_id !== manifest.run_id ||
      relocated.artifact_id !== manifest.artifact_id ||
      relocated.generation !== manifest.generation ||
      relocated.revision !== manifest.revision ||
      relocated.content_sha256 !== beforeSha
    ) {
      throw new Error(`relocation identity drifted ${JSON.stringify(relocated)}`);
    }
    step('download-after');
    const after = await page.request.get(manifest.download_path);
    if (after.status() !== 200) {
      throw new Error(`durable download returned ${after.status()} ${await after.text()}`);
    }
    const afterBody = await after.body();
    const afterSha = sha256(afterBody);
    if (afterSha !== beforeSha || !afterBody.equals(beforeBody)) {
      throw new Error(`durable download sha ${afterSha} does not match the saved body`);
    }
    const report = {
      experiment_revision_id: manifest.experiment_revision_id,
      sample_revision_id: manifest.sample_revision_id,
      run_id: manifest.run_id,
      artifact_id: manifest.artifact_id,
      generation: manifest.generation,
      revision: manifest.revision,
      source_sha256: manifest.source_sha256,
      download_sha256_before: beforeSha,
      download_sha256_after: afterSha,
      relocation_sha256: relocated.content_sha256,
      track_read: trackRead,
      track_sha256: sha256(trackBody),
      track_bytes_match_download: sha256(trackBody) === beforeSha,
      track_source_form: track.source_form,
      same_file: true,
    };
    clearTimeout(hardStop);
    console.log(JSON.stringify(report, null, 2));
  } catch (error) {
    console.error(error instanceof Error ? error.stack ?? error.message : String(error));
    throw error;
  } finally {
    await stopAll(browser);
  }
}

main().catch(async (error) => {
  console.error(error instanceof Error ? error.stack ?? error.message : String(error));
  console.error(`combined work ${work}`);
  await stopAll();
  process.exitCode = 1;
});
