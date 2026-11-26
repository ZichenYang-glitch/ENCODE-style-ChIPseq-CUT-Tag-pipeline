/**
 * Packaged catalog acceptance for a Peri target suggestion.
 *
 * Starts test/browser/catalog_classification_server.py, signs in through the
 * packaged page, and checks the synthetic experiments. When the page exposes
 * Ask Peri / 询问 Peri, Preview selected / 预览所选建议, and
 * Adopt previewed changes / 采用预览中的变更, it generates a suggestion,
 * previews one item, and adopts only that preview. Generation must leave
 * both revisions unchanged. Adoption may advance exactly one. A later quota
 * refusal must not add a loopback call.
 *
 * Offline proposal registration is not treated as a model source.
 * Run: node frontend/scripts/catalog-classification-browser.mjs
 * Python: HELIXWEAVE_PYTHON. Screenshot: CLASSIFICATION_SCREENSHOT.
 */
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { chromium } from '@playwright/test';

const repoRoot = fileURLToPath(new URL('../..', import.meta.url));
const serverPath = fileURLToPath(
  new URL('../../test/browser/catalog_classification_server.py', import.meta.url),
);
const python = process.env.HELIXWEAVE_PYTHON ?? 'python3';
const screenshotPath = process.env.CLASSIFICATION_SCREENSHOT ?? '';

const ASK = /^(Ask Peri|询问 Peri)$/;
const PREVIEW = /^(Preview selected|预览所选建议)$/;
const ADOPT = /^(Adopt previewed changes|采用预览中的变更)$/;
const QUOTA_CODES = new Set([
  'SCOPE_DENIED',
  'QUOTA_EXCEEDED',
  'QUOTA_EXHAUSTED',
  'BUDGET_EXCEEDED',
  'GRANT_EXHAUSTED',
  'MODEL_BUDGET_EXCEEDED',
  'INSUFFICIENT_QUOTA',
]);
const KNOWN_CATALOG_POSTS = new Set([
  'catalogRegisterProposal',
  'catalogPreviewProposal',
  'catalogAdoptProposal',
  'catalogCreateExperiment',
  'catalogCreateTerm',
  'catalogReviseExperiment',
  'catalogReviseTerm',
  'catalogReviseSampleCatalog',
  'catalogArchiveTerm',
  'catalogArchiveExperiment',
  'catalogRegisterRegisteredResult',
  'catalogReviseRegisteredResult',
  'catalogArchiveRegisteredResult',
]);

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
      finish(
        new Error(`classification fixture exited before it was ready (${code})\n${buffer}`),
      );
    });
    child.stderr.on('data', (chunk) => {
      process.stderr.write(chunk);
    });
    child.stdout.on('data', (chunk) => {
      buffer += chunk.toString();
      const line = buffer
        .split('\n')
        .find((item) => item.startsWith('CLASSIFICATION_FIXTURE_READY '));
      if (!line) return;
      try {
        finish(null, JSON.parse(line.slice('CLASSIFICATION_FIXTURE_READY '.length)));
      } catch (error) {
        finish(error instanceof Error ? error : new Error(String(error)));
      }
    });
  });
}

function denialCode(body) {
  if (!body || typeof body !== 'object') return '';
  const nested = body.error && typeof body.error === 'object' ? body.error : {};
  return String(
    body.failure_reason || body.reason_code || body.code || body.error_code || nested.code || '',
  );
}

function isQuotaDenial(status, body) {
  if (status === 402 || status === 429) return true;
  return QUOTA_CODES.has(denialCode(body));
}

async function providerCalls(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`provider counter returned ${response.status}`);
  const body = await response.json();
  if (typeof body.count !== 'number') throw new Error('provider counter has no count');
  return body.count;
}

async function revisionNumbers(page, manifest) {
  const response = await page.request.get(
    `/api/v1/catalog/projects/${manifest.project_id}/experiments?limit=20`,
  );
  if (!response.ok()) throw new Error(`experiment list returned ${response.status()}`);
  const body = await response.json();
  const byId = new Map(
    (body.experiments ?? []).map((item) => [item.experiment_id, item.revision_number]),
  );
  return {
    selected: byId.get(manifest.selected.experiment_id),
    excluded: byId.get(manifest.excluded.experiment_id),
  };
}

function sameRevisions(left, right) {
  return left.selected === right.selected && left.excluded === right.excluded;
}

function generateOperations(spec) {
  const found = [];
  for (const [path, methods] of Object.entries(spec.paths ?? {})) {
    for (const [method, operation] of Object.entries(methods ?? {})) {
      if (!operation || typeof operation !== 'object' || !operation.operationId) continue;
      if (method !== 'post' || !path.includes('/catalog')) continue;
      if (KNOWN_CATALOG_POSTS.has(operation.operationId)) continue;
      const blob = `${operation.operationId} ${operation.summary ?? ''} ${path}`.toLowerCase();
      if (/suggest|generat|classif|peri/.test(blob)) {
        found.push({ operationId: operation.operationId, path });
      }
    }
  }
  return found;
}

async function buttonCount(page, pattern) {
  return page.getByRole('button', { name: pattern }).count();
}

const api = spawn(python, [serverPath], {
  cwd: repoRoot,
  env: {
    ...process.env,
    PYTHONPATH: fileURLToPath(new URL('../../src', import.meta.url)),
  },
  stdio: ['ignore', 'pipe', 'pipe'],
});

const report = {
  tested: {},
  pending: [],
  page_flow: 'pending',
};

try {
  const manifest = await waitForReady(api);
  if (!manifest.database_is_temporary || manifest.offline_proposal_registered) {
    throw new Error('fixture must use a temporary database and must not register offline');
  }
  const browser = await chromium.launch({
    headless: true,
    args: ['--no-proxy-server'],
  });
  const page = await browser.newPage({ baseURL: manifest.api_base });
  page.setDefaultTimeout(20_000);
  try {
    const catalogDocument = page.waitForResponse(
      (response) =>
        response.request().method() === 'GET' &&
        new URL(response.url()).pathname === '/catalog',
    );
    await page.goto('/');
    await page.locator('#login-username').fill(manifest.admin_username);
    await page.locator('#login-password').fill(manifest.admin_password);
    await page.getByRole('button', { name: 'Sign in' }).click();
    await page.waitForURL(/\/workflows$/);
    await page.goto('/catalog');
    const catalogResponse = await catalogDocument;
    if (catalogResponse.status() !== 200) {
      throw new Error(`packaged /catalog returned ${catalogResponse.status()}`);
    }
    const contentType = catalogResponse.headers()['content-type'] ?? '';
    if (!contentType.includes('text/html')) {
      throw new Error(`packaged /catalog content type was ${contentType}`);
    }
    report.tested.packaged_catalog_document = true;
    report.tested.synthetic_admin_login = true;

    const search = page.getByRole('form', { name: 'Find experiments' });
    await search.getByLabel('Project').selectOption(manifest.project_id);
    await search.getByRole('button', { name: 'Find experiments' }).click();
    const experiments = page.getByRole('list', { name: 'Experiments' });
    for (const item of [manifest.selected, manifest.excluded]) {
      const card = experiments.getByRole('button', { name: item.display_name });
      await card.waitFor();
      await card.getByText('修订 1').waitFor();
    }
    await experiments.getByRole('button', { name: manifest.selected.display_name }).click();
    const detail = page.getByRole('region', { name: 'Experiment revision' });
    await detail.getByRole('heading', { name: manifest.selected.display_name }).waitFor();
    await detail.locator('p', { hasText: '当前 1' }).waitFor();
    await detail.getByText('未知').first().waitFor();
    const opened = await revisionNumbers(page, manifest);
    if (opened.selected !== 1 || opened.excluded !== 1) {
      throw new Error(`seed revisions were ${JSON.stringify(opened)}`);
    }
    report.tested.seeded_revisions_visible = true;

    const idleCalls = await providerCalls(manifest.provider_calls);
    if (idleCalls !== 0) {
      throw new Error(`opening the catalog called the provider ${idleCalls} times`);
    }
    report.tested.provider_idle_on_page_load = true;
    report.tested.database_is_temporary = true;

    const openapiResponse = await page.request.get('/openapi.json');
    if (!openapiResponse.ok()) throw new Error(`openapi returned ${openapiResponse.status()}`);
    const openapi = await openapiResponse.json();
    const operationIds = new Set();
    for (const methods of Object.values(openapi.paths ?? {})) {
      for (const operation of Object.values(methods ?? {})) {
        if (operation && operation.operationId) operationIds.add(operation.operationId);
      }
    }
    report.m2 = {
      preview: operationIds.has('catalogPreviewProposal'),
      adopt: operationIds.has('catalogAdoptProposal'),
      offline_register: operationIds.has('catalogRegisterProposal'),
    };
    if (!report.m2.preview || !report.m2.adopt) {
      throw new Error('packaged origin is missing the M2 preview or adopt operation');
    }
    report.tested.m2_preview_and_adopt_operations = true;
    report.generate_operations = generateOperations(openapi);

    let askCount = await buttonCount(page, ASK);
    let previewCount = await buttonCount(page, PREVIEW);
    const language = page.getByLabel(/^(Language|语言)$/);
    report.language_select = (await language.count()) > 0;
    if (askCount > 0 && report.language_select) {
      await language.selectOption({ label: '简体中文' });
      const chineseAsk = await page.getByRole('button', { name: '询问 Peri' }).count();
      const chinesePreview = await page
        .getByRole('button', { name: '预览所选建议' })
        .count();
      if (chineseAsk !== 1 || chinesePreview < 1) {
        throw new Error('language select did not expose the Chinese classification controls');
      }
      report.tested.chinese_controls = true;
      await language.selectOption({ label: 'English' });
      askCount = await buttonCount(page, ASK);
      previewCount = await buttonCount(page, PREVIEW);
      if (askCount === 0) {
        throw new Error('English Ask Peri disappeared after the language switch');
      }
    }
    report.controls = { ask: askCount, preview: previewCount };

    if (askCount === 0) {
      report.page_flow = 'pending';
      report.pending.push(
        'ask_peri_generate',
        'preview_selected_only',
        'adopt_previewed_changes',
        'revision_unchanged_until_adopt',
        'quota_refusal_zero_extra_calls',
      );
      if (!report.language_select) report.pending.push('chinese_controls');
      if (report.generate_operations.length === 0) {
        report.pending.push('public_generate_operation');
      }
    } else {
      const before = await revisionNumbers(page, manifest);
      const beforeCalls = await providerCalls(manifest.provider_calls);
      const generated = page.waitForResponse(
        (response) =>
          response.request().method() === 'POST' &&
          response.url().startsWith(manifest.api_base),
      );
      await page.getByRole('button', { name: ASK }).click();
      const generateResponse = await generated;
      const generateBody = await generateResponse.json().catch(() => ({}));
      const generateRequest = generateResponse.request().postDataJSON();
      if (generateRequest && generateRequest.source === 'offline_suggestion') {
        throw new Error('Ask Peri submitted an offline suggestion');
      }
      const callsAfterGenerate = await providerCalls(manifest.provider_calls);
      if (isQuotaDenial(generateResponse.status(), generateBody)) {
        if (callsAfterGenerate !== beforeCalls) {
          throw new Error('quota refusal reached the loopback provider');
        }
        report.tested.quota_refusal_zero_extra_calls = true;
        report.pending.push(
          'ask_peri_generate',
          'preview_selected_only',
          'adopt_previewed_changes',
          'revision_unchanged_until_adopt',
        );
      } else {
        if (generateResponse.status() >= 400) {
          throw new Error(
            `Ask Peri failed with ${generateResponse.status()} ${JSON.stringify(generateBody)}`,
          );
        }
        if (callsAfterGenerate !== beforeCalls + 1) {
          throw new Error(
            `Ask Peri made ${callsAfterGenerate - beforeCalls} provider calls`,
          );
        }
        const afterGenerate = await revisionNumbers(page, manifest);
        if (!sameRevisions(before, afterGenerate)) {
          throw new Error('generation changed a catalog revision');
        }
        report.tested.generate_leaves_revisions_unchanged = true;
        const boxes = page
          .getByRole('group', { name: /Choose target changes to preview|选择要预览的靶标变更/ })
          .getByRole('checkbox');
        const choiceCount = await boxes.count();
        if (choiceCount < 2) {
          throw new Error(`expected two target choices, saw ${choiceCount}`);
        }
        await boxes.nth(0).check();
        for (let index = 1; index < choiceCount; index += 1) {
          await boxes.nth(index).uncheck();
        }
        const previewed = page.waitForResponse(
          (response) =>
            response.request().method() === 'POST' &&
            new URL(response.url()).pathname.endsWith('/preview'),
        );
        await page.getByRole('button', { name: PREVIEW }).click();
        const previewResponse = await previewed;
        if (previewResponse.status() !== 200) {
          throw new Error(`preview returned ${previewResponse.status()}`);
        }
        const previewBody = await previewResponse.json();
        if (
          previewBody.selected_item_ids?.length !== 1 ||
          previewBody.excluded_item_ids?.length < 1 ||
          previewBody.has_conflicts
        ) {
          throw new Error(`preview was not a single conflict-free subset: ${JSON.stringify(previewBody)}`);
        }
        const afterPreview = await revisionNumbers(page, manifest);
        if (!sameRevisions(before, afterPreview)) {
          throw new Error('preview changed a catalog revision');
        }
        report.tested.preview_selected_only = true;
        if ((await buttonCount(page, ADOPT)) !== 1) {
          throw new Error('preview did not reveal Adopt previewed changes');
        }
        const adopted = page.waitForResponse(
          (response) =>
            response.request().method() === 'POST' &&
            new URL(response.url()).pathname.endsWith('/adoptions'),
        );
        await page.getByRole('button', { name: ADOPT }).click();
        const adoptResponse = await adopted;
        if (adoptResponse.status() !== 200) {
          throw new Error(`adopt returned ${adoptResponse.status()}`);
        }
        const adoptBody = await adoptResponse.json();
        if (adoptBody.selected_item_ids?.length !== 1) {
          throw new Error('adoption did not keep the previewed subset');
        }
        const afterAdopt = await revisionNumbers(page, manifest);
        const advanced = [afterAdopt.selected, afterAdopt.excluded].filter(
          (revision) => revision === 2,
        );
        const held = [afterAdopt.selected, afterAdopt.excluded].filter(
          (revision) => revision === 1,
        );
        if (advanced.length !== 1 || held.length !== 1) {
          throw new Error(`adoption revisions were ${JSON.stringify(afterAdopt)}`);
        }
        if ((await providerCalls(manifest.provider_calls)) !== callsAfterGenerate) {
          throw new Error('preview or adoption called the provider');
        }
        report.tested.adopt_changes_one_revision = true;
        report.page_flow = 'passed';

        if ((await buttonCount(page, ASK)) > 0) {
          const denied = page.waitForResponse(
            (response) =>
              response.request().method() === 'POST' &&
              response.url().startsWith(manifest.api_base) &&
              !new URL(response.url()).pathname.endsWith('/preview') &&
              !new URL(response.url()).pathname.endsWith('/adoptions'),
          );
          await page.getByRole('button', { name: ASK }).click();
          const deniedResponse = await denied;
          const deniedBody = await deniedResponse.json().catch(() => ({}));
          const callsAfterDenial = await providerCalls(manifest.provider_calls);
          if (!isQuotaDenial(deniedResponse.status(), deniedBody)) {
            report.pending.push('quota_refusal_zero_extra_calls');
          } else if (callsAfterDenial !== callsAfterGenerate) {
            throw new Error('quota refusal made an extra provider call');
          } else {
            report.tested.quota_refusal_zero_extra_calls = true;
          }
        } else {
          report.pending.push('quota_refusal_zero_extra_calls');
        }
      }
    }

    if (screenshotPath) {
      await page.screenshot({ path: screenshotPath, fullPage: true });
      report.screenshot = true;
    }
    console.log(JSON.stringify(report, null, 2));
  } finally {
    await browser.close();
  }
} finally {
  api.kill('SIGTERM');
}
