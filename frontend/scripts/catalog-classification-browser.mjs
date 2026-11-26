// Run directly to start the isolated pytest fixture, or receive its manifest on stdin.
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { chromium, expect } from '@playwright/test';

const catalog = '/api/v1/catalog';

async function exercise({ origin, locale, admin, project, selected, other, target, label }) {
  const browser = await chromium.launch({ args: ['--no-proxy-server'] });
  try {
    const page = await browser.newPage({ baseURL: origin, locale });
    page.setDefaultTimeout(12_000);
    const english = locale === 'en-US';
    const panel = page.getByRole('region', {
      name: english ? 'Peri target suggestions' : 'Peri 靶标建议',
    });
    const ask = panel.getByRole('button', {
      name: english ? 'Ask Peri' : '询问 Peri', exact: true,
    });
    const previewButton = panel.getByRole('button', {
      name: english ? 'Preview selected' : '预览所选建议',
    });
    const adoptButton = panel.getByRole('button', {
      name: english ? 'Adopt previewed changes' : '采用预览中的变更',
    });
    const writes = [];
    page.on('request', (req) => {
      const url = new URL(req.url());
      if (req.method() === 'POST' && url.origin === origin && url.pathname.startsWith(catalog + '/')) {
        writes.push(url.pathname);
      }
    });

    async function get(path) {
      const response = await page.request.get(path);
      assert.equal(response.status(), 200, path);
      return response.json();
    }

    async function heads() {
      const { experiments } = await get(catalog + '/projects/' + project + '/experiments');
      assert.equal(experiments.length, 2);
      return Object.fromEntries(experiments.map((row) => [row.experiment_id, row]));
    }

    async function providerCount() {
      return (await get('/__fixture/provider-count')).count;
    }

    async function clickPost(button, path, status = 200) {
      const [response] = await Promise.all([
        page.waitForResponse((r) => r.request().method() === 'POST' && r.url() === origin + path),
        button.click(),
      ]);
      const body = await response.json();
      assert.equal(response.status(), status, JSON.stringify(body));
      return { request: response.request(), body };
    }

    const generationPath = (row) =>
      catalog + '/projects/' + project + '/experiments/' + row.experiment_id + '/classification-suggestions';

    async function prepare(row) {
      const experiments = page.getByRole('list', { name: /^(Experiments|实验)$/ });
      await experiments.getByRole('button', { name: new RegExp(row.display_name) }).click();
      const detail = page.getByRole('region', { name: /^(Experiment revision|实验修订)$/ });
      await expect(detail.getByRole('heading', { name: row.display_name })).toBeVisible();
      await expect(ask).toBeDisabled();
      await panel.getByLabel(english ? 'Target clue' : '靶标线索').fill(label);
      await panel.getByLabel(english ? 'Candidate target' : '候选靶标').selectOption(target.term_revision_id);
    }

    await page.goto('/');
    await page.locator('#login-username').fill(admin.username);
    await page.locator('#login-password').fill(admin.password);
    await page.getByRole('button', { name: /^(Sign in|登录)$/ }).click();
    await page.waitForURL(/\/workflows$/);
    const document = await page.goto('/catalog');
    assert.equal(document.status(), 200);
    assert.match(document.headers()['content-type'], /text\/html/);
    const search = page.getByRole('form', { name: /^(Find experiments|查找实验)$/ });
    await search.getByLabel(/^(Project|课题|项目)/).selectOption(project);
    await search.getByRole('button', { name: /^(Find experiments|查找实验)$/ }).click();

    const before = await heads();
    for (const row of [selected, other]) {
      assert.equal(before[row.experiment_id].experiment_revision_id, row.experiment_revision_id);
      assert.equal(before[row.experiment_id].revision_number, 1);
      assert.equal(before[row.experiment_id].target.status, 'unknown');
    }
    await prepare(selected);
    assert.equal(await providerCount(), 0);
    const generated = await clickPost(ask, generationPath(selected));
    assert.deepEqual(generated.request.postDataJSON(), {
      experiment_revision_id: selected.experiment_revision_id,
      field: 'target',
      reason_code: 'TARGET_TERM_MATCH',
      clue: label,
      vocabulary_term_revision_id: target.term_revision_id,
    });
    const generateKey = await generated.request.headerValue('Idempotency-Key');
    assert.ok(generateKey);
    assert.equal(generated.body.outcome, 'proposal');
    const proposal = generated.body.proposal;
    assert.equal(proposal.project_id, project);
    assert.equal(proposal.source, 'peri_bound_model');
    assert.equal(proposal.proposal_state, 'pending');
    assert.equal(proposal.items.length, 1);
    const [item] = proposal.items;
    assert.equal(item.target_id, selected.experiment_id);
    assert.equal(item.base_revision_id, selected.experiment_revision_id);
    assert.deepEqual(item.allowed_fields, ['target']);
    assert.equal(item.after.target.term_revision_id, target.term_revision_id);
    assert.deepEqual(item.citations.reason_codes, ['TARGET_TERM_MATCH']);
    assert.deepEqual(Object.fromEntries(item.citations.source_refs.map((ref) => [ref.kind, ref.ref_id])), {
      experiment_revision: selected.experiment_revision_id,
      term_revision: target.term_revision_id,
    });
    assert.equal(await providerCount(), 1);
    assert.deepEqual(await heads(), before);

    const replay = await page.request.post(generationPath(selected), {
      headers: {
        'Content-Type': 'application/json',
        'Idempotency-Key': generateKey,
        'X-CSRF-Token': await generated.request.headerValue('X-CSRF-Token'),
      },
      data: generated.request.postData(),
    });
    assert.equal(replay.status(), 200);
    assert.deepEqual(await replay.json(), generated.body);
    assert.equal(await providerCount(), 1);

    const proposalPath = catalog + '/proposals/' + proposal.proposal_id;
    const choice = panel.getByRole('checkbox', { name: new RegExp(label) });
    await expect(choice).not.toBeChecked();
    await expect(previewButton).toBeDisabled();
    await expect(adoptButton).toHaveCount(0);
    await choice.check();
    await expect(previewButton).toBeEnabled();
    await choice.uncheck();
    await expect(previewButton).toBeDisabled();
    await expect(adoptButton).toHaveCount(0);
    await choice.check();
    assert.deepEqual(writes, [generationPath(selected)]);
    const previewed = await clickPost(previewButton, proposalPath + '/preview');
    const ids = [item.item_id];
    assert.deepEqual(previewed.request.postDataJSON(), { selected_item_ids: ids });
    const preview = previewed.body;
    assert.equal(preview.has_conflicts, false);
    assert.deepEqual(preview.selected_item_ids, ids);
    assert.deepEqual(preview.items.map((change) => change.target_id), [selected.experiment_id]);
    assert.deepEqual(preview.expected_revisions, { [selected.experiment_id]: 1 });
    assert.ok(preview.preview_digest);
    await expect(panel.getByRole('listitem')).toHaveText(new RegExp(label));
    await expect(adoptButton).toBeEnabled();
    assert.deepEqual(await heads(), before);
    assert.equal(await providerCount(), 1);
    const screenshot = process.env.CLASSIFICATION_SCREENSHOT;
    if (screenshot && (english || screenshot.includes('{locale}'))) {
      await panel.screenshot({ path: screenshot.replace('{locale}', locale) });
    }

    const adoption = await clickPost(adoptButton, proposalPath + '/adoptions');
    assert.deepEqual(adoption.request.postDataJSON(), {
      selected_item_ids: ids,
      expected_revisions: preview.expected_revisions,
      preview_digest: preview.preview_digest,
    });
    const receipt = adoption.body;
    assert.ok(receipt.receipt_id);
    assert.equal(receipt.proposal_id, proposal.proposal_id);
    assert.equal(receipt.preview_digest, preview.preview_digest);
    assert.deepEqual(receipt.selected_item_ids, ids);
    assert.deepEqual(receipt.applied_targets.map((change) => change.target_id), [selected.experiment_id]);
    const after = await heads();
    assert.deepEqual(after[other.experiment_id], before[other.experiment_id]);
    const revised = after[selected.experiment_id];
    assert.notEqual(revised.experiment_revision_id, selected.experiment_revision_id);
    assert.equal(revised.experiment_revision_id, receipt.applied_targets[0].result_revision_id);
    assert.equal(revised.revision_number, 2);
    assert.equal(revised.target.status, 'known');
    assert.equal(revised.target.term_revision_id, target.term_revision_id);
    await expect(page.getByRole('status').filter({ hasText: receipt.receipt_id }).first()).toBeVisible();
    assert.equal(await providerCount(), 1);

    await prepare(other);
    const denied = await clickPost(ask, generationPath(other), 409);
    const deniedKey = await denied.request.headerValue('Idempotency-Key');
    assert.ok(deniedKey && deniedKey !== generateKey);
    assert.equal(denied.request.postDataJSON().experiment_revision_id, other.experiment_revision_id);
    assert.equal(denied.body.issues[0].code, 'catalog_model_quota_exhausted');
    await expect(panel.getByRole('alert')).toContainText(english ? /remaining calls|quota/i : /额度/);
    assert.equal(await providerCount(), 1);
    assert.deepEqual(await heads(), after);
    assert.deepEqual(writes, [
      generationPath(selected), proposalPath + '/preview', proposalPath + '/adoptions', generationPath(other),
    ]);
    console.log(JSON.stringify({
      locale,
      source: proposal.source,
      projection: generated.body.projection,
      providerCalls: await providerCount(),
      selectedRevision: revised.revision_number,
      otherRevision: after[other.experiment_id].revision_number,
      denial: denied.body.issues[0].code,
      page_flow: 'passed',
    }));
  } finally {
    await browser.close();
  }
}

if (process.argv.includes('--fixture-stdin')) {
  await exercise(JSON.parse(readFileSync(0, 'utf8')));
} else {
  const root = fileURLToPath(new URL('../..', import.meta.url));
  const result = spawnSync(process.env.HELIXWEAVE_PYTHON ?? 'python3', [
    '-I', '-S', 'scripts/checkout_bootstrap.py', '--repository-root', '.', 'pytest', '-qs',
    'test/api/catalog_classification_browser.py',
  ], { cwd: root, stdio: 'inherit' });
  if (result.error) throw result.error;
  process.exitCode = result.status ?? 1;
}
