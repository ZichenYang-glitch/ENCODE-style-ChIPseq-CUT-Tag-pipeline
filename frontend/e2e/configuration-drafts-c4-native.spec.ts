import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdirSync, readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

import { expect, test, type Page, type Response } from '@playwright/test';

import { login } from './auth';

const runtimeRoot = process.env.ENCODE_PIPELINE_E2E_ROOT!;
const evidenceRoot = process.env.HELIXWEAVE_C4_EVIDENCE!;
const binaryPath = process.env.HELIXWEAVE_PERI_BINARY!;
const acceptedDigest =
  '5cc82caf6e04578f19ff5790fbc4eb8c9e115793ff6950a699c9e2ad28c50f47';

test.describe.configure({ mode: 'serial' });

test.beforeAll(() => {
  const digest = createHash('sha256').update(readFileSync(binaryPath)).digest('hex');
  expect(digest).toBe(acceptedDigest);
  mkdirSync(evidenceRoot, { recursive: true });
});

function workflowId(): string {
  const manifest = JSON.parse(readFileSync(join(runtimeRoot, 'runtime.json'), 'utf8')) as {
    workflowId: string;
  };
  return manifest.workflowId;
}

function seedProject(): void {
  execFileSync('python3', [join(process.cwd(), 'e2e/seed_authoring_fixture.py')], {
    env: process.env,
    stdio: 'pipe',
  });
}

function record(name: string, value: unknown): void {
  writeFileSync(join(evidenceRoot, name), `${JSON.stringify(value, null, 2)}\n`);
}

async function openDrafts(page: Page) {
  const heading = page.getByRole('heading', { level: 1, name: 'Configuration drafts' });
  const username = page.getByLabel('Username');
  await page.goto(`/workflows/${workflowId()}/configuration-drafts`);
  const alreadyThere = await heading
    .waitFor({ state: 'visible', timeout: 5000 })
    .then(() => true)
    .catch(() => false);
  if (!alreadyThere) {
    await expect(username).toBeVisible();
    await login(page);
    await page.goto(`/workflows/${workflowId()}/configuration-drafts`);
  }
  await expect(heading).toBeVisible();
}

async function createDraft(page: Page): Promise<void> {
  const created = page.waitForResponse(
    (response) =>
      response.url().includes('/api/v1/agent/authoring-sessions') &&
      response.request().method() === 'POST' &&
      !response.url().includes('/draft') &&
      !response.url().includes('/validate') &&
      !response.url().includes('/connection') &&
      !response.url().includes('/turns'),
  );
  await page.getByRole('button', { name: 'New draft' }).click();
  expect((await created).ok()).toBeTruthy();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
}

async function saveCompleteDraft(page: Page): Promise<Response> {
  const threads = page.getByRole('spinbutton', { name: /thread/i });
  await threads.click();
  await threads.fill('');
  await threads.pressSequentially('4');
  const saved = page.waitForResponse(
    (response) => response.url().includes('/draft') && response.request().method() === 'PUT',
  );
  await page.getByRole('button', { name: 'Save edits' }).click();
  expect((await saved).ok()).toBeTruthy();

  await page.getByLabel('Add a sample revision to the plan').selectOption({ label: 'Browser sample' });
  await page.getByLabel('Browser sample assay').selectOption('chipseq');
  await page.getByLabel('Browser sample target').fill('CTCF');
  await page.getByLabel('Browser sample peak_mode').selectOption('narrow');
  await page.getByLabel('Browser sample layout').selectOption('SE');
  await page.getByLabel('Browser sample role').selectOption('treatment');
  await page.getByLabel('Browser sample read 1').selectOption({ label: 'browser-sample-r1' });
  const reference = page.getByRole('combobox', { name: 'Reference profile' });
  const referenceValue = await reference.locator('option').nth(1).getAttribute('value');
  expect(referenceValue).toBeTruthy();
  await reference.selectOption(referenceValue!);
  const selectionSaved = page.waitForResponse(
    (response) => response.url().includes('/draft') && response.request().method() === 'PUT',
  );
  await page.getByRole('button', { name: 'Save edits' }).click();
  const response = await selectionSaved;
  expect(response.ok()).toBeTruthy();
  return response;
}

test('a disconnected assistant still saves and validates the current revision @desktop', async ({
  page,
}) => {
  seedProject();
  await openDrafts(page);
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  await createDraft(page);
  const saved = await saveCompleteDraft(page);
  const savedBody = (await saved.json()) as {
    draft_revision: number;
    draft: {
      config: { threads?: number };
      workflow_samples: {
        assay: string;
        target: string;
        reads: { input_file_revision_id: string }[];
      }[];
      reference_profile_revision_id: string | null;
    };
  };
  expect(savedBody.draft.config.threads).toBe(4);
  expect(savedBody.draft.workflow_samples[0]?.assay).toBe('chipseq');
  expect(savedBody.draft.workflow_samples[0]?.target).toBe('CTCF');
  expect(savedBody.draft.workflow_samples[0]?.reads[0]?.input_file_revision_id).toMatch(
    /^inpfr_[0-9a-f]{32}$/,
  );
  expect(savedBody.draft.reference_profile_revision_id).toMatch(/^refpr_[0-9a-f]{32}$/);

  await page.reload();
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  await expect(page.getByRole('spinbutton', { name: /thread/i })).toHaveValue('4');
  await expect(page.getByLabel('Browser sample target')).toHaveValue('CTCF');
  await expect(page.getByText('Assistant available')).toHaveCount(0);

  const validated = page.waitForResponse(
    (response) => response.url().includes('/validate') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Validate draft' }).click();
  const validation = await validated;
  expect(validation.ok()).toBeTruthy();
  const validationBody = (await validation.json()) as {
    ok: boolean;
    draft_revision: number;
    workflow_id?: string;
    snapshot?: {
      snapshot_id?: string;
      schema_version?: string;
      sample_revision_ids?: string[];
      reference_profile?: { revision_id?: string } | null;
    } | null;
  };
  expect(validationBody.ok).toBe(true);
  expect(validationBody.draft_revision).toBe(savedBody.draft_revision);
  expect(validationBody.workflow_id).toBe(workflowId());
  expect(validationBody.snapshot?.snapshot_id).toMatch(/^vsnap_[0-9a-f]{32}$/);
  expect(validationBody.snapshot?.schema_version).toBeTruthy();
  expect(validationBody.snapshot?.sample_revision_ids?.length).toBeGreaterThan(0);
  expect(validationBody.snapshot?.reference_profile?.revision_id).toBe(
    savedBody.draft.reference_profile_revision_id,
  );
  record('validation.json', {
    draft_revision: validationBody.draft_revision,
    workflow_id: validationBody.workflow_id,
    snapshot_id: validationBody.snapshot?.snapshot_id,
    schema_version: validationBody.snapshot?.schema_version,
    sample_revision_ids: validationBody.snapshot?.sample_revision_ids,
    reference_profile_revision_id: validationBody.snapshot?.reference_profile?.revision_id,
  });
  await page.screenshot({ path: join(evidenceRoot, 'validated-draft.png') });
});

test('an explicit connection adopts one native proposal and survives a new login @desktop', async ({
  page,
}) => {
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  const connection = page.waitForResponse(
    (response) => response.url().includes('/connection') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Connect assistant' }).click();
  expect((await connection).ok()).toBeTruthy();
  await expect(page.getByText('Assistant available')).toBeVisible();
  const scope = admittedScopeEvidence();
  expect(scope.owner_outside).toBe(true);
  expect(scope.memory_max).toBe('536870912');
  expect(scope.swap_max).toBe('0');
  expect(scope.pids_max).toBe('128');
  expect(scope.cpu_max).toBe('100000 100000');
  record('native-scope.json', { ...scope, visible_processes: processesInNativeScope() });

  await page.getByRole('radio', { name: 'Suggest values' }).check();
  await page.getByRole('checkbox', { name: 'Threads' }).check();
  await page.getByLabel('Preferred value').fill('2');
  const turn = page.waitForResponse(
    (response) => response.url().includes('/turns') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Request suggestion' }).click();
  const turnResponse = await turn;
  const turnText = await turnResponse.text();
  record('turn-response.json', { status: turnResponse.status(), body: turnText });
  expect(turnResponse.ok()).toBeTruthy();
  const turnBody = JSON.parse(turnText) as {
    status?: string;
    reason_code?: string | null;
    result_kind?: string | null;
    proposal_digest?: string;
    proposal_id?: string;
    adoptable?: boolean;
  };
  expect(turnBody.result_kind).toBe('proposal');
  expect(turnBody.proposal_digest).toMatch(/^[0-9a-f]{64}$/);
  expect(turnBody.adoptable).toBe(true);
  const toolLog = readFileSync(join(runtimeRoot, 'logs', 'peri-owner.log'), 'utf8');
  expect(toolLog).toContain('peri-authoring-tool get_authoring_schema');
  expect(toolLog).toContain('peri-authoring-tool get_draft_summary');
  record('proposal.json', {
    result_kind: turnBody.result_kind,
    proposal_id: turnBody.proposal_id,
    proposal_digest: turnBody.proposal_digest,
    tools: ['get_authoring_schema', 'get_draft_summary'],
  });

  const adopted = page.waitForResponse(
    (response) => response.url().includes('/adopt') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Adopt suggestion' }).click();
  const adoption = await adopted;
  expect(adoption.ok()).toBeTruthy();
  const adoptionRequest = adoption.request().postDataJSON() as {
    proposal_digest?: string;
    expected_revision?: number;
  };
  expect(adoptionRequest.proposal_digest).toBe(turnBody.proposal_digest);
  const adoptionBody = (await adoption.json()) as {
    draft_revision: number;
    adoption_receipt?: {
      proposal_id?: string;
      base_revision?: number;
      result_revision?: number;
      changed?: boolean;
    };
  };
  expect(adoptionBody.adoption_receipt?.proposal_id).toBe(turnBody.proposal_id);
  expect(adoptionBody.adoption_receipt?.changed).toBe(true);
  expect(adoptionBody.adoption_receipt?.result_revision).toBeGreaterThan(
    adoptionBody.adoption_receipt?.base_revision ?? -1,
  );
  record('adoption.json', adoptionBody);
  await expect(page.getByText('Suggestion adopted.')).toBeVisible();
  const adoptedValidation = page.waitForResponse(
    (response) => response.url().includes('/validate') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Validate draft' }).click();
  const adoptedValidationResponse = await adoptedValidation;
  expect(adoptedValidationResponse.ok()).toBeTruthy();
  const adoptedValidationBody = (await adoptedValidationResponse.json()) as {
    ok: boolean;
    draft_revision: number;
    snapshot?: { snapshot_id?: string; schema_version?: string } | null;
  };
  expect(adoptedValidationBody.ok).toBe(true);
  expect(adoptedValidationBody.draft_revision).toBe(adoptionBody.draft_revision);
  expect(adoptedValidationBody.draft_revision).toBe(
    adoptionBody.adoption_receipt?.result_revision,
  );
  expect(adoptedValidationBody.snapshot?.snapshot_id).toMatch(/^vsnap_[0-9a-f]{32}$/);
  expect(adoptedValidationBody.snapshot?.schema_version).toBeTruthy();
  record('adopted-validation.json', adoptedValidationBody);
  await expect(page.getByRole('article', { name: 'Adopted record' })).toBeVisible();
  await expect(page.getByRole('spinbutton', { name: /thread/i })).toHaveValue('2');
  await page.screenshot({ path: join(evidenceRoot, 'adopted-proposal.png') });

  await page.getByTestId('current-username').click();
  const loggedOut = page.waitForResponse(
    (response) =>
      response.url().includes('/api/v1/auth/logout') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Sign out' }).click();
  expect((await loggedOut).ok()).toBeTruthy();
  await page.goto('/login');
  await expect(page.getByLabel('Username')).toBeVisible();
  await login(page);
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  await expect(page.getByRole('article', { name: 'Adopted record' })).toBeVisible();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  await expect(page.getByText(/assistant is not connected/i)).toBeVisible();
  await page.screenshot({ path: join(evidenceRoot, 'receipt-after-login.png') });
});

test('an explanation and a closed connection do not create another adoptable proposal @desktop', async ({
  page,
}) => {
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  await page.getByRole('button', { name: 'Connect assistant' }).click();
  await expect(page.getByText('Assistant available')).toBeVisible();
  await page.getByRole('radio', { name: 'Explain a validation issue' }).check();
  const explained = page.waitForResponse(
    (response) => response.url().includes('/turns') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Request suggestion' }).click();
  const explanation = await explained;
  expect(explanation.ok()).toBeTruthy();
  const explanationBody = (await explanation.json()) as { result_kind?: string; adoptable?: boolean };
  expect(explanationBody.result_kind).toBe('explanation');
  expect(explanationBody.adoptable).toBe(false);
  await expect(page.getByText(/cannot be adopted/i)).toBeVisible();
  await expect(page.getByRole('button', { name: 'Adopt suggestion' })).toHaveCount(0);

  await page.getByRole('button', { name: 'Close assistant' }).click();
  await expect(page.getByText('Assistant available')).toHaveCount(0);
  let turnsAfterClose = 0;
  page.on('request', (request) => {
    if (request.method() === 'POST' && request.url().includes('/turns')) turnsAfterClose += 1;
  });
  await expect(page.getByRole('button', { name: 'Request suggestion' })).toBeDisabled();
  await expect(page.getByText(/Connect the assistant before requesting/i)).toBeVisible();
  expect(turnsAfterClose).toBe(0);
});

test('unsaved edits survive cancel, a revision conflict, and a completed turn whose response is hidden @desktop', async ({
  page,
}) => {
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  const threads = page.getByRole('spinbutton', { name: /thread/i });
  await threads.click();
  await threads.fill('');
  await threads.pressSequentially('7');
  await expect(page.getByText(/you have unsaved edits/i)).toBeVisible();
  await page.getByRole('link', { name: 'Runs' }).click();
  await expect(page.getByRole('dialog', { name: 'Unsaved edits' })).toBeVisible();
  await page.getByRole('button', { name: 'Cancel' }).click();
  await expect(threads).toHaveValue('7');

  const saved = page.waitForResponse(
    (response) => response.url().includes('/draft') && response.request().method() === 'PUT',
  );
  await page.getByRole('button', { name: 'Save edits' }).click();
  const saveResponse = await saved;
  expect(saveResponse.ok()).toBeTruthy();
  const saveBody = (await saveResponse.json()) as { draft_revision: number };
  const staleBody = saveResponse.request().postDataJSON() as {
    expected_revision: number;
    config: { threads?: number };
  };
  staleBody.expected_revision = saveBody.draft_revision;
  staleBody.config.threads = 6;
  const csrf = (await page.context().cookies()).find((cookie) => cookie.name === 'helixweave_csrf');
  const conflicted = await page.request.put(saveResponse.url(), {
    data: staleBody,
    headers: csrf ? { 'X-CSRF-Token': csrf.value } : {},
  });
  expect(conflicted.ok()).toBeTruthy();
  await threads.click();
  await threads.fill('');
  await threads.pressSequentially('9');
  await expect(page.getByRole('button', { name: 'Save edits' })).toBeEnabled();
  const rejected = page.waitForResponse(
    (response) => response.url().includes('/draft') && response.request().method() === 'PUT',
  );
  await page.getByRole('button', { name: 'Save edits' }).click();
  expect((await rejected).status()).toBe(409);
  await expect(page.getByText(/nothing was overwritten/i)).toBeVisible();
  await expect(threads).toHaveValue('9');
  await page.reload();
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  await expect(page.getByRole('spinbutton', { name: /thread/i })).toHaveValue('6');
  await page.getByRole('button', { name: 'Connect assistant' }).click();
  await expect(page.getByText('Assistant available')).toBeVisible();
  let hidden = false;
  let lostTurnId = '';
  let completedTurn = '';
  const beforeTools = toolCounts();
  await page.route('**/turns', async (route) => {
    if (hidden || route.request().method() !== 'POST') {
      await route.continue();
      return;
    }
    hidden = true;
    lostTurnId = String(
      (route.request().postDataJSON() as { client_turn_id?: string }).client_turn_id ?? '',
    );
    const completed = await route.fetch();
    completedTurn = await completed.text();
    await route.abort('failed');
  });
  await page.getByRole('radio', { name: 'Suggest values' }).check();
  await page.getByRole('checkbox', { name: 'Threads' }).check();
  await page.getByLabel('Preferred value').fill('3');
  await page.getByRole('button', { name: 'Request suggestion' }).click();
  await expect(page.getByText(/suggestion response was not confirmed/i).first()).toBeVisible();
  const completedBody = JSON.parse(completedTurn) as {
    client_turn_id?: string;
    proposal_digest?: string;
    proposal_id?: string;
  };
  const afterCommit = toolCounts();
  expect(afterCommit.schema).toBe(beforeTools.schema + 1);
  expect(afterCommit.summary).toBe(beforeTools.summary + 1);
  expect(completedBody.client_turn_id).toBe(lostTurnId);
  expect(completedBody.proposal_digest).toMatch(/^[0-9a-f]{64}$/);
  const replayed = page.waitForResponse(
    (response) => response.url().includes('/turns') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Check this request' }).click();
  const replay = await replayed;
  expect(replay.ok()).toBeTruthy();
  const replayRequest = replay.request().postDataJSON() as { client_turn_id?: string };
  const replayBody = (await replay.json()) as {
    client_turn_id?: string;
    proposal_digest?: string;
    proposal_id?: string;
  };
  expect(replayRequest.client_turn_id).toBe(lostTurnId);
  expect(replayBody.client_turn_id).toBe(completedBody.client_turn_id);
  expect(replayBody.proposal_id).toBe(completedBody.proposal_id);
  expect(replayBody.proposal_digest).toBe(completedBody.proposal_digest);
  expect(toolCounts()).toEqual(afterCommit);
  record('replayed-turn.json', {
    client_turn_id: replayBody.client_turn_id,
    proposal_id: replayBody.proposal_id,
    proposal_digest: replayBody.proposal_digest,
    native_schema_reads_before_retry: afterCommit.schema,
    native_schema_reads_after_retry: toolCounts().schema,
    native_summary_reads_before_retry: afterCommit.summary,
    native_summary_reads_after_retry: toolCounts().summary,
  });
});

test('rejecting a proposal leaves the saved revision unchanged @desktop', async ({ page }) => {
  await openDrafts(page);
  await page.getByRole('button', { name: /revision /i }).first().click();
  const connect = page.getByRole('button', { name: 'Connect assistant' });
  if (await connect.isEnabled()) {
    await connect.click();
  }
  await expect(page.getByText('Assistant available')).toBeVisible();
  const threads = page.getByRole('spinbutton', { name: /thread/i });
  if (!(await page.getByRole('button', { name: 'Reject suggestion' }).isVisible())) {
    await page.getByRole('radio', { name: 'Suggest values' }).check();
    await page.getByRole('checkbox', { name: 'Threads' }).check();
    await page.getByLabel('Preferred value').fill('5');
    const suggested = page.waitForResponse(
      (response) => response.url().includes('/turns') && response.request().method() === 'POST',
    );
    await page.getByRole('button', { name: 'Request suggestion' }).click();
    expect((await suggested).ok()).toBeTruthy();
  }
  const before = await threads.inputValue();
  const rejected = page.waitForResponse(
    (response) => response.url().includes('/reject') && response.request().method() === 'POST',
  );
  await page.getByRole('button', { name: 'Reject suggestion' }).click();
  const rejection = await rejected;
  expect(rejection.ok()).toBeTruthy();
  const rejectionBody = (await rejection.json()) as {
    proposal_state?: string;
    adoptable?: boolean;
    draft_revision?: number;
  };
  expect(rejectionBody.proposal_state).toBe('rejected');
  expect(rejectionBody.adoptable).toBe(false);
  await expect(page.getByText('Suggestion rejected. The draft was not changed.')).toBeVisible();
  await expect(threads).toHaveValue(before);
  record('rejection.json', { ...rejectionBody, threads: before });
  record('acceptance-scope.json', {
    browser_exercised: [
      'save_and_reread',
      'explicit_connect',
      'two_read_tools',
      'proposal_and_explanation',
      'adopt',
      'validate_adopted_revision',
      'receipt_after_login',
      'unsaved_navigation',
      'revision_conflict_409',
      'completed_turn_with_hidden_response',
      'reject',
    ],
    browser_not_exercised: {
      draft_ttl: 'AuthoringService._require_active and _adoptable; draft TTL is 7 days',
      proposal_ttl: 'AuthoringService._adoptable and repository.reject; proposal TTL is 30 minutes',
    },
  });
});

function toolCounts(): { schema: number; summary: number } {
  const log = readFileSync(join(runtimeRoot, 'logs', 'peri-owner.log'), 'utf8');
  return {
    schema: log.split('peri-authoring-tool get_authoring_schema').length - 1,
    summary: log.split('peri-authoring-tool get_draft_summary').length - 1,
  };
}

function admittedScopeEvidence(): {
  unit: string;
  pid: string;
  path: string;
  memory_max: string;
  swap_max: string;
  pids_max: string;
  cpu_max: string;
  owner_outside: boolean;
} {
  const log = readFileSync(join(runtimeRoot, 'logs', 'peri-owner.log'), 'utf8');
  const line = log.split('\n').find((item) => item.includes('peri-authoring-scope '));
  expect(line).toBeTruthy();
  const text = line ?? '';
  return {
    unit: scopeField(text, 'unit', 'pid'),
    pid: scopeField(text, 'pid', 'path'),
    path: scopeField(text, 'path', 'memory_max'),
    memory_max: scopeField(text, 'memory_max', 'swap_max'),
    swap_max: scopeField(text, 'swap_max', 'pids_max'),
    pids_max: scopeField(text, 'pids_max', 'cpu_max'),
    cpu_max: scopeField(text, 'cpu_max', 'owner_outside'),
    owner_outside: scopeField(text, 'owner_outside') === 'True',
  };
}

function scopeField(line: string, key: string, next?: string): string {
  const marker = `${key}=`;
  const start = line.indexOf(marker);
  if (start < 0) return '';
  const from = start + marker.length;
  if (!next) return line.slice(from).trim();
  const end = line.indexOf(` ${next}=`, from);
  return line.slice(from, end < 0 ? undefined : end);
}

function processesInNativeScope(): { pid: string; comm: string; path: string }[] {
  let names: string[] = [];
  try {
    names = readdirSync('/proc');
  } catch {
    return [];
  }
  const hits: { pid: string; comm: string; path: string }[] = [];
  for (const name of names) {
    if (!/^\d+$/.test(name)) continue;
    try {
      const cgroup = readFileSync(`/proc/${name}/cgroup`, 'utf8');
      const relative = cgroup
        .split('\n')
        .find((line) => line.startsWith('0::') && line.includes('hw-peri-a1-'))
        ?.slice('0::'.length);
      if (!relative) continue;
      hits.push({
        pid: name,
        comm: readFileSync(`/proc/${name}/comm`, 'utf8').trim(),
        path: join('/sys/fs/cgroup', relative.replace(/^\//, '')),
      });
    } catch {
      continue;
    }
  }
  return hits;
}
