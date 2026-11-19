import { expect, test, type Route } from '@playwright/test';
import { createAuthoringSchemaFixture, WORKFLOW_ID } from '../src/features/input-workbench/test-fixtures';

const PROJECT = `prj_${'1'.repeat(32)}`;
const SID = `aas_${'1'.repeat(32)}`;
const REV_A = `smpr_${'a'.repeat(32)}`;
const REV_B = `smpr_${'b'.repeat(32)}`;
const REV_OTHER = `smpr_${'c'.repeat(32)}`;
const EXPR = `expr_${'e'.repeat(32)}`;
const REFERENCE = `refpr_${'1'.repeat(32)}`;
const SNAPSHOT_ID = 'vsnap_0123456789abcdef0123456789abcdef';
const PINNED = [REV_A, REV_B];

const PROFILE = {
  profile_id: `refp_${'1'.repeat(32)}`,
  revision_id: REFERENCE,
  revision_number: 1,
  display_name: 'Synthetic reference',
  organism: 'synthetic',
  assembly: 'test',
  identity_sha256: 'ab'.repeat(32),
};

function setsEqual(left: readonly string[], right: readonly string[]): boolean {
  if (left.length === 0 || left.length !== right.length) return false;
  const pending = [...right];
  for (const value of left) {
    const index = pending.indexOf(value);
    if (index < 0) return false;
    pending.splice(index, 1);
  }
  return pending.length === 0;
}

async function fulfill(route: Route, status: number, body: unknown) {
  await route.fulfill({
    status,
    contentType: 'application/json',
    body: JSON.stringify(body),
  });
}

test('synthetic draft carries context and binds only an equal sample set', async ({ page }) => {
  const savedIds = [REV_A, REV_B];
  const validateBodies: Array<Record<string, unknown>> = [];
  const forbidden: string[] = [];
  let frozenBindings = 0;

  await page.route('**/api/v1/**', async (route) => {
    const request = route.request();
    const url = new URL(request.url());
    const path = url.pathname;
    const method = request.method();
    if (
      (method === 'PUT' && path.endsWith('/draft')) ||
      (method === 'POST' && path.endsWith(`/authoring-sessions/${SID}/validate`)) ||
      (method === 'POST' && path.endsWith('/runs'))
    ) {
      forbidden.push(`${method} ${path}`);
      await fulfill(route, 500, { ok: false });
      return;
    }
    if (path === '/api/v1/auth/session' && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        setup_required: false,
        authenticated: true,
        principal: { user_id: 'user-1', username: 'ada', role: 'admin' },
      });
      return;
    }
    if (path === '/api/v1/projects' && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        projects: [{ project_id: PROJECT, display_name: 'Synthetic project' }],
      });
      return;
    }
    if (path === '/api/v1/agent/authoring-sessions' && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        sessions: [{
          sid: SID,
          project_id: PROJECT,
          workflow_id: WORKFLOW_ID,
          draft_state: 'active',
          draft_revision: 3,
          expires_at: '2026-11-17T00:00:00Z',
          purge_after: '2026-11-24T00:00:00Z',
          connection: { state: 'closed' },
          requires_new_connection: true,
        }],
        next_cursor: null,
      });
      return;
    }
    if (path === `/api/v1/agent/authoring-sessions/${SID}` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        sid: SID,
        project_id: PROJECT,
        workflow_id: WORKFLOW_ID,
        draft_state: 'active',
        draft_revision: 3,
        expires_at: '2026-11-17T00:00:00Z',
        purge_after: '2026-11-24T00:00:00Z',
        schema_digest: 'a'.repeat(64),
        projection_digest: 'b'.repeat(64),
        draft: {
          config: { threads: 4, outdir: 'results' },
          options: { strict_inputs: false },
          samples: [{ sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' }],
          sample_revision_ids: savedIds,
          input_selections: [],
          reference_profile_revision_id: REFERENCE,
        },
        authoring_capability: {
          available: true,
          projection_version: 'hw-a1/1',
          allowed_intents: ['suggest', 'explain_issues'],
          allowed_issue_codes: [],
          fields: [{
            token: 'resource.threads',
            label: 'Threads',
            type: 'integer',
            constraints: { type: 'integer', minimum: 1, maximum: 16 },
          }],
        },
        connection: { state: 'closed' },
        requires_new_connection: true,
      });
      return;
    }
    if (path === `/api/v1/agent/authoring-sessions/${SID}/receipts` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        sid: SID,
        draft_state: 'active',
        receipts: [],
        next_cursor: null,
      });
      return;
    }
    if (path === `/api/v1/projects/${PROJECT}/samples` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        project_id: PROJECT,
        samples: [],
        next_cursor: null,
      });
      return;
    }
    if (path === `/api/v1/projects/${PROJECT}/input-revisions` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        project_id: PROJECT,
        revisions: [],
        next_cursor: null,
      });
      return;
    }
    if (path === `/api/v1/workflows/${WORKFLOW_ID}` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        workflow_id: WORKFLOW_ID,
        workflow: {
          metadata: {
            workflow_id: WORKFLOW_ID,
            name: 'Synthetic workflow',
            version: '0.0.0',
            description: 'Synthetic workflow for experiment binding.',
            engines: [],
            tags: [],
          },
          schema_version: '1.1.0',
          capabilities: { supports: [] },
          upstream_identity: null,
          availability: {
            authoring: 'available',
            execution: 'available',
            reason_code: 'WORKFLOW_EXECUTION_READY',
          },
        },
        issues: [],
      });
      return;
    }
    if (path === `/api/v1/workflows/${WORKFLOW_ID}/schema` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        workflow_id: WORKFLOW_ID,
        schema: createAuthoringSchemaFixture(),
        issues: [],
      });
      return;
    }
    if (path === `/api/v1/workflows/${WORKFLOW_ID}/reference-profiles` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        workflow_id: WORKFLOW_ID,
        profiles: [PROFILE],
        issues: [],
      });
      return;
    }
    if (path === `/api/v1/catalog/projects/${PROJECT}/experiments` && method === 'GET') {
      await fulfill(route, 200, {
        ok: true,
        project_id: PROJECT,
        revision_view: 'current',
        experiments: [{
          experiment_id: `exp_${'1'.repeat(32)}`,
          experiment_revision_id: EXPR,
          revision_number: 1,
          display_name: 'Synthetic experiment',
          samples: [
            { sample_revision_id: REV_A },
            { sample_revision_id: REV_B },
          ],
        }],
        next_cursor: null,
      });
      return;
    }
    if (path === `/api/v1/workflows/${WORKFLOW_ID}/validate` && method === 'POST') {
      const body = request.postDataJSON() as Record<string, unknown>;
      validateBodies.push(body);
      const ids = Array.isArray(body.sample_revision_ids)
        ? body.sample_revision_ids.filter((id): id is string => typeof id === 'string')
        : [];
      const experimentId = body.experiment_revision_id;
      if (typeof experimentId === 'string') {
        if (experimentId !== EXPR || !setsEqual(ids, PINNED)) {
          await fulfill(route, 409, {
            ok: false,
            issues: [{
              code: 'catalog_sample_set_mismatch',
              message: 'The sample revisions do not match the experiment.',
              severity: 'error',
            }],
          });
          return;
        }
        frozenBindings += 1;
        await fulfill(route, 200, validationEnvelope({
          status: 'experiment_snapshot',
          experiment_revision_id: experimentId,
          payload_digest: 'd'.repeat(64),
        }));
        return;
      }
      await fulfill(route, 200, validationEnvelope({ status: 'unresolved' }));
      return;
    }
    await fulfill(route, 404, { ok: false, path, method });
  });

  await page.goto(`/workflows/${WORKFLOW_ID}/configuration-drafts`);
  await page.getByRole('button', { name: /Synthetic project/ }).click();
  await page.getByRole('link', { name: 'Open experiment validation' }).click();
  await expect(page.getByTestId('carried-project-id')).toContainText(PROJECT);
  await expect(page.getByTestId('carried-reference-profile')).toContainText(REFERENCE);
  await expect(page.getByTestId('draft-review-json')).toContainText('"threads": 4');
  await expect(page.getByTestId('draft-review-json')).toContainText('S1');
  await expect(page.getByLabel('Submitted sample revision 1')).toHaveValue(REV_A);
  await expect(page.getByLabel('Submitted sample revision 2')).toHaveValue(REV_B);
  await expect(page.getByTestId('experiment-binding-status')).toHaveText(
    'No experiment binding. Choosing a revision or opening this page does not create one.',
  );
  expect(validateBodies).toHaveLength(0);

  await page.getByRole('button', { name: 'Reverse sample revision order' }).click();
  await expect(page.getByLabel('Submitted sample revision 1')).toHaveValue(REV_B);
  await expect(page.getByLabel('Submitted sample revision 2')).toHaveValue(REV_A);
  const experiment = page.locator('#experiment-revision-choice');
  await experiment.selectOption({ label: 'Synthetic experiment' });
  await expect(page.getByTestId('experiment-selection-preview')).toContainText(
    'Order does not matter. Validation can freeze this revision.',
  );
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  await expect(page.getByTestId('experiment-binding-status')).toHaveText(
    `Bound to experiment revision ${EXPR} for this validation snapshot.`,
  );
  expect(validateBodies.at(-1)?.sample_revision_ids).toEqual([REV_B, REV_A]);
  expect(validateBodies.at(-1)?.experiment_revision_id).toBe(EXPR);
  expect(frozenBindings).toBe(1);

  await experiment.selectOption({ label: 'Omit experiment' });
  await expect(page.getByTestId('experiment-binding-status')).toHaveText(
    'No experiment binding. Choosing a revision or opening this page does not create one.',
  );
  // The previous snapshot cannot describe this changed selection.
  await expect(page.getByTestId('create-validated-run-button')).toBeDisabled();
  expect.soft(await page.getByTestId('validation-feedback').innerText())
    .not.toContain('This exact draft can create one run.');
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  await expect(page.getByTestId('experiment-binding-status')).toHaveText(
    'Not bound. This validation did not select an experiment revision.',
  );
  expect(validateBodies.at(-1)?.experiment_revision_id).toBeUndefined();
  expect(frozenBindings).toBe(1);

  await page.getByLabel('Submitted sample revision 1').fill(REV_OTHER);
  await experiment.selectOption({ label: 'Synthetic experiment' });
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  await expect(page.getByTestId('experiment-binding-status')).toHaveText(
    'No experiment binding was saved.',
  );
  await expect(page.getByText('catalog_sample_set_mismatch')).toBeVisible();
  expect(validateBodies.at(-1)?.sample_revision_ids).toEqual([REV_OTHER, REV_A]);
  expect(frozenBindings).toBe(1);
  expect(forbidden).toEqual([]);
  expect(savedIds).toEqual([REV_A, REV_B]);

  // A lost response leaves the client unable to confirm persistence.
  await page.getByLabel('Submitted sample revision 1').fill(REV_B);
  await page.route('**/validate', route => route.abort('connectionreset'));
  await page.getByRole('button', { name: 'Validate current inputs' }).click();
  await expect(page.getByText('VALIDATION_UNAVAILABLE')).toBeVisible();
  expect.soft(await page.getByTestId('experiment-binding-status').innerText())
    .not.toContain('No experiment binding was saved.');

  // Inspect the persisted draft in Config and then return to Review.
  await page.getByRole('tab', { name: /Config/ }).click();
  await page.waitForURL('**/new-run?step=config');
  await page.getByRole('tab', { name: /Review/ }).click();
  await page.waitForURL('**/new-run?step=review');
  expect.soft(await experiment.count(), 'experiment selector after tab navigation').toBe(1);
  expect.soft(await page.getByTestId('draft-review-json').allTextContents()).toEqual([
    expect.stringContaining('S1'),
  ]);

  // Enter ordinary validation directly, with an active project available.
  await page.goto(`/workflows/${WORKFLOW_ID}/new-run`);
  await expect(page.getByRole('heading', { name: 'Input workbench' })).toBeVisible();
  expect.soft(await experiment.count(), 'experiment selector at ordinary entry').toBe(1);

});

function validationEnvelope(binding: { status: 'unresolved' } | {
  status: 'experiment_snapshot';
  experiment_revision_id: string;
  payload_digest: string;
}) {
  return {
    ok: true,
    workflow_id: WORKFLOW_ID,
    value: null,
    snapshot: {
      snapshot_id: SNAPSHOT_ID,
      workflow_id: WORKFLOW_ID,
      schema_version: '1.1.0',
      adapter_version: '0.3.0',
      payload_digest: 'a'.repeat(64),
      validated_at: '2026-07-14T00:00:00.000Z',
      expires_at: '2026-07-14T00:30:00.000Z',
      project_id: PROJECT,
      binding_mode: 'bound_v1',
      provenance: 'resolved',
      sample_revision_ids: PINNED,
      binding_digest: 'b'.repeat(64),
      input_binding: {
        mode: 'compatibility_unresolved_v1',
        adapter_contract_version: null,
        digest: 'c'.repeat(64),
        fully_managed: false,
        input_uses: [],
      },
      reference_profile: PROFILE,
      experiment_binding: binding,
    },
    issues: [],
  };
}
