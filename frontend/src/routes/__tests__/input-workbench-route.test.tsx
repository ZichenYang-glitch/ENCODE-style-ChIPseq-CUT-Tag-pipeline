import { act, fireEvent, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import type { SchemaResponse } from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { createStubWorkflowClient } from '../../api/client';
import { appRoutes } from '../../app/router';
import { createAuthoringSchemaFixture, WORKFLOW_ID } from '../../features/input-workbench/test-fixtures';
import { renderWithRouter } from '../../test/test-utils';

const generatedMocks = vi.hoisted(() => ({
  createRun: vi.fn(),
  getWorkflowSchema: vi.fn(),
  validateWorkflow: vi.fn(),
}));

beforeAll(async () => {
  await import('../workflows/new-run');
});

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: generatedMocks.getWorkflowSchema,
  validateWorkflow: generatedMocks.validateWorkflow,
}));

vi.mock('../../api/generated/runs/runs', () => ({
  createRun: generatedMocks.createRun,
}));

vi.mock('@uiw/react-codemirror', () => ({
  default: ({
    value,
    onChange,
    'aria-label': ariaLabel,
  }: {
    value: string;
    onChange: (value: string) => void;
    'aria-label'?: string;
  }) => (
    <textarea
      aria-label={ariaLabel ?? 'Advanced config YAML'}
      value={value}
      onChange={(event) => onChange(event.target.value)}
    />
  ),
}));

function successResponse(): SchemaResponse {
  return {
    ok: true,
    workflow_id: WORKFLOW_ID,
    schema: createAuthoringSchemaFixture(),
    issues: [],
  };
}

function validatedSnapshotResponse() {
  return {
    ok: true,
    workflow_id: WORKFLOW_ID,
    value: null,
    snapshot: {
      snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
      workflow_id: WORKFLOW_ID,
      schema_version: '1.1.0',
      adapter_version: '0.3.0',
      payload_digest: 'a'.repeat(64),
      validated_at: '2026-07-14T00:00:00.000Z',
      expires_at: '2026-07-14T00:30:00.000Z',
    },
    issues: [],
  };
}

function createdRunResponse() {
  return {
    ok: true,
    run: {
      run_id: 'run-snapshot-1',
      workflow_id: WORKFLOW_ID,
      inputs: {
        config: { outdir: 'results', threads: 8, qc: {} },
        samples: [
          { sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' },
        ],
        options: { strict_inputs: false },
      },
      status: 'created',
      created_at: '2026-07-14T00:01:00.000Z',
      updated_at: '2026-07-14T00:01:00.000Z',
      started_at: null,
      ended_at: null,
      current_stage: null,
      cancellation_reason: null,
      error: null,
      tags: {},
    },
    issues: [],
  };
}

function deferred<T>() {
  let resolve!: (value: T) => void;
  const promise = new Promise<T>((resolvePromise) => {
    resolve = resolvePromise;
  });
  return { promise, resolve };
}

function textFile(contents: string, name: string, type = 'text/tab-separated-values') {
  const file = new File([contents], name, { type });
  const text = vi.fn().mockResolvedValue(contents);
  Object.defineProperty(file, 'text', { value: text });
  return { file, text };
}

async function authorValidDraft(user: ReturnType<typeof userEvent.setup>) {
  await user.click(await screen.findByRole('tab', { name: 'Samples' }));
  await user.upload(
    screen.getByLabelText('Import samples TSV'),
    textFile(
      'sample\tfastq_1\tlayout\nS1\t/data/S1.fastq.gz\tSE\n',
      'samples.encode.tsv',
    ).file,
  );
  await screen.findAllByDisplayValue('S1');
  await user.click(screen.getByRole('tab', { name: 'Review' }));
  expect(screen.getByText(/structurally ready for backend validation/i)).toBeVisible();
}

describe('schema input workbench route', () => {
  beforeEach(() => {
    generatedMocks.createRun.mockReset();
    generatedMocks.createRun.mockResolvedValue(createdRunResponse());
    generatedMocks.getWorkflowSchema.mockReset();
    generatedMocks.getWorkflowSchema.mockResolvedValue(successResponse());
    generatedMocks.validateWorkflow.mockReset();
    generatedMocks.validateWorkflow.mockResolvedValue(
      validatedSnapshotResponse(),
    );
  });

  it('loads the generated schema operation and renders a draft-only workbench', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });

    expect(
      await screen.findByRole('heading', { name: 'Input workbench' }),
    ).toBeInTheDocument();
    expect(generatedMocks.getWorkflowSchema).toHaveBeenCalledWith(WORKFLOW_ID);
    expect(screen.getByText(/Draft only/i)).toBeInTheDocument();
    expect(screen.getByText(/refresh clears this draft/i)).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: /Start run/i })).not.toBeInTheDocument();
    expect(screen.getByText('Replicate analysis')).toBeVisible();
    expect(screen.getByText('ChIP-seq IDR')).toBeVisible();
    expect(
      screen.getByRole('checkbox', { name: 'Replicate analysis enabled' }),
    ).toBeChecked();
    expect(
      screen.getByRole('checkbox', { name: 'ChIP-seq IDR enabled' }),
    ).not.toBeChecked();
    expect(document.body.textContent).not.toMatch(/stage4b|stage5/);
    await userEvent.setup().click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('draft-review-json')).not.toHaveTextContent(
      /stage4b|stage5/,
    );
    expect(screen.getByRole('button', { name: 'Validate current inputs' })).toBeDisabled();
  });

  it('loads bulk authoring from the generated schema operation, not stub detail hints', async () => {
    const response = successResponse();
    response.workflow_id = 'bulk-rnaseq';
    response.schema!.schema_version = '1.0.0';
    generatedMocks.getWorkflowSchema.mockResolvedValue(response);

    renderWithRouter(appRoutes, {
      initialEntries: ['/workflows/bulk-rnaseq/new-run'],
    });

    expect(
      await screen.findByRole('heading', { name: 'Input workbench' }),
    ).toBeVisible();
    expect(generatedMocks.getWorkflowSchema).toHaveBeenCalledWith(
      'bulk-rnaseq',
    );
    await userEvent.setup().click(screen.getByRole('tab', { name: 'Options' }));
    expect(screen.getByRole('heading', { name: 'Adapter options' })).toBeVisible();
  });

  it('prefills and locks an adapter-declared constant sample cell', async () => {
    const response = successResponse();
    const items = (response.schema!.sample_schema as {
      items: {
        properties: Record<string, unknown>;
        required: string[];
      };
    }).items;
    items.properties.platform = {
      type: 'string',
      title: 'Sequencing platform',
      const: 'ILLUMINA',
    };
    items.required.push('platform');
    generatedMocks.getWorkflowSchema.mockResolvedValue(response);
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });

    await user.click(await screen.findByRole('button', { name: 'Add sample row' }));

    const platformCells = screen.getAllByLabelText('Sample 1 platform');
    expect(platformCells).toHaveLength(2);
    for (const platform of platformCells) {
      expect(platform).toHaveValue('ILLUMINA');
      expect(platform).toHaveAttribute('readonly');
      await user.type(platform, 'OTHER');
      expect(platform).toHaveValue('ILLUMINA');
    }
  });

  it('keeps validation available but disables run creation when execution is not configured', async () => {
    generatedMocks.validateWorkflow.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      value: { validated: true },
      snapshot: null,
      issues: [],
    });
    const workflowClient = createStubWorkflowClient();
    workflowClient.getWorkflow = vi.fn(async (workflowId) => {
      const response = await createStubWorkflowClient().getWorkflow(workflowId);
      if (response.workflow) {
        response.workflow = {
          ...response.workflow,
          availability: {
            ...response.workflow.availability,
            execution: 'not_configured',
            reason_code: 'WORKFLOW_EXECUTION_NOT_CONFIGURED',
          },
        };
      }
      return response;
    });
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
      clients: { workflowClient },
    });

    await authorValidDraft(user);
    const validate = screen.getByRole('button', {
      name: 'Validate current inputs',
    });
    expect(validate).toBeEnabled();
    await user.click(validate);

    expect(await screen.findByText(/execution runtime is not configured/i)).toBeVisible();
    expect(screen.getByText(/backend validation succeeded/i)).toBeVisible();
    expect(screen.queryByText(/VALIDATION_NOT_CONFIRMED/i)).not.toBeInTheDocument();
    expect(
      screen.getByRole('button', { name: 'Create run from validated inputs' }),
    ).toBeDisabled();
    expect(generatedMocks.createRun).not.toHaveBeenCalled();
  });

  it('keeps a successful deprecated-alias warning non-terminal and create-ready', async () => {
    const user = userEvent.setup();
    generatedMocks.validateWorkflow.mockResolvedValue({
      ...validatedSnapshotResponse(),
      issues: [
        {
          code: 'ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED',
          message: 'Deprecated compatibility fields should be removed.',
          severity: 'warning',
          path: 'config',
        },
      ],
    });
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);

    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));

    expect(
      await screen.findByText('ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED'),
    ).toBeVisible();
    expect(
      screen.getByText('ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED').closest('[role="status"]'),
    ).not.toBeNull();
    expect(
      screen.getByText('ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED').closest('[role="alert"]'),
    ).toBeNull();
    expect(
      screen.getByRole('button', { name: 'Create run from validated inputs' }),
    ).toBeEnabled();
  });

  it('preserves opposite legacy switches through a Form edit, Review, and validate', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const editor = screen.getByRole('textbox', { name: 'Advanced config YAML' });
    await user.clear(editor);
    await user.type(
      editor,
      'outdir: results\nthreads: 4\nstage4b: false\nstage5: true',
    );
    await user.click(screen.getByRole('button', { name: 'Form mode' }));
    fireEvent.change(screen.getByRole('spinbutton', { name: /threads/i }), {
      target: { value: '12' },
    });
    await waitFor(() =>
      expect(screen.getByTestId('draft-review-json')).toHaveTextContent(
        '"threads": 12',
      ),
    );
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const updatedEditor = screen.getByRole('textbox', {
      name: 'Advanced config YAML',
    });
    expect((updatedEditor as HTMLTextAreaElement).value).toContain('stage4b: false');
    expect((updatedEditor as HTMLTextAreaElement).value).toContain('stage5: true');
    expect((updatedEditor as HTMLTextAreaElement).value).toContain('threads: 12');
    expect((updatedEditor as HTMLTextAreaElement).value).not.toMatch(
      /replicate_analysis|chipseq_idr/,
    );

    await authorValidDraft(user);
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent('stage4b');
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent('stage5');
    expect(screen.getByTestId('draft-review-json')).not.toHaveTextContent(
      /replicate_analysis|chipseq_idr/,
    );
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    expect(generatedMocks.validateWorkflow).toHaveBeenCalledWith(
      WORKFLOW_ID,
      expect.objectContaining({
        config: {
          outdir: 'results',
          threads: 12,
          stage4b: false,
          stage5: true,
        },
      }),
    );
  });

  it('validates the exact draft, creates from only its snapshot, and navigates to the run', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);

    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    await waitFor(() => expect(generatedMocks.validateWorkflow).toHaveBeenCalledTimes(1));
    expect(generatedMocks.validateWorkflow).toHaveBeenCalledWith(
      WORKFLOW_ID,
      expect.objectContaining({
        config: expect.objectContaining({ outdir: 'results', threads: 8 }),
        samples: [
          expect.objectContaining({
            sample: 'S1',
            fastq_1: '/data/S1.fastq.gz',
            layout: 'SE',
          }),
        ],
        options: expect.objectContaining({ strict_inputs: false }),
      }),
    );
    expect(
      await screen.findByText(/This exact draft can create one run/i),
    ).toBeVisible();

    await user.click(
      screen.getByRole('button', { name: 'Create run from validated inputs' }),
    );
    await waitFor(() => expect(generatedMocks.createRun).toHaveBeenCalledTimes(1));
    expect(generatedMocks.createRun).toHaveBeenCalledWith(WORKFLOW_ID, {
      snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
    });
    await waitFor(() =>
      expect(router.state.location.pathname).toBe('/runs/run-snapshot-1'),
    );
  });

  it('invalidates the validated snapshot immediately after a semantic edit', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    expect(
      await screen.findByRole('button', {
        name: 'Create run from validated inputs',
      }),
    ).toBeEnabled();

    await user.click(screen.getByRole('tab', { name: 'Samples' }));
    await user.type(screen.getAllByLabelText('Sample 1 sample')[0], '-changed');
    await user.click(screen.getByRole('tab', { name: 'Review' }));

    expect(
      screen.getByRole('button', { name: 'Create run from validated inputs' }),
    ).toBeDisabled();
    expect(
      await screen.findByText(/Inputs changed after validation/i),
    ).toBeVisible();
  });

  it('discards a validation response when the draft changed in flight', async () => {
    const pending = deferred<ReturnType<typeof validatedSnapshotResponse>>();
    generatedMocks.validateWorkflow.mockReturnValue(pending.promise);
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));

    await user.click(screen.getByRole('tab', { name: 'Samples' }));
    await user.type(screen.getAllByLabelText('Sample 1 sample')[0], '-new');
    await act(async () => {
      pending.resolve(validatedSnapshotResponse());
      await pending.promise;
    });
    await user.click(screen.getByRole('tab', { name: 'Review' }));

    expect(
      screen.getByRole('button', { name: 'Create run from validated inputs' }),
    ).toBeDisabled();
    expect(
      screen.getByText(/Inputs changed while validation was running/i),
    ).toBeVisible();
  });

  it('does not navigate when the draft changes while create is in flight', async () => {
    const pending = deferred<ReturnType<typeof createdRunResponse>>();
    generatedMocks.createRun.mockReturnValue(pending.promise);
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    await user.click(
      await screen.findByRole('button', {
        name: 'Create run from validated inputs',
      }),
    );

    await user.click(screen.getByRole('tab', { name: 'Samples' }));
    await user.type(screen.getAllByLabelText('Sample 1 sample')[0], '-new');
    await act(async () => {
      pending.resolve(createdRunResponse());
      await pending.promise;
    });

    expect(router.state.location.pathname).toBe(
      `/workflows/${WORKFLOW_ID}/new-run`,
    );
    expect(
      await screen.findByText(/Inputs changed while run creation was running/i),
    ).toBeVisible();
  });

  it('retries an unconfirmed create with the same snapshot id', async () => {
    generatedMocks.createRun
      .mockRejectedValueOnce(new Error('private transport details'))
      .mockResolvedValueOnce(createdRunResponse());
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    const create = await screen.findByRole('button', {
      name: 'Create run from validated inputs',
    });
    await user.click(create);

    expect(await screen.findByRole('alert')).toHaveTextContent(
      /could not be confirmed/i,
    );
    expect(screen.queryByText(/private transport details/i)).not.toBeInTheDocument();
    expect(create).toBeEnabled();
    await user.click(create);
    await waitFor(() => expect(router.state.location.pathname).toBe('/runs/run-snapshot-1'));
    expect(generatedMocks.createRun).toHaveBeenNthCalledWith(1, WORKFLOW_ID, {
      snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
    });
    expect(generatedMocks.createRun).toHaveBeenNthCalledWith(2, WORKFLOW_ID, {
      snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
    });
  });

  it('requires revalidation after the API reports an expired snapshot', async () => {
    generatedMocks.createRun.mockRejectedValue(
      new ApiError(
        409,
        'VALIDATED_SNAPSHOT_EXPIRED',
        'Validated input snapshot expired.',
        [
          {
            code: 'VALIDATED_SNAPSHOT_EXPIRED',
            message: 'Validated inputs expired. Validate the current draft again.',
          },
        ],
      ),
    );
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await authorValidDraft(user);
    await user.click(screen.getByRole('button', { name: 'Validate current inputs' }));
    const create = await screen.findByRole('button', {
      name: 'Create run from validated inputs',
    });
    await user.click(create);

    expect(await screen.findByRole('alert')).toHaveTextContent(
      'VALIDATED_SNAPSHOT_EXPIRED',
    );
    expect(create).toBeDisabled();
  });

  it('keeps a stable loading region while schema loading is pending', async () => {
    const pending = deferred<SchemaResponse>();
    generatedMocks.getWorkflowSchema.mockReturnValue(pending.promise);
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });

    expect(screen.getByTestId('input-workbench-loading')).toBeInTheDocument();
    await act(async () => pending.resolve(successResponse()));
    expect(
      await screen.findByRole('heading', { name: 'Input workbench' }),
    ).toBeInTheDocument();
  });

  it('handles a null schema and a network failure with a retry', async () => {
    generatedMocks.getWorkflowSchema
      .mockRejectedValueOnce(new Error('/private/internal traceback'))
      .mockResolvedValueOnce(successResponse());
    const user = userEvent.setup();
    const first = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });

    expect(await screen.findByRole('alert')).toHaveTextContent(
      'Unable to load the workflow authoring contract.',
    );
    expect(screen.queryByText(/private|traceback/i)).not.toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: /Retry schema/i }));
    expect(
      await screen.findByRole('heading', { name: 'Input workbench' }),
    ).toBeInTheDocument();
    first.unmount();

    generatedMocks.getWorkflowSchema.mockResolvedValue({
      ok: false,
      workflow_id: WORKFLOW_ID,
      schema: null,
      issues: [],
    });
    const second = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    expect(await screen.findByRole('alert')).toHaveTextContent(
      'This workflow does not publish an authoring contract.',
    );
    second.unmount();
  });

  it('preserves unknown config fields across YAML and Form edits', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });

    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const yaml = screen.getByLabelText('Advanced config YAML');
    await user.clear(yaml);
    await user.type(
      yaml,
      'advanced_field: keep-me\noutdir: results\nthreads: 8\n',
    );
    await user.click(screen.getByRole('button', { name: 'Form mode' }));

    const threads = screen.getByRole('spinbutton', { name: /threads/i });
    fireEvent.change(threads, { target: { value: '12' } });
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    expect(
      (screen.getByLabelText('Advanced config YAML') as HTMLTextAreaElement).value,
    ).toContain('advanced_field: keep-me');
  });

  it('shows adapter schema errors next to form fields after blur', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });

    const threads = screen.getByRole('spinbutton', { name: /threads/i });
    fireEvent.change(threads, { target: { value: '0' } });
    fireEvent.blur(threads);

    expect(await screen.findByText(/must be >= 1/i)).toBeInTheDocument();
  });

  it('rejects an unsafe integer emitted by the RJSF Form without previewing stale config', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });

    const threads = screen.getByRole('spinbutton', { name: /threads/i });
    fireEvent.change(threads, {
      target: { value: String(Number.MAX_SAFE_INTEGER) },
    });
    expect(screen.queryByTestId('config-form-safety-issue')).not.toBeInTheDocument();
    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent(
      String(Number.MAX_SAFE_INTEGER),
    );
    await user.click(screen.getByRole('tab', { name: 'Config' }));

    fireEvent.change(screen.getByRole('spinbutton', { name: /threads/i }), {
      target: { value: String(Number.MAX_SAFE_INTEGER + 1) },
    });

    expect(await screen.findByTestId('config-form-safety-issue')).toHaveTextContent(
      /cannot be represented safely as JSON/i,
    );
    expect(screen.getByRole('spinbutton', { name: /threads/i })).toHaveValue(
      Number.MAX_SAFE_INTEGER,
    );
    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent(
      /preview is unavailable/i,
    );
    expect(screen.getByTestId('draft-review-json')).not.toHaveTextContent(
      String(Number.MAX_SAFE_INTEGER + 1),
    );
    await user.click(screen.getByRole('tab', { name: 'Config' }));
    await user.click(
      screen.getByRole('button', { name: 'Use previous safe config' }),
    );
    expect(screen.queryByTestId('config-form-safety-issue')).not.toBeInTheDocument();
    expect(screen.getByRole('spinbutton', { name: /threads/i })).toHaveValue(
      Number.MAX_SAFE_INTEGER,
    );
  });

  it('resets unsafe options Form state to the previous canonical value', async () => {
    const response = successResponse();
    response.schema!.option_schema = {
      $schema: 'https://json-schema.org/draft/2020-12/schema',
      type: 'object',
      properties: {
        safe_limit: {
          type: 'integer',
          default: Number.MAX_SAFE_INTEGER,
        },
      },
      additionalProperties: false,
    };
    generatedMocks.getWorkflowSchema.mockResolvedValue(response);
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=options`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });

    fireEvent.change(screen.getByRole('spinbutton', { name: /safe_limit/i }), {
      target: { value: String(Number.MAX_SAFE_INTEGER + 1) },
    });
    expect(await screen.findByRole('alert')).toHaveTextContent(
      /Options form data cannot be represented safely as JSON/i,
    );
    expect(screen.getByRole('spinbutton', { name: /safe_limit/i })).toHaveValue(
      Number.MAX_SAFE_INTEGER,
    );
    await user.click(
      screen.getByRole('button', { name: 'Use previous safe options' }),
    );
    expect(screen.queryByRole('alert')).not.toBeInTheDocument();
    expect(screen.getByRole('spinbutton', { name: /safe_limit/i })).toHaveValue(
      Number.MAX_SAFE_INTEGER,
    );
  });

  it('retains invalid YAML, blocks stale Review, and recovers in place', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const yaml = screen.getByLabelText('Advanced config YAML');
    fireEvent.change(yaml, { target: { value: 'threads: [' } });

    expect(screen.getByRole('alert')).toHaveTextContent(/Config YAML/i);
    expect(screen.getByRole('button', { name: 'Form mode' })).toBeDisabled();
    expect(screen.getByRole('tab', { name: 'Review' })).toBeDisabled();

    fireEvent.change(yaml, { target: { value: 'threads: 4\n' } });
    expect(screen.queryByRole('alert')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Form mode' })).toBeEnabled();
  });

  it('never exposes the last-known-good preview when history activates Review with invalid YAML', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const yaml = screen.getByLabelText('Advanced config YAML');
    fireEvent.change(yaml, {
      target: { value: 'threads: 4\nlast_known_good: must-not-leak\n' },
    });
    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent(
      'last_known_good',
    );

    await act(async () => router.navigate(-1));
    fireEvent.change(screen.getByLabelText('Advanced config YAML'), {
      target: { value: 'threads: [' },
    });
    await act(async () => router.navigate(1));

    expect(screen.getByRole('tab', { name: 'Review' })).toHaveAttribute(
      'data-state',
      'active',
    );
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent(
      /preview is unavailable/i,
    );
    expect(screen.getByTestId('draft-review-json')).not.toHaveTextContent(
      'last_known_good',
    );
  });

  it('imports and edits quoted CRLF TSV rows in adapter column order', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });

    const { file } = textFile(
      'layout\tsample\tfastq_1\r\n' +
        'SE\t"S1"\t/data/S1.fastq.gz\r\n',
      'samples.encode.tsv',
    );
    await user.upload(screen.getByLabelText('Import samples TSV'), file);

    expect((await screen.findAllByLabelText('Sample 1 sample'))[0]).toHaveValue('S1');
    expect(screen.getAllByLabelText('Sample 1 fastq_1')[0]).toHaveValue(
      '/data/S1.fastq.gz',
    );
    await user.clear(screen.getAllByLabelText('Sample 1 sample')[0]);
    await user.type(screen.getAllByLabelText('Sample 1 sample')[0], 'edited-S1');
    await user.click(screen.getByRole('button', { name: 'Add sample row' }));
    expect(screen.getAllByLabelText('Sample 2 sample')[0]).toHaveValue('');
    await user.click(screen.getAllByRole('button', { name: 'Remove sample 2' })[0]);

    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent('edited-S1');
    const reviewText = screen.getByTestId('draft-review-json').textContent ?? '';
    expect(reviewText.indexOf('sample')).toBeLessThan(reviewText.indexOf('fastq_1'));
  });

  it('blocks Review for an overlong manually edited sample cell without echoing its value', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'Add sample row' }));
    const overlong = 'x'.repeat(4_097);
    fireEvent.change(screen.getAllByLabelText('Sample 1 sample')[0], {
      target: { value: overlong },
    });

    const issue = await screen.findByTestId('sample-transport-issue');
    expect(issue).toHaveTextContent(/exceeds the workflow authoring limit/i);
    expect(issue).toHaveTextContent(/Row 1.*Column sample/i);
    expect(issue).not.toHaveTextContent(overlong);
    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByText(/needs attention before backend validation/i)).toBeVisible();
    expect(screen.getByTestId('draft-review-json')).toHaveTextContent(
      /preview is unavailable/i,
    );
  });

  it('describes multiple invalid sample cells without assigning the first cell coordinates to every control', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'Add sample row' }));
    const sample = screen.getAllByLabelText('Sample 1 sample')[0];
    const fastq = screen.getAllByLabelText('Sample 1 fastq_1')[0];

    await act(async () => {
      fireEvent.change(sample, { target: { value: 'x'.repeat(4_097) } });
      fireEvent.change(fastq, { target: { value: 'y'.repeat(4_097) } });
    });

    expect(sample).toHaveAttribute('aria-invalid', 'true');
    expect(fastq).toHaveAttribute('aria-invalid', 'true');
    expect(sample).toHaveAttribute(
      'aria-describedby',
      'sample-cell-invalid-description',
    );
    expect(fastq).toHaveAttribute(
      'aria-describedby',
      'sample-cell-invalid-description',
    );
    const description = document.getElementById(
      'sample-cell-invalid-description',
    );
    expect(description).toHaveTextContent(/cannot be transported/i);
    expect(description).not.toHaveTextContent(/Row 1|Column sample/i);
    const issue = screen.getByTestId('sample-transport-issue');
    expect(issue).toHaveTextContent(/Row 1.*Column sample/i);
    expect(issue).not.toHaveTextContent(/fastq_1/i);
  });

  it('keeps a manual cell controlled and focused through repeated invalid edits', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'Add sample row' }));
    const sample = screen.getAllByLabelText('Sample 1 sample')[0];
    const overlong = 'x'.repeat(4_097);
    await user.click(sample);
    act(() => fireEvent.change(sample, { target: { value: overlong } }));
    expect(sample).toHaveFocus();

    (sample as HTMLInputElement).setSelectionRange(
      overlong.length,
      overlong.length,
    );
    await user.keyboard('yz');
    expect(sample).toHaveValue(`${overlong}yz`);
    expect(sample).toHaveFocus();
    expect(screen.getByTestId('sample-transport-issue')).toBeInTheDocument();

    await user.clear(sample);
    expect(sample).toHaveFocus();
    await user.keyboard('recovered');
    expect(sample).toHaveValue('recovered');
    expect(sample).toHaveFocus();
    expect(screen.queryByTestId('sample-transport-issue')).not.toBeInTheDocument();
  });

  it('wraps a maximum-length adapter column in the Review issue', async () => {
    const longColumn = 'c'.repeat(128);
    const response = successResponse();
    const items = (response.schema?.sample_schema as {
      items: { properties: Record<string, unknown>; required: string[] };
    }).items;
    items.properties = {
      [longColumn]: { type: 'string', maxLength: 4_096 },
    };
    items.required = [longColumn];
    generatedMocks.getWorkflowSchema.mockResolvedValue(response);
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'Add sample row' }));
    fireEvent.change(screen.getAllByLabelText(`Sample 1 ${longColumn}`)[0], {
      target: { value: 'x'.repeat(4_097) },
    });
    await user.click(screen.getByRole('tab', { name: 'Review' }));

    const column = screen.getByText(longColumn, { selector: 'code' });
    expect(column).toHaveClass('break-all');
    expect(column.closest('li')).toHaveClass('min-w-0', 'break-words');
  });

  it('keeps imported rows after a malformed replacement file', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    const input = screen.getByLabelText('Import samples TSV');
    await user.upload(
      input,
      textFile(
        'sample\tfastq_1\tlayout\nS1\t/data/S1.fastq.gz\tSE\n',
        'valid.tsv',
      ).file,
    );
    expect((await screen.findAllByLabelText('Sample 1 sample'))[0]).toHaveValue('S1');

    await user.upload(
      input,
      textFile('sample\tunknown\nS2\tx\n', 'bad.tsv').file,
    );
    expect(await screen.findByRole('alert')).toHaveTextContent(
      'The sample file contains a column not declared by this workflow. Row 1. Column unknown.',
    );
    expect(screen.getAllByLabelText('Sample 1 sample')[0]).toHaveValue('S1');

    await user.upload(
      input,
      textFile(
        `sample\tfastq_1\tlayout\nS2\t${'x'.repeat(4_097)}\tSE\n`,
        'overlong.tsv',
      ).file,
    );
    expect(screen.getByRole('alert')).toHaveTextContent(
      'A sample cell exceeds the workflow authoring limit. Row 2. Column fastq_1.',
    );
    expect(screen.getAllByLabelText('Sample 1 sample')[0]).toHaveValue('S1');
  });

  it('keeps the newest TSV selection when an older file read finishes later', async () => {
    const firstRead = deferred<string>();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    const input = screen.getByLabelText('Import samples TSV');
    const firstFile = new File(['pending'], 'first.tsv', {
      type: 'text/tab-separated-values',
    });
    Object.defineProperty(firstFile, 'text', {
      value: vi.fn().mockReturnValue(firstRead.promise),
    });
    const secondFile = textFile(
      'sample\tfastq_1\tlayout\nnewest\t/data/newest.fastq.gz\tSE\n',
      'second.tsv',
    ).file;

    fireEvent.change(input, { target: { files: [firstFile] } });
    fireEvent.change(input, { target: { files: [secondFile] } });
    expect((await screen.findAllByLabelText('Sample 1 sample'))[0]).toHaveValue(
      'newest',
    );

    await act(async () => {
      firstRead.resolve(
        'sample\tfastq_1\tlayout\nolder\t/data/older.fastq.gz\tSE\n',
      );
      await firstRead.promise;
    });
    expect(screen.getAllByLabelText('Sample 1 sample')[0]).toHaveValue('newest');
  });

  it('drives steps from the URL and keeps the in-memory draft through history navigation', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'Add sample row' }));
    await user.type(screen.getAllByLabelText('Sample 1 sample')[0], 'history-S1');
    await user.click(screen.getByRole('tab', { name: 'Options' }));
    expect(router.state.location.search).toContain('step=options');

    await act(async () => router.navigate(-1));
    expect(await screen.findByRole('tab', { name: 'Samples' })).toHaveAttribute(
      'data-state',
      'active',
    );
    expect(screen.getAllByLabelText('Sample 1 sample')[0]).toHaveValue('history-S1');
  });

  it('rejects an oversized file before reading it', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    const file = new File([new Uint8Array(2_097_153)], 'large.tsv', {
      type: 'text/tab-separated-values',
    });
    const textSpy = vi.fn().mockResolvedValue('should not be read');
    Object.defineProperty(file, 'text', { value: textSpy });
    fireEvent.change(screen.getByLabelText('Import samples TSV'), {
      target: { files: [file] },
    });
    expect(await screen.findByRole('alert')).toHaveTextContent(/larger/i);
    expect(textSpy).not.toHaveBeenCalled();
  });
});
