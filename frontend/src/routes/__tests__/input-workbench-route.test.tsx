import { act, fireEvent, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import type { SchemaResponse } from '../../api/generated/models';
import { appRoutes } from '../../app/router';
import { createAuthoringSchemaFixture, WORKFLOW_ID } from '../../features/input-workbench/test-fixtures';
import { renderWithRouter } from '../../test/test-utils';

const generatedMocks = vi.hoisted(() => ({
  getWorkflowSchema: vi.fn(),
}));

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: generatedMocks.getWorkflowSchema,
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

describe('schema input workbench route', () => {
  beforeEach(() => {
    generatedMocks.getWorkflowSchema.mockReset();
    generatedMocks.getWorkflowSchema.mockResolvedValue(successResponse());
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
    expect(screen.queryByRole('button', { name: /Create run|Start run|Validate/i })).not.toBeInTheDocument();
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
    expect(await screen.findByRole('alert')).toHaveTextContent(/not declared/i);
    expect(screen.getAllByLabelText('Sample 1 sample')[0]).toHaveValue('S1');
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
    const schema = createAuthoringSchemaFixture();
    schema.limits.max_request_bytes = 4;
    generatedMocks.getWorkflowSchema.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      schema,
      issues: [],
    });
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run?step=samples`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    const { file, text: textSpy } = textFile('sample\nS1\n', 'large.tsv');
    fireEvent.change(screen.getByLabelText('Import samples TSV'), {
      target: { files: [file] },
    });
    expect(await screen.findByRole('alert')).toHaveTextContent(/larger/i);
    expect(textSpy).not.toHaveBeenCalled();
  });
});
