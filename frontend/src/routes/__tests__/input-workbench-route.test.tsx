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
