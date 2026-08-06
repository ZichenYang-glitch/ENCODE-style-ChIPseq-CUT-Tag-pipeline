import { beforeEach, describe, expect, it, vi } from 'vitest';
import { fireEvent, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { appRoutes } from './app/router';
import { renderWithRouter } from './test/test-utils';
import { createAuthoringSchemaFixture } from './features/input-workbench/test-fixtures';

const generatedMocks = vi.hoisted(() => ({
  createRun: vi.fn(),
  getWorkflowSchema: vi.fn(),
  validateWorkflow: vi.fn(),
}));
const referenceProfileMocks = vi.hoisted(() => ({
  listCompatibleReferenceProfiles: vi.fn(),
}));

vi.mock('./api/generated/workflows/workflows', () => ({
  getWorkflowSchema: generatedMocks.getWorkflowSchema,
  listCompatibleReferenceProfiles:
    referenceProfileMocks.listCompatibleReferenceProfiles,
  validateWorkflow: generatedMocks.validateWorkflow,
}));
vi.mock('./api/generated/runs/runs', () => ({
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

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';
const GRCH38_PROFILE = {
  profile_id: `refp_${'1'.repeat(32)}`,
  revision_id: `refpr_${'1'.repeat(32)}`,
  revision_number: 1,
  display_name: 'Human GRCh38',
  organism: 'Homo sapiens',
  assembly: 'GRCh38',
  identity_sha256: '1'.repeat(64),
};
const MM10_PROFILE = {
  profile_id: `refp_${'2'.repeat(32)}`,
  revision_id: `refpr_${'2'.repeat(32)}`,
  revision_number: 1,
  display_name: 'Mouse mm10',
  organism: 'Mus musculus',
  assembly: 'mm10',
  identity_sha256: '2'.repeat(64),
};

function validationResponse(referenceProfile = GRCH38_PROFILE) {
  return {
    ok: true,
    workflow_id: WORKFLOW_ID,
    value: null,
    snapshot: {
      snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
      workflow_id: WORKFLOW_ID,
      schema_version: '1.0.0',
      adapter_version: '0.3.0',
      payload_digest: 'a'.repeat(64),
      validated_at: '2026-08-03T00:00:00.000Z',
      expires_at: '2026-08-03T00:30:00.000Z',
      reference_profile: referenceProfile,
    },
    issues: [],
  };
}

async function authorValidDraft(user: ReturnType<typeof userEvent.setup>) {
  await user.click(await screen.findByRole('tab', { name: 'Samples' }));
  const contents =
    'sample\tfastq_1\tlayout\nS1\t/data/S1.fastq.gz\tSE\n';
  const file = new File([contents], 'samples.tsv', {
    type: 'text/tab-separated-values',
  });
  Object.defineProperty(file, 'text', {
    value: vi.fn().mockResolvedValue(contents),
  });
  await user.upload(screen.getByLabelText('Import samples TSV'), file);
  await screen.findAllByDisplayValue('S1');
  await user.click(screen.getByRole('tab', { name: 'Review' }));
}

describe('App shell', () => {
  beforeEach(() => {
    generatedMocks.getWorkflowSchema.mockReset();
    generatedMocks.getWorkflowSchema.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      schema: createAuthoringSchemaFixture(),
      issues: [],
    });
    generatedMocks.validateWorkflow.mockReset();
    generatedMocks.validateWorkflow.mockResolvedValue(validationResponse());
    generatedMocks.createRun.mockReset();
    referenceProfileMocks.listCompatibleReferenceProfiles.mockReset();
    referenceProfileMocks.listCompatibleReferenceProfiles.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      profiles: [GRCH38_PROFILE],
      issues: [],
    });
  });

  it('renders the HelixWeave heading', async () => {
    renderWithRouter(appRoutes, { initialEntries: ['/workflows'] });
    expect(
      await screen.findByRole('heading', { name: 'HelixWeave' }),
    ).toBeInTheDocument();
  });

  it('shows the stub workflow catalog', async () => {
    renderWithRouter(appRoutes, { initialEntries: ['/workflows'] });
    expect(await screen.findByText(/ENCODE-style ChIP-seq/i)).toBeInTheDocument();
  });

  it('loads schema hints when a workflow is selected', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, { initialEntries: ['/workflows'] });
    await user.click(await screen.findByText(/ENCODE-style ChIP-seq/i));
    expect(await screen.findByText(/Config schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Sample schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Options schema/i)).toBeInTheDocument();
  });

  it('offers the reference-bound Input Workbench from workflow detail', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });

    expect(await screen.findByRole('link', { name: 'Author inputs' })).toBeVisible();
    expect(screen.getByTestId('validate-button')).toBeInTheDocument();
    await user.click(screen.getByRole('link', { name: 'Author inputs' }));

    await waitFor(() =>
      expect(router.state.location.pathname).toBe(
        `/workflows/${WORKFLOW_ID}/new-run`,
      ),
    );
    expect(await screen.findByRole('heading', { name: 'Input workbench' })).toBeVisible();
    expect(screen.getByRole('combobox', { name: 'Reference profile' })).toHaveValue(
      GRCH38_PROFILE.revision_id,
    );
  });

  it('renders backend validation issues in the reference-bound workbench', async () => {
    generatedMocks.validateWorkflow.mockResolvedValue({
      ok: false,
      workflow_id: WORKFLOW_ID,
      value: null,
      snapshot: null,
      issues: [
        {
          code: 'ENCODE_SAMPLES_INVALID',
          message: 'Sample sheet is invalid.',
          severity: 'error',
          source: 'adapter',
        },
      ],
    });
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await authorValidDraft(user);
    await user.click(screen.getByTestId('validate-draft-button'));

    expect(await screen.findByText(/Sample sheet is invalid/i)).toBeVisible();
    expect(screen.getByText('ENCODE_SAMPLES_INVALID')).toBeVisible();
  });

  it('renders a local YAML error without calling the backend', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await screen.findByRole('heading', { name: 'Input workbench' });
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const editor = screen.getByLabelText('Advanced config YAML');
    fireEvent.change(editor, { target: { value: 'threads: [' } });

    expect(await screen.findByRole('alert')).toHaveTextContent(/Config YAML/i);
    expect(screen.getByRole('tab', { name: 'Review' })).toBeDisabled();
    expect(generatedMocks.validateWorkflow).not.toHaveBeenCalled();
  });

  it('keeps the existing read-only validation agent on workflow detail', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });

    expect(await screen.findByRole('link', { name: 'Author inputs' })).toBeVisible();
    expect(screen.getByText(/Validation Assistant/i)).toBeInTheDocument();
  });

  it('keeps existing run controls on workflow detail', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });

    expect(await screen.findByRole('link', { name: 'Author inputs' })).toBeVisible();
    expect(screen.getByText(/Run progress/i)).toBeInTheDocument();
    expect(screen.getByTestId('create-run-button')).toBeDisabled();
  });

  it('enables run creation only after validation confirms the selected reference revision', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    await authorValidDraft(user);
    await user.click(screen.getByTestId('validate-draft-button'));

    expect(
      await screen.findByTestId('create-validated-run-button'),
    ).toBeEnabled();
    expect(generatedMocks.validateWorkflow).toHaveBeenCalledWith(
      WORKFLOW_ID,
      expect.objectContaining({
        reference_profile_revision_id: GRCH38_PROFILE.revision_id,
      }),
    );
  });

  it('invalidates a confirmed snapshot when the user switches reference revisions', async () => {
    referenceProfileMocks.listCompatibleReferenceProfiles.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      profiles: [GRCH38_PROFILE, MM10_PROFILE],
      issues: [],
    });
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
    });
    const reference = await screen.findByRole('combobox', {
      name: 'Reference profile',
    });
    await user.selectOptions(reference, GRCH38_PROFILE.revision_id);
    await authorValidDraft(user);
    await user.click(screen.getByTestId('validate-draft-button'));
    expect(
      await screen.findByTestId('create-validated-run-button'),
    ).toBeEnabled();

    await user.selectOptions(reference, MM10_PROFILE.revision_id);

    expect(screen.getByTestId('create-validated-run-button')).toBeDisabled();
    expect(await screen.findByText(/Inputs changed after validation/i)).toBeVisible();
  });
});
