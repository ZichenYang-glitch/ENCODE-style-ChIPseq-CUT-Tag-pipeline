import { beforeEach, describe, expect, it, vi } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { RouteObject } from 'react-router-dom';
import { renderWithRouter } from '../../test/test-utils';
import { fetchSessionState } from '../../api/authClient';
import { authoringListProjects } from '../../api/generated/authoring/authoring';
import {
  catalogCreateExperiment,
  catalogGenerateClassificationSuggestion,
  catalogCreateTerm,
  catalogListExperimentRevisionRunResults,
  catalogListExperimentRevisions,
  catalogListExperiments,
  catalogListSampleRevisionRunResults,
  catalogListSamples,
  catalogListTerms,
  catalogReviseExperiment,
} from '../../api/generated/catalog/catalog';
import type {
  CatalogExperimentResponse,
  CatalogPublishedRunResponse,
  CatalogTermSummaryResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { ExperimentCatalogPage } from './ExperimentCatalogPage';

vi.mock('../../api/generated/catalog/catalog', () => ({
  catalogListExperiments: vi.fn(),
  catalogListExperimentRevisions: vi.fn(),
  catalogListExperimentRevisionRunResults: vi.fn(),
  catalogListSampleRevisionRunResults: vi.fn(),
  catalogListTerms: vi.fn(),
  catalogListSamples: vi.fn(),
  catalogCreateTerm: vi.fn(),
  catalogGenerateClassificationSuggestion: vi.fn(),
  catalogPreviewProposal: vi.fn(),
  catalogAdoptProposal: vi.fn(),
  catalogCreateExperiment: vi.fn(),
  catalogReviseExperiment: vi.fn(),
}));

vi.mock('../../api/generated/authoring/authoring', () => ({
  authoringListProjects: vi.fn(),
}));

const PROJECT_ID = `prj_${'1'.repeat(32)}`;
const ASSAY_ID = `term_${'a'.repeat(32)}`;
const ASSAY_REVISION = `trmv_${'a'.repeat(32)}`;
const CELL_ID = `term_${'c'.repeat(32)}`;
const CELL_REVISION = `trmv_${'c'.repeat(32)}`;
const TARGET_ID = `term_${'d'.repeat(32)}`;
const TARGET_REVISION = `trmv_${'d'.repeat(32)}`;
const ROLE_ID = `term_${'e'.repeat(32)}`;
const ROLE_REVISION = `trmv_${'e'.repeat(32)}`;
const SAMPLE_ID = `smp_${'2'.repeat(32)}`;
const SAMPLE_REVISION = `smpr_${'3'.repeat(32)}`;
const EXPERIMENT_ID = `exp_${'4'.repeat(32)}`;
const OTHER_EXPERIMENT_ID = `exp_${'7'.repeat(32)}`;
const REV_1 = `expr_${'5'.repeat(32)}`;
const REV_2 = `expr_${'6'.repeat(32)}`;

const routes: RouteObject[] = [
  { path: '/catalog', element: <ExperimentCatalogPage /> },
];

function term(
  termType: string,
  displayName: string,
  termId: string,
  termRevisionId: string,
): CatalogTermSummaryResponse {
  return {
    term_id: termId,
    term_type: termType,
    term_revision_id: termRevisionId,
    revision_number: 1,
    display_name: displayName,
    archived_at: null,
  };
}

function experiment(input: {
  experimentId: string;
  revisionId: string;
  revisionNumber: number;
  displayName: string;
}): CatalogExperimentResponse {
  return {
    experiment_id: input.experimentId,
    experiment_revision_id: input.revisionId,
    revision_number: input.revisionNumber,
    display_name: input.displayName,
    assay: {
      status: 'known',
      term_id: ASSAY_ID,
      term_revision_id: ASSAY_REVISION,
      display_name: 'CUT&Tag',
    },
    target: {
      status: 'known',
      term_id: TARGET_ID,
      term_revision_id: TARGET_REVISION,
      display_name: 'H3K27me3',
    },
    batch: { status: 'unknown' },
    samples: [
      {
        sample_id: SAMPLE_ID,
        sample_revision_id: SAMPLE_REVISION,
        ordinal: 0,
        replicate_role: {
          status: 'known',
          term_id: ROLE_ID,
          term_revision_id: ROLE_REVISION,
          display_name: 'treatment',
        },
        species: { status: 'unknown' },
        tissue: { status: 'unknown' },
        cell_type: {
          status: 'known',
          values: [
            {
              term_id: CELL_ID,
              term_revision_id: CELL_REVISION,
              display_name: 'regulatory T cell',
            },
          ],
        },
        cell_line: { status: 'not_applicable', values: [] },
        treatment: { status: 'not_applicable', values: [] },
        tag: { status: 'not_applicable', values: [] },
      },
    ],
  };
}

const currentRevision = experiment({
  experimentId: EXPERIMENT_ID,
  revisionId: REV_2,
  revisionNumber: 2,
  displayName: 'Treg CUT&Tag v2',
});
const olderRevision = experiment({
  experimentId: EXPERIMENT_ID,
  revisionId: REV_1,
  revisionNumber: 1,
  displayName: 'Treg CUT&Tag v1',
});
const pageTwo = experiment({
  experimentId: OTHER_EXPERIMENT_ID,
  revisionId: `expr_${'8'.repeat(32)}`,
  revisionNumber: 1,
  displayName: 'CUT&Tag page 20',
});

const publishedRun: CatalogPublishedRunResponse = {
  run_id: 'run_treg_h3k27me3_01',
  provenance: 'resolved',
  publication_evidence: 'present',
  artifact_association_kind: 'run_input_scope',
  workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
  artifact_id: 'peaks',
  output_type: 'peaks',
  artifact_generation: `artifactgen-${'a'.repeat(64)}`,
  artifact_revision: `artifactrev-${'b'.repeat(64)}`,
  published_at: '2026-11-18T00:00:00Z',
  current_artifact_generation: `artifactgen-${'a'.repeat(64)}`,
  generation_status: 'current',
  run_sample_binding: {
    binding_mode: 'bound_v1',
    provenance: 'resolved',
    associated_run_samples: [
      {
        sample_id: SAMPLE_ID,
        sample_revision_id: SAMPLE_REVISION,
        revision_number: 2,
        ordinal: 0,
      },
    ],
  },
};

function installCatalogReads() {
  vi.mocked(catalogListExperiments).mockClear();
  vi.mocked(catalogListSamples).mockClear();
  vi.mocked(catalogCreateExperiment).mockClear();
  vi.mocked(catalogCreateTerm).mockClear();
  vi.mocked(catalogReviseExperiment).mockClear();
  vi.mocked(authoringListProjects).mockResolvedValue({
    ok: true,
    projects: [{ project_id: PROJECT_ID, display_name: 'Treg project' }],
    next_cursor: null,
  });
  vi.mocked(catalogListTerms).mockResolvedValue({
    terms: [
      term('assay', 'CUT&Tag', ASSAY_ID, ASSAY_REVISION),
      term('assay', 'ATAC-seq', `term_${'b'.repeat(32)}`, `trmv_${'b'.repeat(32)}`),
      term('cell_type', 'regulatory T cell', CELL_ID, CELL_REVISION),
      term('target', 'H3K27me3', TARGET_ID, TARGET_REVISION),
      term('batch', 'batch-1', `term_${'9'.repeat(32)}`, `trmv_${'9'.repeat(32)}`),
      term('replicate_role', 'treatment', ROLE_ID, ROLE_REVISION),
    ],
    next_cursor: null,
  });
  vi.mocked(catalogListSamples).mockResolvedValue({
    project_id: PROJECT_ID,
    samples: [
      {
        sample_id: SAMPLE_ID,
        sample_revision_id: SAMPLE_REVISION,
        revision_number: 2,
        display_name: 'treg',
        created_at: '2026-11-18T00:00:00Z',
        catalog_state: 'absent',
      },
    ],
    next_cursor: null,
  });
  vi.mocked(catalogListExperiments).mockImplementation(async (_projectId, params) => {
    if (params?.cursor === 'catcur_page2') {
      return {
        ok: true,
        project_id: PROJECT_ID,
        revision_view: 'linked',
        experiments: [pageTwo],
        next_cursor: null,
      };
    }
    return {
      ok: true,
      project_id: PROJECT_ID,
      revision_view: 'linked',
      experiments: [currentRevision],
      next_cursor: 'catcur_page2',
    };
  });
  vi.mocked(catalogListExperimentRevisions).mockResolvedValue({
    experiment_revisions: [olderRevision, currentRevision],
    next_cursor: null,
  });
  vi.mocked(catalogListExperimentRevisionRunResults).mockImplementation(
    async (_projectId, experimentRevisionId) => ({
      ok: true,
      project_id: PROJECT_ID,
      experiment_revision_id: experimentRevisionId,
      trace_status: 'unresolved',
      reason: 'experiment_snapshot_binding_absent',
      results: [],
      next_cursor: null,
    }),
  );
  vi.mocked(catalogListSampleRevisionRunResults).mockResolvedValue({
    ok: true,
    project_id: PROJECT_ID,
    sample_revision_id: SAMPLE_REVISION,
    results: [publishedRun],
    next_cursor: null,
  });
  vi.mocked(catalogReviseExperiment).mockResolvedValue({
    ...currentRevision,
    experiment_revision_id: `expr_${'f'.repeat(32)}`,
    revision_number: 3,
    display_name: 'Treg CUT&Tag v3',
  });
}

async function search(user: ReturnType<typeof userEvent.setup>) {
  await waitFor(() => {
    expect(screen.getByLabelText('Project')).toHaveValue(PROJECT_ID);
  });
  await user.selectOptions(screen.getByLabelText('Assay'), ASSAY_ID);
  await user.selectOptions(screen.getByLabelText('Cell type'), CELL_ID);
  await user.click(screen.getByRole('button', { name: 'Find experiments' }));
}

describe('Experiment catalog page', () => {
  beforeEach(() => {
    vi.mocked(fetchSessionState).mockResolvedValue({
      ok: true,
      setup_required: false,
      authenticated: true,
      principal: {
        user_id: `usr_${'0'.repeat(32)}`,
        username: 'test-admin',
        role: 'administrator',
      },
    });
    installCatalogReads();
  });

  it('sends one authorized target clue for the selected revision and reuses its request key after a lost response', async () => {
    const user = userEvent.setup();
    vi.mocked(catalogGenerateClassificationSuggestion).mockReset()
      .mockRejectedValueOnce(new Error('lost response'))
      .mockResolvedValueOnce({
        ok: true, outcome: 'insufficient_evidence', field: 'target', reason_code: 'TARGET_TERM_MATCH',
        experiment_id: EXPERIMENT_ID, experiment_revision_id: REV_2, experiment_revision_number: 2,
        catalog_revision_changed: false, decision: 'insufficient_evidence', proposal: null,
        projection: { field: 'experiment.target', reason_code: 'TARGET_TERM_MATCH', clue: 'CTCF',
          candidates: [{ token: 'ctcf', label: 'CTCF' }] },
      });
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));
    const ask = await screen.findByRole('button', { name: 'Ask Peri' });
    expect(ask).toBeDisabled();
    await user.type(screen.getByRole('textbox', { name: 'Target clue' }), 'CTCF');
    await user.selectOptions(screen.getByLabelText('Candidate target'), TARGET_REVISION);
    await user.click(ask);
    expect(await screen.findByText(/could not retrieve the suggestion/i)).toBeInTheDocument();
    await user.click(ask);
    expect(await screen.findByText(/insufficient evidence to suggest a target/i)).toBeInTheDocument();
    const [first, retry] = vi.mocked(catalogGenerateClassificationSuggestion).mock.calls;
    expect(first).toEqual([PROJECT_ID, EXPERIMENT_ID, {
      experiment_revision_id: REV_2, field: 'target', reason_code: 'TARGET_TERM_MATCH',
      clue: 'CTCF', vocabulary_term_revision_id: TARGET_REVISION,
    }, { 'Idempotency-Key': expect.any(String) }]);
    expect(retry).toEqual(first);
    expect(screen.queryByRole('button', { name: 'Preview selected' })).toBeNull();
    expect(catalogReviseExperiment).not.toHaveBeenCalled();
  });

  it('searches, pages, opens an older revision, and shows unbound run scope', async () => {
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });

    expect(
      await screen.findByRole('heading', { name: 'Experiment catalog' }),
    ).toBeInTheDocument();
    expect(screen.getByText(/record and revise experiment classification/i)).toBeInTheDocument();
    await search(user);

    expect(catalogListExperiments).toHaveBeenCalledWith(PROJECT_ID, {
      limit: 20,
      assay_term_id: [ASSAY_ID],
      cell_type_term_id: [CELL_ID],
    });
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));
    expect(
      await screen.findByRole('heading', { name: 'Treg CUT&Tag v2' }),
    ).toBeInTheDocument();
    const trace = screen.getByRole('region', { name: 'Experiment trace' });
    expect(
      await within(trace).findByRole('heading', { name: '未绑定' }),
    ).toBeInTheDocument();
    expect(trace).toHaveTextContent('experiment_snapshot_binding_absent');
    expect(trace).toHaveTextContent('也不把运行样本写成文件贡献者');

    const scope = screen.getByRole('region', { name: 'Run input scope' });
    expect(scope).toHaveTextContent('运行输入范围');
    expect(scope).toHaveTextContent('不是文件贡献者');
    expect(scope).toHaveTextContent('run_input_scope');
    expect(scope).toHaveTextContent(`运行样本 ${SAMPLE_ID}`);
    expect(within(scope).queryByRole('heading', { name: '文件贡献者' })).not.toBeInTheDocument();

    await user.selectOptions(screen.getByLabelText('Revision'), REV_1);
    expect(
      await screen.findByRole('heading', { name: 'Treg CUT&Tag v1' }),
    ).toBeInTheDocument();
    expect(screen.getByRole('region', { name: 'Experiment revision' })).toHaveTextContent(
      '旧版 1',
    );

    await user.click(screen.getByRole('button', { name: 'Next page' }));
    expect(
      await screen.findByRole('button', { name: /CUT&Tag page 20/ }),
    ).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: /ATAC excluded/ })).not.toBeInTheDocument();
    expect(catalogListExperiments).toHaveBeenLastCalledWith(
      PROJECT_ID,
      expect.objectContaining({ cursor: 'catcur_page2', limit: 20 }),
    );
  });

  it('lets an administrator revise the current head without rewriting the old revision', async () => {
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));
    const revise = await screen.findByRole('form', {
      name: 'Revise experiment classification',
    });
    await waitFor(() => {
      expect(
        within(revise).getByRole('button', { name: 'Revise classification' }),
      ).toBeEnabled();
    });
    const name = within(revise).getByLabelText('Display name');
    await user.clear(name);
    await user.type(name, 'Treg CUT&Tag v3');
    await user.click(within(revise).getByRole('button', { name: 'Revise classification' }));
    expect(
      await screen.findByText('Revised classification to Treg CUT&Tag v3.'),
    ).toBeInTheDocument();

    await waitFor(() => {
      expect(catalogReviseExperiment).toHaveBeenCalledWith(
        EXPERIMENT_ID,
        expect.objectContaining({
          expected_revision_number: 2,
          display_name: 'Treg CUT&Tag v3',
          assay: {
            status: 'known',
            term_id: ASSAY_ID,
            term_revision_id: ASSAY_REVISION,
          },
          samples: [
            expect.objectContaining({
              sample_revision_id: SAMPLE_REVISION,
              replicate_role: {
                status: 'known',
                term_id: ROLE_ID,
                term_revision_id: ROLE_REVISION,
              },
            }),
          ],
        }),
        expect.objectContaining({
          'Idempotency-Key': expect.stringMatching(/^cat[A-Za-z0-9]+$/),
        }),
      );
    });
    expect(catalogCreateExperiment).not.toHaveBeenCalled();
    expect(catalogCreateTerm).not.toHaveBeenCalled();
  });

  it('hides classification writes from a member', async () => {
    vi.mocked(fetchSessionState).mockResolvedValue({
      ok: true,
      setup_required: false,
      authenticated: true,
      principal: {
        user_id: `usr_${'a'.repeat(32)}`,
        username: 'lab-member',
        role: 'member',
      },
    });
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });

    expect(
      await screen.findByText(/limited to administrators/i),
    ).toBeInTheDocument();
    expect(
      screen.queryByRole('form', { name: 'Record an experiment' }),
    ).not.toBeInTheDocument();
    expect(
      screen.queryByRole('form', { name: 'Revise experiment classification' }),
    ).not.toBeInTheDocument();
    await search(user);
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));
    expect(screen.queryByRole('button', { name: 'Ask Peri' })).toBeNull();
    expect(catalogListExperiments).toHaveBeenCalled();
    expect(catalogListSamples).not.toHaveBeenCalled();
  });

  it('pages revision history and edits the list head', async () => {
    const history = Array.from({ length: 24 }, (_, index) => {
      const number = index + 1;
      return experiment({
        experimentId: EXPERIMENT_ID,
        revisionId: `expr_${number.toString(16).padStart(32, '0')}`,
        revisionNumber: number,
        displayName: `History ${number}`,
      });
    });
    const headRevision = history[23];
    if (headRevision === undefined) throw new Error('missing head revision');
    vi.mocked(catalogListExperiments).mockResolvedValue({
      ok: true,
      project_id: PROJECT_ID,
      revision_view: 'linked',
      experiments: [headRevision],
      next_cursor: null,
    });
    vi.mocked(catalogListExperimentRevisions).mockImplementation(async (_id, params) => {
      if (params?.cursor === 'revcur_page2') {
        return { experiment_revisions: history.slice(20), next_cursor: null };
      }
      return {
        experiment_revisions: history.slice(0, 20),
        next_cursor: 'revcur_page2',
      };
    });
    vi.mocked(catalogReviseExperiment).mockResolvedValue({
      ...headRevision,
      experiment_revision_id: `expr_${'f'.repeat(32)}`,
      revision_number: 25,
      display_name: 'History 25',
    });

    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /History 24/ }));

    expect(await screen.findByRole('heading', { name: 'History 24' })).toBeInTheDocument();
    const detail = screen.getByRole('region', { name: 'Experiment revision' });
    expect(detail).toHaveTextContent('当前 24');
    expect(screen.getByText('Revision page 1')).toBeInTheDocument();
    expect(detail).toHaveTextContent('旧版 1: History 1');
    expect(detail.querySelectorAll('option')).toHaveLength(21);

    await user.click(screen.getByRole('button', { name: 'Next revision page' }));
    expect(await screen.findByText('Revision page 2')).toBeInTheDocument();
    expect(detail).toHaveTextContent('旧版 21: History 21');
    expect(catalogListExperimentRevisions).toHaveBeenCalledWith(EXPERIMENT_ID, {
      limit: 20,
      cursor: 'revcur_page2',
    });

    await user.selectOptions(
      screen.getByLabelText('Revision'),
      history[20]?.experiment_revision_id ?? '',
    );
    expect(await screen.findByRole('heading', { name: 'History 21' })).toBeInTheDocument();

    const revise = screen.getByRole('form', { name: 'Revise experiment classification' });
    const name = within(revise).getByLabelText('Display name');
    expect(name).toHaveValue('History 24');
    await user.clear(name);
    await user.type(name, 'History 25');
    await user.click(within(revise).getByRole('button', { name: 'Revise classification' }));
    expect(await screen.findByText('Revised classification to History 25.')).toBeInTheDocument();
    expect(catalogReviseExperiment).toHaveBeenCalledWith(
      EXPERIMENT_ID,
      expect.objectContaining({
        expected_revision_number: 24,
        display_name: 'History 25',
      }),
      expect.any(Object),
    );
  });

  it('reports a concurrent revision conflict without writing', async () => {
    const winner = experiment({
      experimentId: EXPERIMENT_ID,
      revisionId: `expr_${'f'.repeat(32)}`,
      revisionNumber: 3,
      displayName: 'Winner',
    });
    let conflicted = false;
    vi.mocked(catalogReviseExperiment).mockImplementation(async () => {
      conflicted = true;
      throw new ApiError(
        409,
        'catalog_revision_conflict',
        'The catalog revision changed before this write.',
      );
    });
    vi.mocked(catalogListExperiments).mockImplementation(async () => ({
      ok: true,
      project_id: PROJECT_ID,
      revision_view: 'linked',
      experiments: conflicted ? [winner] : [currentRevision],
      next_cursor: null,
    }));
    vi.mocked(catalogListExperimentRevisions).mockImplementation(async () => ({
      experiment_revisions: conflicted
        ? [olderRevision, currentRevision, winner]
        : [olderRevision, currentRevision],
      next_cursor: null,
    }));
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));
    const revise = await screen.findByRole('form', {
      name: 'Revise experiment classification',
    });
    await user.clear(within(revise).getByLabelText('Display name'));
    await user.type(within(revise).getByLabelText('Display name'), 'Losing name');
    await user.click(within(revise).getByRole('button', { name: 'Revise classification' }));

    expect(await screen.findByRole('alert')).toHaveTextContent('冲突。这次提交没有写入。');
    expect(await screen.findByRole('heading', { name: 'Winner' })).toBeInTheDocument();
    expect(screen.getByRole('region', { name: 'Experiment revision' })).toHaveTextContent(
      '当前 3',
    );
    expect(screen.queryByText('Revised classification to Losing name.')).not.toBeInTheDocument();
    expect(catalogReviseExperiment).toHaveBeenCalledWith(
      EXPERIMENT_ID,
      expect.objectContaining({ expected_revision_number: 2 }),
      expect.any(Object),
    );
  });

  it('keeps the list head editable while later revision pages remain', async () => {
    const head = experiment({
      experimentId: EXPERIMENT_ID,
      revisionId: `expr_${(401).toString(16).padStart(32, '0')}`,
      revisionNumber: 401,
      displayName: 'History 401',
    });
    vi.mocked(catalogListExperiments).mockResolvedValue({
      ok: true,
      project_id: PROJECT_ID,
      revision_view: 'linked',
      experiments: [head],
      next_cursor: null,
    });
    vi.mocked(catalogListExperimentRevisions).mockResolvedValue({
      experiment_revisions: [
        experiment({
          experimentId: EXPERIMENT_ID,
          revisionId: `expr_${'1'.padStart(32, '0')}`,
          revisionNumber: 1,
          displayName: 'History 1',
        }),
      ],
      next_cursor: 'revcur_more',
    });
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /History 401/ }));

    expect(await screen.findByRole('heading', { name: 'History 401' })).toBeInTheDocument();
    const detail = screen.getByRole('region', { name: 'Experiment revision' });
    expect(detail).toHaveTextContent('当前 401');
    expect(screen.getByRole('button', { name: 'Next revision page' })).toBeEnabled();
    expect(
      within(screen.getByRole('form', { name: 'Revise experiment classification' })).getByRole(
        'button',
        { name: 'Revise classification' },
      ),
    ).toBeEnabled();
    expect(catalogListExperimentRevisions).toHaveBeenCalledWith(EXPERIMENT_ID, {
      limit: 20,
    });
    expect(catalogListExperimentRevisions).not.toHaveBeenCalledWith(
      EXPERIMENT_ID,
      expect.objectContaining({ cursor: 'revcur_more' }),
    );
  });

  it('keeps the opened experiment when the list moves to another page', async () => {
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));
    expect(
      await screen.findByRole('heading', { name: 'Treg CUT&Tag v2' }),
    ).toBeInTheDocument();

    await user.click(screen.getByRole('button', { name: 'Next page' }));
    expect(
      await screen.findByRole('button', { name: /CUT&Tag page 20/ }),
    ).toBeInTheDocument();
    expect(screen.getByRole('heading', { name: 'Treg CUT&Tag v2' })).toBeInTheDocument();
    expect(screen.getByRole('region', { name: 'Experiment revision' })).toBeInTheDocument();
  });

  it('shows a sample load failure instead of an empty project', async () => {
    const user = userEvent.setup();
    vi.mocked(catalogListSamples).mockRejectedValue(
      new ApiError(503, 'catalog_unavailable', 'Sample revisions could not be loaded.'),
    );
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);

    const form = await screen.findByRole('form', { name: 'Record an experiment' });
    expect(await within(form).findByRole('alert')).toHaveTextContent(
      'Sample revisions could not be loaded.',
    );
    expect(within(form).queryByText(/no sample revisions/i)).not.toBeInTheDocument();
    expect(within(form).getByRole('button', { name: 'Record experiment' })).toBeDisabled();
  });

  it('shows an empty project only after the sample list succeeds with no rows', async () => {
    const user = userEvent.setup();
    vi.mocked(catalogListSamples).mockResolvedValue({
      project_id: PROJECT_ID,
      samples: [],
      next_cursor: null,
    });
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);

    const form = await screen.findByRole('form', { name: 'Record an experiment' });
    expect(await within(form).findByText(/no sample revisions/i)).toBeInTheDocument();
    expect(within(form).queryByRole('alert')).not.toBeInTheDocument();
    expect(within(form).getByRole('button', { name: 'Record experiment' })).toBeDisabled();
  });

  it('follows every sample-run page and shows a shared run once', async () => {
    const secondSampleId = `smp_${'8'.repeat(32)}`;
    const secondSampleRevision = `smpr_${'9'.repeat(32)}`;
    const firstSample = currentRevision.samples[0];
    if (firstSample === undefined) throw new Error('fixture sample missing');
    const twoSamples = {
      ...currentRevision,
      samples: [
        firstSample,
        {
          ...firstSample,
          sample_id: secondSampleId,
          sample_revision_id: secondSampleRevision,
          ordinal: 1,
        },
      ],
    };
    const laterRun: CatalogPublishedRunResponse = {
      ...publishedRun,
      artifact_id: 'signal',
      artifact_generation: `artifactgen-${'c'.repeat(64)}`,
    };
    const unpublished = {
      run_id: 'run_shared_unpublished',
      provenance: 'resolved',
      publication_evidence: 'absent' as const,
    };
    vi.mocked(catalogListExperiments).mockResolvedValue({
      ok: true,
      project_id: PROJECT_ID,
      revision_view: 'linked',
      experiments: [twoSamples],
      next_cursor: null,
    });
    vi.mocked(catalogListExperimentRevisions).mockResolvedValue({
      experiment_revisions: [twoSamples],
      next_cursor: null,
    });
    vi.mocked(catalogListSampleRevisionRunResults).mockImplementation(
      async (_projectId, _sampleRevisionId, params) => {
        if (params?.cursor === 'runcur_2') {
          return {
            ok: true,
            project_id: PROJECT_ID,
            sample_revision_id: _sampleRevisionId,
            results: [laterRun],
            next_cursor: null,
          };
        }
        return {
          ok: true,
          project_id: PROJECT_ID,
          sample_revision_id: _sampleRevisionId,
          results: [publishedRun, unpublished],
          next_cursor: 'runcur_2',
        };
      },
    );
    const user = userEvent.setup();
    renderWithRouter(routes, { initialEntries: ['/catalog'] });
    await search(user);
    await user.click(await screen.findByRole('button', { name: /Treg CUT&Tag v2/ }));

    const scope = await screen.findByRole('region', { name: 'Run input scope' });
    expect(await within(scope).findByText('signal')).toBeInTheDocument();
    expect(within(scope).getAllByText('peaks')).toHaveLength(1);
    expect(within(scope).getAllByText(/run_shared_unpublished/)).toHaveLength(1);
    expect(catalogListSampleRevisionRunResults).toHaveBeenCalledWith(
      PROJECT_ID,
      SAMPLE_REVISION,
      { limit: 20, cursor: 'runcur_2' },
    );
    expect(catalogListSampleRevisionRunResults).toHaveBeenCalledWith(
      PROJECT_ID,
      secondSampleRevision,
      { limit: 20, cursor: 'runcur_2' },
    );
  });
});
