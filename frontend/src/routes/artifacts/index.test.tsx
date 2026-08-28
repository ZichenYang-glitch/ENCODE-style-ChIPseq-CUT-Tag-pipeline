import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import type {
  ArtifactPublicationDetailResponse,
  ArtifactPublicationListResponse,
  ArtifactPublicationResponse,
} from '../../api/generated/models';
import {
  getArtifactPublication,
  listArtifactPublications,
} from '../../api/generated/artifact-publications/artifact-publications';
import { downloadRunArtifact } from '../../api/generated/artifacts/artifacts';
import { ApiError } from '../../api/fetcher';
import {
  safeArtifactDownloadFilename,
  saveArtifactBlob,
} from '../../features/run-artifacts/artifactDownload';
import { appRoutes } from '../../app/router';
import { renderWithRouter } from '../../test/test-utils';

vi.mock('../../api/generated/artifact-publications/artifact-publications', () => ({
  listArtifactPublications: vi.fn(),
  getArtifactPublication: vi.fn(),
}));
vi.mock('../../api/generated/artifacts/artifacts', () => ({
  downloadRunArtifact: vi.fn(),
}));
vi.mock('../../features/run-artifacts/artifactDownload', () => ({
  safeArtifactDownloadFilename: vi.fn(() => 'artifact-a.download'),
  saveArtifactBlob: vi.fn(),
}));

const listPublicationsMock = vi.mocked(listArtifactPublications);
const getPublicationMock = vi.mocked(getArtifactPublication);
const downloadArtifactMock = vi.mocked(downloadRunArtifact);
const safeFilenameMock = vi.mocked(safeArtifactDownloadFilename);
const saveBlobMock = vi.mocked(saveArtifactBlob);

const PROJECT_ID = `prj_${'1'.repeat(32)}`;
const SAMPLE_ID = `smp_${'2'.repeat(32)}`;
const SAMPLE_REVISION_ID = `smpr_${'3'.repeat(32)}`;
const GENERATION_A = `artifactgen-${'a'.repeat(64)}`;
const GENERATION_B = `artifactgen-${'b'.repeat(64)}`;
const REVISION = `artifactrev-${'c'.repeat(64)}`;

function publication(
  overrides: Partial<ArtifactPublicationResponse> = {},
): ArtifactPublicationResponse {
  return {
    run_id: 'run-a',
    project_id: PROJECT_ID,
    workflow_id: 'workflow-a',
    artifact_id: 'artifact-a',
    output_type: 'counts',
    artifact_generation: GENERATION_A,
    artifact_revision: REVISION,
    published_at: '2026-08-07T08:00:00Z',
    current_artifact_generation: GENERATION_A,
    generation_status: 'current',
    run_sample_binding: {
      binding_mode: 'bound_v1',
      provenance: 'resolved',
      associated_run_samples: [
        {
          sample_id: SAMPLE_ID,
          sample_revision_id: SAMPLE_REVISION_ID,
          revision_number: 4,
          ordinal: 0,
        },
      ],
    },
    ...overrides,
  };
}

function listResponse(
  publications: ArtifactPublicationResponse[],
  nextCursor: string | null = null,
): ArtifactPublicationListResponse {
  return { ok: true, publications, next_cursor: nextCursor, issues: [] };
}

function detailResponse(value: ArtifactPublicationResponse): ArtifactPublicationDetailResponse {
  return { ok: true, publication: value, issues: [] };
}

beforeAll(async () => {
  await import('./index');
});

beforeEach(() => {
  listPublicationsMock.mockReset();
  getPublicationMock.mockReset();
  downloadArtifactMock.mockReset();
  safeFilenameMock.mockClear();
  saveBlobMock.mockClear();
  listPublicationsMock.mockResolvedValue(listResponse([publication()]));
  getPublicationMock.mockResolvedValue(detailResponse(publication()));
  downloadArtifactMock.mockResolvedValue(new Blob(['artifact']));
});

describe('artifact publication routes', () => {
  it('applies every public filter, clears old rows, and uses the exact API params', async () => {
    const user = userEvent.setup();
    listPublicationsMock.mockImplementation(async (params) =>
      listResponse([
        publication({
          run_id: params?.run_id ?? 'run-initial',
          artifact_id: params?.run_id ? 'artifact-filtered' : 'artifact-initial',
        }),
      ]),
    );
    renderWithRouter(appRoutes, { initialEntries: ['/artifacts'] });
    expect(await screen.findAllByRole('link', {
      name: 'Open artifact ID artifact-initial',
    })).not.toHaveLength(0);
    expect(screen.getAllByText('Associated run samples')).not.toHaveLength(0);
    expect(screen.getAllByRole('button', {
      name: `Copy full sample revision ID ${SAMPLE_REVISION_ID}`,
    })).not.toHaveLength(0);
    expect(screen.getByText(
      'These are input samples associated with the run; they are not per-artifact sample attribution.',
    )).toBeInTheDocument();

    await user.click(screen.getByText('Filters (0 active)'));
    await user.type(screen.getByLabelText('Project'), PROJECT_ID);
    await user.type(screen.getByLabelText('Run'), 'run-filtered');
    await user.type(screen.getByLabelText('Workflow'), 'workflow-filtered');
    await user.type(screen.getByLabelText('Artifact type'), 'counts');
    await user.type(
      screen.getByLabelText('Associated run sample revision'),
      SAMPLE_REVISION_ID,
    );
    await user.type(
      screen.getByLabelText('Published from'),
      '2026-08-07T16:00:00.123456+08:00',
    );
    await user.type(
      screen.getByLabelText('Published before'),
      '2026-08-07T16:00:00.123457+08:00',
    );
    await user.click(screen.getByLabelText('Current publications only'));
    await user.selectOptions(screen.getByLabelText('Rows per page'), '100');
    await user.click(screen.getByRole('button', { name: 'Apply filters' }));

    expect(screen.queryAllByRole('link', {
      name: 'Open artifact ID artifact-initial',
    })).toHaveLength(0);
    expect(await screen.findAllByRole('link', {
      name: 'Open artifact ID artifact-filtered',
    })).not.toHaveLength(0);
    expect(listPublicationsMock).toHaveBeenLastCalledWith({
      project_id: PROJECT_ID,
      run_id: 'run-filtered',
      workflow_id: 'workflow-filtered',
      output_type: 'counts',
      associated_run_sample_revision_id: SAMPLE_REVISION_ID,
      published_from: '2026-08-07T08:00:00.123456Z',
      published_before: '2026-08-07T08:00:00.123457Z',
      current_only: false,
      after: undefined,
      limit: 100,
    });
  });

  it('blocks repeated unsafe URL values and redacts exception text', async () => {
    const first = renderWithRouter(appRoutes, {
      initialEntries: ['/artifacts?run_id=run-a&run_id=run-b'],
    });
    expect(
      await screen.findByRole('heading', { name: 'Artifact filters could not be used' }),
    ).toBeInTheDocument();
    expect(listPublicationsMock).not.toHaveBeenCalled();
    first.unmount();

    listPublicationsMock.mockRejectedValue(new Error('file:///private/results SECRET_TOKEN=x'));
    renderWithRouter(appRoutes, { initialEntries: ['/artifacts'] });
    expect(
      await screen.findByRole('heading', { name: 'Artifact publications could not be loaded' }),
    ).toBeInTheDocument();
    expect(screen.queryByText(/private\/results|SECRET_TOKEN/)).not.toBeInTheDocument();
  });

  it('renders run-level sample binding semantics and downloads current revision-bound bytes', async () => {
    const user = userEvent.setup();
    const blob = new Blob(['current artifact']);
    downloadArtifactMock.mockResolvedValue(blob);
    const current = renderWithRouter(appRoutes, {
      initialEntries: [`/artifacts/run-a/artifact-a?generation=${GENERATION_A}`],
    });

    expect(await screen.findByRole('heading', { name: 'Associated run samples' }))
      .toBeInTheDocument();
    expect(screen.getByText(
      'These are input samples associated with the run; they are not per-artifact sample attribution.',
    )).toBeInTheDocument();
    expect(screen.getAllByRole('button', {
      name: `Copy full sample revision ID ${SAMPLE_REVISION_ID}`,
    })).not.toHaveLength(0);

    await user.click(screen.getByRole('button', { name: 'Download artifact artifact-a' }));
    await waitFor(() => expect(saveBlobMock).toHaveBeenCalledWith(blob, 'artifact-a.download'));
    expect(downloadArtifactMock).toHaveBeenCalledWith('run-a', 'artifact-a', {
      generation: GENERATION_A,
      revision: REVISION,
    });
    expect(safeFilenameMock).toHaveBeenCalledWith('', 'artifact-a');
    current.unmount();

    getPublicationMock.mockResolvedValue(
      detailResponse(
        publication({
          project_id: `prj_${'0'.repeat(32)}`,
          run_sample_binding: {
            binding_mode: 'legacy_v1',
            provenance: 'unresolved',
            associated_run_samples: [],
          },
        }),
      ),
    );
    renderWithRouter(appRoutes, {
      initialEntries: [`/artifacts/run-a/artifact-a?generation=${GENERATION_A}`],
    });
    expect(
      await screen.findByText(
        'No associated run samples were recorded for this run binding.',
      ),
    ).toBeInTheDocument();
  });

  it('keeps superseded publications metadata-only and labels an HTTP 409 stale', async () => {
    const sampleRevisionB = `smpr_${'4'.repeat(32)}`;
    const sampleRevisionC = `smpr_${'5'.repeat(32)}`;
    const superseded = publication({
      current_artifact_generation: GENERATION_B,
      generation_status: 'superseded',
      run_sample_binding: {
        binding_mode: 'bound_v1',
        provenance: 'resolved',
        associated_run_samples: [
          {
            sample_id: SAMPLE_ID,
            sample_revision_id: SAMPLE_REVISION_ID,
            revision_number: 4,
            ordinal: 0,
          },
          {
            sample_id: `smp_${'4'.repeat(32)}`,
            sample_revision_id: sampleRevisionB,
            revision_number: 1,
            ordinal: 1,
          },
          {
            sample_id: `smp_${'5'.repeat(32)}`,
            sample_revision_id: sampleRevisionC,
            revision_number: 2,
            ordinal: 2,
          },
        ],
      },
    });
    listPublicationsMock.mockResolvedValue(listResponse([superseded]));
    const list = renderWithRouter(appRoutes, { initialEntries: ['/artifacts'] });
    expect(await screen.findAllByText(/Current:/)).not.toHaveLength(0);
    expect(screen.getAllByRole('button', {
      name: `Copy full current artifact generation ${GENERATION_B}`,
    })).not.toHaveLength(0);
    expect(screen.getAllByText('+1 more associated run samples')).not.toHaveLength(0);
    expect(screen.queryByText(sampleRevisionC)).not.toBeInTheDocument();
    list.unmount();

    getPublicationMock.mockResolvedValue(detailResponse(superseded));
    const first = renderWithRouter(appRoutes, {
      initialEntries: [`/artifacts/run-a/artifact-a?generation=${GENERATION_A}`],
    });
    expect(await screen.findByText('Superseded · historical metadata only')).toBeInTheDocument();
    expect(screen.getByText('Current generation')).toBeInTheDocument();
    expect(screen.getByRole('button', {
      name: `Copy full current artifact generation ${GENERATION_B}`,
    })).toBeInTheDocument();
    expect(screen.getByRole('button', {
      name: `Copy full sample revision ID ${sampleRevisionC}`,
    })).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: /Download artifact/ })).not.toBeInTheDocument();
    first.unmount();

    getPublicationMock.mockResolvedValue(detailResponse(publication()));
    downloadArtifactMock.mockRejectedValue(
      new ApiError(409, 'RUN_ARTIFACT_DOWNLOAD_CONFLICT', '/private/path', [
        { code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT', message: 'private', path: 'generation' },
      ]),
    );
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/artifacts/run-a/artifact-a?generation=${GENERATION_A}`],
    });
    const download = await screen.findByRole('button', {
      name: 'Download artifact artifact-a',
    });
    await user.click(download);
    expect(await screen.findByText('Stale generation')).toBeInTheDocument();
    expect(download).toBeDisabled();
    expect(screen.queryByText('/private/path')).not.toBeInTheDocument();
    await waitFor(() => expect(getPublicationMock.mock.calls.length).toBeGreaterThan(1));
  });

  it.each(['revision', 'artifact_id'])(
    'does not label a %s download conflict as a stale generation',
    async (path) => {
      downloadArtifactMock.mockRejectedValue(
        new ApiError(409, 'RUN_ARTIFACT_DOWNLOAD_CONFLICT', '/private/path', [
          { code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT', message: 'private', path },
        ]),
      );
      const user = userEvent.setup();
      renderWithRouter(appRoutes, {
        initialEntries: [`/artifacts/run-a/artifact-a?generation=${GENERATION_A}`],
      });
      const download = await screen.findByRole('button', {
        name: 'Download artifact artifact-a',
      });
      await user.click(download);

      expect(await screen.findByText(
        'Download could not be completed. Retry when the current artifact is available.',
      )).toBeInTheDocument();
      expect(screen.queryByText('Stale generation')).not.toBeInTheDocument();
      expect(download).toBeEnabled();
      expect(screen.queryByText('/private/path')).not.toBeInTheDocument();
      expect(getPublicationMock).toHaveBeenCalledTimes(1);
    },
  );

  it('fails closed on detail identity mismatch without rendering private fields', async () => {
    getPublicationMock.mockResolvedValue(
      detailResponse({
        ...publication({ run_id: 'run-other' }),
        // Runtime additions must never be rendered by the publication projection.
        uri: 'file:///private/artifact',
        filename: 'secret.tsv',
      } as ArtifactPublicationResponse),
    );
    renderWithRouter(appRoutes, {
      initialEntries: [`/artifacts/run-a/artifact-a?generation=${GENERATION_A}`],
    });
    expect(
      await screen.findByRole('heading', { name: 'Artifact publication could not be loaded' }),
    ).toBeInTheDocument();
    expect(screen.queryByText(/private\/artifact|secret\.tsv/)).not.toBeInTheDocument();
  });
});
