import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import type { ComponentProps } from 'react';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type {
  ArtifactReferenceResponse,
  RunArtifactDetailResponse,
  RunArtifactsResponse,
} from '../../api/generated/models';
import {
  downloadRunArtifact,
  getRunArtifact,
  listRunArtifacts,
} from '../../api/generated/artifacts/artifacts';
import { ApiError } from '../../api/fetcher';
import { ArtifactBrowser } from './ArtifactBrowser';

vi.mock('../../api/generated/artifacts/artifacts', () => ({
  listRunArtifacts: vi.fn(),
  getRunArtifact: vi.fn(),
  downloadRunArtifact: vi.fn(),
}));

const listArtifactsMock = vi.mocked(listRunArtifacts);
const getArtifactMock = vi.mocked(getRunArtifact);
const downloadArtifactMock = vi.mocked(downloadRunArtifact);
const GENERATION_A = `artifactgen-${'a'.repeat(64)}`;
const GENERATION_B = `artifactgen-${'b'.repeat(64)}`;
const CURSOR_A = `artifactcur_${'a'.repeat(64)}`;

function artifact(id: string): ArtifactReferenceResponse {
  return {
    artifact_id: id,
    run_id: 'run-1',
    artifact_type: 'file',
    name: `${id}.tsv`,
    uri: `run://runs/run-1/artifacts/${id}`,
    mime_type: 'text/tab-separated-values',
    produced_at: '2026-07-12T12:00:00Z',
    revision: `artifactrev-${'a'.repeat(64)}`,
    relative_path: `results/${id}.tsv`,
    output_type: id,
    size_bytes: 42,
    metadata: { scope: 'project' },
  };
}

function page(
  artifacts: ArtifactReferenceResponse[],
  nextCursor: string | null = null,
  generation = GENERATION_A,
): RunArtifactsResponse {
  return {
    ok: true,
    run_id: 'run-1',
    artifact_generation: generation,
    artifacts,
    next_cursor: nextCursor,
    issues: [],
  };
}

function detail(
  value: ArtifactReferenceResponse,
  generation = GENERATION_A,
): RunArtifactDetailResponse {
  return {
    ok: true,
    run_id: 'run-1',
    artifact_generation: generation,
    artifact: value,
    issues: [],
  };
}

function renderBrowser(
  overrides: Partial<ComponentProps<typeof ArtifactBrowser>> = {},
) {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false }, mutations: { retry: false } },
  });
  const props: ComponentProps<typeof ArtifactBrowser> = {
    runId: 'run-1',
    runStatus: 'succeeded',
    outcome: { kind: 'indexed', count: 1, generation: GENERATION_A },
    selectedArtifactId: null,
    onSelectArtifact: vi.fn(),
    onRefreshStatus: vi.fn(),
    ...overrides,
  };
  return {
    props,
    queryClient,
    ...render(
      <QueryClientProvider client={queryClient}>
        <ArtifactBrowser {...props} />
      </QueryClientProvider>,
    ),
  };
}

beforeEach(() => {
  listArtifactsMock.mockReset();
  getArtifactMock.mockReset();
  downloadArtifactMock.mockReset();
});

afterEach(() => {
  vi.clearAllMocks();
});

describe('ArtifactBrowser honest states', () => {
  it('does not query before success or while indexing', () => {
    const first = renderBrowser({ runStatus: 'running', outcome: { kind: 'pending' } });
    expect(screen.getByText('Artifacts are not available yet')).toBeInTheDocument();
    first.unmount();
    renderBrowser({ outcome: { kind: 'pending' } });
    expect(screen.getByText('Indexing artifacts')).toBeInTheDocument();
    expect(listArtifactsMock).not.toHaveBeenCalled();
  });

  it('shows unconfirmed and failed states with canonical status refresh', async () => {
    const user = userEvent.setup();
    const onRefreshStatus = vi.fn();
    const first = renderBrowser({ outcome: { kind: 'unconfirmed' }, onRefreshStatus });
    await user.click(screen.getByRole('button', { name: 'Refresh status' }));
    expect(onRefreshStatus).toHaveBeenCalledTimes(1);
    first.unmount();

    renderBrowser({
      outcome: { kind: 'failed', reasonCode: 'ARTIFACT_EXTRACTION_ADAPTER_FAILED' },
      onRefreshStatus,
    });
    expect(screen.getByText('Artifact indexing failed')).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: /retry extraction/i })).not.toBeInTheDocument();
  });

  it('shows true empty only after indexed zero and a successful empty response', async () => {
    listArtifactsMock.mockResolvedValue(page([]));
    renderBrowser({
      outcome: { kind: 'indexed', count: 0, generation: GENERATION_A },
    });
    expect(await screen.findByText('No indexed artifacts')).toBeInTheDocument();
  });

  it('keeps nonzero empty responses honest', async () => {
    listArtifactsMock.mockResolvedValue(page([]));
    renderBrowser({
      outcome: { kind: 'indexed', count: 1, generation: GENERATION_A },
    });
    expect(await screen.findByText('Artifact index is not ready')).toBeInTheDocument();
    expect(screen.queryByText('No indexed artifacts')).not.toBeInTheDocument();
  });
});

describe('ArtifactBrowser queries', () => {
  it('loads explicit pages of 50 and stops a repeated cursor', async () => {
    const user = userEvent.setup();
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], CURSOR_A))
      .mockResolvedValueOnce(page([artifact('artifact-b')], CURSOR_A));
    renderBrowser({
      outcome: { kind: 'indexed', count: 2, generation: GENERATION_A },
    });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));
    expect(await screen.findAllByText('artifact-b')).not.toHaveLength(0);
    expect(listArtifactsMock).toHaveBeenNthCalledWith(
      1,
      'run-1',
      { after: undefined, generation: GENERATION_A, limit: 50 },
    );
    expect(listArtifactsMock).toHaveBeenNthCalledWith(
      2,
      'run-1',
      { after: CURSOR_A, generation: GENERATION_A, limit: 50 },
    );
    expect(screen.queryByRole('button', { name: 'Load more artifacts' })).not.toBeInTheDocument();
  });

  it('preserves existing rows when another page fails and offers retry', async () => {
    const user = userEvent.setup();
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], CURSOR_A))
      .mockRejectedValueOnce(new Error('/private/path'));
    renderBrowser({
      outcome: { kind: 'indexed', count: 2, generation: GENERATION_A },
    });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));
    expect(await screen.findByText(/Existing rows are preserved/)).toBeInTheDocument();
    expect(screen.getAllByText('artifact-a').length).toBeGreaterThan(0);
    expect(screen.queryByText('/private/path')).not.toBeInTheDocument();
  });

  it('discards generation-conflicted pages and refreshes status without retrying the old generation', async () => {
    const user = userEvent.setup();
    const changed = new ApiError(
      409,
      'RUN_ARTIFACT_GENERATION_CHANGED',
      'Artifact generation changed.',
    );
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], CURSOR_A))
      .mockRejectedValueOnce(changed)
      .mockResolvedValueOnce(
        page([artifact('artifact-b')], null, GENERATION_B),
      );
    let rendered: ReturnType<typeof renderBrowser>;
    const onRefreshStatus = vi
      .fn<[], Promise<void>>()
      .mockRejectedValueOnce(new Error('/private/status'))
      .mockImplementationOnce(async () => {
        rendered.rerender(
          <QueryClientProvider client={rendered.queryClient}>
            <ArtifactBrowser
              {...rendered.props}
              outcome={{
                kind: 'indexed',
                count: 1,
                generation: GENERATION_B,
              }}
            />
          </QueryClientProvider>,
        );
      });
    rendered = renderBrowser({
      outcome: { kind: 'indexed', count: 2, generation: GENERATION_A },
      onRefreshStatus,
    });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));

    expect(
      await screen.findByText('Artifact generation changed'),
    ).toBeInTheDocument();
    expect(screen.queryAllByText('artifact-a')).toHaveLength(0);
    expect(onRefreshStatus).toHaveBeenCalledTimes(1);
    expect(
      listArtifactsMock.mock.calls.filter(
        ([, parameters]) => parameters?.generation === GENERATION_A,
      ),
    ).toHaveLength(2);

    await user.click(screen.getByRole('button', { name: 'Refresh status' }));
    expect(onRefreshStatus).toHaveBeenCalledTimes(2);
    expect(await screen.findAllByText('artifact-b')).not.toHaveLength(0);
    expect(listArtifactsMock).toHaveBeenLastCalledWith('run-1', {
      after: undefined,
      generation: GENERATION_B,
      limit: 50,
    });
  });

  it('fails closed before starting canonical refresh for a generation-conflicted page', async () => {
    const user = userEvent.setup();
    const changed = new ApiError(
      409,
      'RUN_ARTIFACT_GENERATION_CHANGED',
      'Artifact generation changed.',
    );
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], CURSOR_A))
      .mockRejectedValueOnce(changed);
    let staleArtifactWasVisibleAtRefresh = true;
    const onRefreshStatus = vi.fn(async () => {
      staleArtifactWasVisibleAtRefresh =
        screen.queryAllByText('artifact-a').length > 0;
    });
    renderBrowser({
      outcome: { kind: 'indexed', count: 2, generation: GENERATION_A },
      onRefreshStatus,
    });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));

    await waitFor(() => expect(onRefreshStatus).toHaveBeenCalledTimes(1));
    expect(staleArtifactWasVisibleAtRefresh).toBe(false);
    expect(screen.queryAllByText('artifact-a')).toHaveLength(0);
  });

  it('keeps a conflicted generation quarantined across unconfirmed status snapshots', async () => {
    const user = userEvent.setup();
    const changed = new ApiError(
      409,
      'RUN_ARTIFACT_GENERATION_CHANGED',
      'Artifact generation changed.',
    );
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], CURSOR_A))
      .mockRejectedValueOnce(changed)
      .mockResolvedValueOnce(
        page([artifact('artifact-b')], null, GENERATION_B),
      );
    const onRefreshStatus = vi.fn().mockResolvedValue(undefined);
    const rendered = renderBrowser({
      outcome: { kind: 'indexed', count: 2, generation: GENERATION_A },
      onRefreshStatus,
    });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));
    expect(await screen.findByText('Artifact generation changed')).toBeInTheDocument();

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <ArtifactBrowser {...rendered.props} outcome={{ kind: 'unconfirmed' }} />
      </QueryClientProvider>,
    );
    expect(
      await screen.findByText('Artifact status could not be confirmed'),
    ).toBeInTheDocument();

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <ArtifactBrowser
          {...rendered.props}
          outcome={{ kind: 'indexed', count: 2, generation: GENERATION_A }}
        />
      </QueryClientProvider>,
    );
    await act(async () => {
      await new Promise((resolve) => setTimeout(resolve, 0));
    });
    expect(
      listArtifactsMock.mock.calls.filter(
        ([, parameters]) => parameters?.generation === GENERATION_A,
      ),
    ).toHaveLength(2);
    expect(screen.getByText('Artifact generation changed')).toBeInTheDocument();

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <ArtifactBrowser
          {...rendered.props}
          outcome={{ kind: 'indexed', count: 1, generation: GENERATION_B }}
        />
      </QueryClientProvider>,
    );
    expect(await screen.findAllByText('artifact-b')).not.toHaveLength(0);
    expect(listArtifactsMock).toHaveBeenLastCalledWith('run-1', {
      after: undefined,
      generation: GENERATION_B,
      limit: 50,
    });
  });

  it('retries an initial redacted list failure', async () => {
    const user = userEvent.setup();
    listArtifactsMock
      .mockRejectedValueOnce(new Error('/private/path'))
      .mockRejectedValueOnce(new Error('/private/path'))
      .mockResolvedValueOnce(page([artifact('artifact-a')]));
    renderBrowser();
    expect(await screen.findByText('Artifacts could not be loaded')).toBeInTheDocument();
    expect(screen.queryByText('/private/path')).not.toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Retry artifact list' }));
    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
  });

  it('restores valid deep-link detail and rejects unsafe IDs without a request', async () => {
    const selected = artifact('artifact-a');
    listArtifactsMock.mockResolvedValue(page([selected]));
    getArtifactMock.mockResolvedValue(detail(selected));
    const first = renderBrowser({ selectedArtifactId: selected.artifact_id });
    expect(await screen.findByText(selected.uri)).toBeInTheDocument();
    expect(getArtifactMock).toHaveBeenCalledWith('run-1', 'artifact-a', {
      generation: GENERATION_A,
    });
    first.unmount();

    listArtifactsMock.mockResolvedValue(page([selected]));
    renderBrowser({ selectedArtifactId: '../escape' });
    expect(await screen.findByText('This artifact link is invalid or unavailable.')).toBeInTheDocument();
    expect(getArtifactMock).toHaveBeenCalledTimes(1);
  });

  it('preserves confirmed detail when a background detail refresh fails', async () => {
    const selected = artifact('artifact-a');
    listArtifactsMock.mockResolvedValue(page([selected]));
    getArtifactMock
      .mockResolvedValueOnce(detail(selected))
      .mockRejectedValueOnce(new Error('/private/path'))
      .mockRejectedValueOnce(new Error('/private/path'));
    const { queryClient } = renderBrowser({
      selectedArtifactId: selected.artifact_id,
    });

    expect(await screen.findByText(selected.uri)).toBeInTheDocument();
    await queryClient.invalidateQueries({
      queryKey: ['run-artifact', 'run-1', GENERATION_A, 'artifact-a'],
    });

    expect(
      await screen.findByText(
        'Artifact detail refresh failed. The last confirmed detail is preserved.',
      ),
    ).toBeInTheDocument();
    expect(screen.getByText(selected.uri)).toBeInTheDocument();
    expect(screen.getAllByText(selected.relative_path).length).toBeGreaterThan(0);
    expect(screen.queryByText('/private/path')).not.toBeInTheDocument();
  });

  it('downloads through the generated operation and reports redacted failures', async () => {
    const user = userEvent.setup();
    const selected = artifact('artifact-a');
    const onRefreshStatus = vi.fn();
    listArtifactsMock.mockResolvedValue(page([selected]));
    getArtifactMock.mockResolvedValue(detail(selected));
    downloadArtifactMock
      .mockResolvedValueOnce(new Blob(['artifact bytes']))
      .mockRejectedValueOnce(new Error('/private/workspace/source.tsv'));
    const createObjectURL = vi.fn().mockReturnValue('blob:artifact');
    const revokeObjectURL = vi.fn();
    vi.stubGlobal('URL', { ...URL, createObjectURL, revokeObjectURL });
    const click = vi
      .spyOn(HTMLAnchorElement.prototype, 'click')
      .mockImplementation(() => undefined);
    renderBrowser({
      selectedArtifactId: selected.artifact_id,
      onRefreshStatus,
    });

    const button = await screen.findByRole('button', {
      name: 'Download artifact-a.tsv',
    });
    await user.click(button);
    expect(downloadArtifactMock).toHaveBeenCalledWith('run-1', 'artifact-a', {
      generation: GENERATION_A,
      revision: selected.revision,
    });
    expect(createObjectURL).toHaveBeenCalledTimes(1);
    expect(click).toHaveBeenCalledTimes(1);
    expect(revokeObjectURL).toHaveBeenCalledWith('blob:artifact');
    expect(await screen.findByText('Download prepared successfully.')).toBeInTheDocument();

    await user.click(button);
    expect(await screen.findByText(/Download could not be completed/)).toBeInTheDocument();
    expect(screen.queryByText('/private/workspace')).not.toBeInTheDocument();
    expect(onRefreshStatus).not.toHaveBeenCalled();
    expect(screen.getByText(selected.uri)).toBeInTheDocument();
  });

  it('quarantines the viewed generation when download reports a generation conflict', async () => {
    const user = userEvent.setup();
    const selected = artifact('artifact-a');
    const current = {
      ...selected,
      revision: `artifactrev-${'b'.repeat(64)}`,
      size_bytes: 84,
    };
    const onRefreshStatus = vi.fn();
    const changed = new ApiError(
      409,
      'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
      'Artifact content is no longer available as indexed.',
      [
        {
          code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
          message: 'Artifact content is no longer available as indexed.',
          path: 'generation',
        },
      ],
    );
    listArtifactsMock
      .mockResolvedValueOnce(page([selected]))
      .mockResolvedValueOnce(page([current], null, GENERATION_B));
    getArtifactMock
      .mockResolvedValueOnce(detail(selected))
      .mockResolvedValueOnce(detail(current, GENERATION_B));
    downloadArtifactMock.mockRejectedValue(changed);
    const rendered = renderBrowser({
      selectedArtifactId: selected.artifact_id,
      onRefreshStatus,
    });

    const button = await screen.findByRole('button', {
      name: 'Download artifact-a.tsv',
    });
    expect(screen.getByText(selected.uri)).toBeInTheDocument();
    await user.click(button);

    expect(
      await screen.findByText('Artifact generation changed'),
    ).toBeInTheDocument();
    expect(screen.queryByText(selected.uri)).not.toBeInTheDocument();
    expect(screen.queryAllByText(selected.relative_path)).toHaveLength(0);
    expect(onRefreshStatus).toHaveBeenCalledTimes(1);
    await waitFor(() => {
      expect(
        rendered.queryClient.getQueryCache().find({
          queryKey: ['run-artifacts', 'run-1', GENERATION_A],
          exact: true,
        })?.state.data,
      ).toBeUndefined();
      expect(
        rendered.queryClient.getQueryCache().find({
          queryKey: [
            'run-artifact',
            'run-1',
            GENERATION_A,
            selected.artifact_id,
          ],
          exact: true,
        })?.state.data,
      ).toBeUndefined();
    });

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <ArtifactBrowser
          {...rendered.props}
          outcome={{ kind: 'indexed', count: 1, generation: GENERATION_B }}
        />
      </QueryClientProvider>,
    );
    expect(await screen.findAllByText('84 B')).not.toHaveLength(0);
    expect(
      screen.queryByText('Artifact generation changed'),
    ).not.toBeInTheDocument();
    expect(onRefreshStatus).toHaveBeenCalledTimes(1);
  });

  it('preserves confirmed evidence for a non-generation download conflict', async () => {
    const user = userEvent.setup();
    const selected = artifact('artifact-a');
    const onRefreshStatus = vi.fn();
    const changed = new ApiError(
      409,
      'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
      'Artifact content is no longer available as indexed.',
      [
        {
          code: 'RUN_ARTIFACT_DOWNLOAD_CONFLICT',
          message: 'Artifact content is no longer available as indexed.',
          path: 'artifact_id',
        },
      ],
    );
    listArtifactsMock.mockResolvedValue(page([selected]));
    getArtifactMock.mockResolvedValue(detail(selected));
    downloadArtifactMock.mockRejectedValue(changed);
    renderBrowser({
      selectedArtifactId: selected.artifact_id,
      onRefreshStatus,
    });

    await user.click(
      await screen.findByRole('button', { name: 'Download artifact-a.tsv' }),
    );

    expect(await screen.findByText(/Download could not be completed/)).toBeInTheDocument();
    expect(screen.getByText(selected.uri)).toBeInTheDocument();
    expect(onRefreshStatus).not.toHaveBeenCalled();
    expect(
      screen.queryByText('Artifact generation changed'),
    ).not.toBeInTheDocument();
  });

  it('clears equal-count list and detail caches when generation changes', async () => {
    const oldArtifact = artifact('artifact-a');
    const currentArtifact = {
      ...oldArtifact,
      revision: `artifactrev-${'b'.repeat(64)}`,
      size_bytes: 84,
    };
    listArtifactsMock
      .mockResolvedValueOnce(page([oldArtifact], null, GENERATION_A))
      .mockResolvedValueOnce(page([currentArtifact], null, GENERATION_B));
    getArtifactMock
      .mockResolvedValueOnce(detail(oldArtifact, GENERATION_A))
      .mockResolvedValueOnce(detail(currentArtifact, GENERATION_B));
    const rendered = renderBrowser({ selectedArtifactId: 'artifact-a' });
    expect(await screen.findAllByText('42 B')).not.toHaveLength(0);

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <ArtifactBrowser
          {...rendered.props}
          outcome={{ kind: 'indexed', count: 1, generation: GENERATION_B }}
        />
      </QueryClientProvider>,
    );

    await waitFor(() => expect(screen.queryByText('42 B')).not.toBeInTheDocument());
    expect(await screen.findAllByText('84 B')).not.toHaveLength(0);
    expect(
      rendered.queryClient
        .getQueryCache()
        .findAll({ queryKey: ['run-artifacts', 'run-1'] })
        .map((query) => query.queryKey),
    ).toEqual([['run-artifacts', 'run-1', GENERATION_B]]);
    expect(
      rendered.queryClient
        .getQueryCache()
        .findAll({ queryKey: ['run-artifact', 'run-1'] })
        .map((query) => query.queryKey),
    ).toEqual([['run-artifact', 'run-1', GENERATION_B, 'artifact-a']]);
  });

  it('does not preserve a deleted detail from an older generation', async () => {
    const selected = artifact('artifact-a');
    listArtifactsMock
      .mockResolvedValueOnce(page([selected], null, GENERATION_A))
      .mockResolvedValueOnce(page([], null, GENERATION_B));
    getArtifactMock.mockResolvedValueOnce(detail(selected, GENERATION_A));
    const rendered = renderBrowser({ selectedArtifactId: 'artifact-a' });
    expect(await screen.findByText(selected.uri)).toBeInTheDocument();

    rendered.rerender(
      <QueryClientProvider client={rendered.queryClient}>
        <ArtifactBrowser
          {...rendered.props}
          outcome={{ kind: 'indexed', count: 0, generation: GENERATION_B }}
        />
      </QueryClientProvider>,
    );

    expect(await screen.findByText('No indexed artifacts')).toBeInTheDocument();
    expect(screen.queryByText(selected.uri)).not.toBeInTheDocument();
  });

  it('fails closed for malformed successful envelopes', async () => {
    listArtifactsMock.mockResolvedValue({
      ...page([]),
      run_id: 'run-other',
    });
    renderBrowser();
    expect(await screen.findByText('Artifacts could not be loaded')).toBeInTheDocument();
    expect(screen.queryByText('No indexed artifacts')).not.toBeInTheDocument();
  });
});
