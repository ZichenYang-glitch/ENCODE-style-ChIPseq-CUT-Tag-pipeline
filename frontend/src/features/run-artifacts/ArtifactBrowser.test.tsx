import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import type { ComponentProps } from 'react';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type {
  ArtifactReferenceResponse,
  RunArtifactDetailResponse,
  RunArtifactsResponse,
} from '../../api/generated/models';
import {
  getRunArtifact,
  listRunArtifacts,
} from '../../api/generated/artifacts/artifacts';
import { ArtifactBrowser } from './ArtifactBrowser';

vi.mock('../../api/generated/artifacts/artifacts', () => ({
  listRunArtifacts: vi.fn(),
  getRunArtifact: vi.fn(),
}));

const listArtifactsMock = vi.mocked(listRunArtifacts);
const getArtifactMock = vi.mocked(getRunArtifact);

function artifact(id: string): ArtifactReferenceResponse {
  return {
    artifact_id: id,
    run_id: 'run-1',
    artifact_type: 'file',
    name: `${id}.tsv`,
    uri: `run://runs/run-1/artifacts/${id}`,
    mime_type: 'text/tab-separated-values',
    produced_at: '2026-07-12T12:00:00Z',
    relative_path: `results/${id}.tsv`,
    output_type: id,
    size_bytes: 42,
    metadata: { scope: 'project' },
  };
}

function page(
  artifacts: ArtifactReferenceResponse[],
  nextCursor: string | null = null,
): RunArtifactsResponse {
  return {
    ok: true,
    run_id: 'run-1',
    artifacts,
    next_cursor: nextCursor,
    issues: [],
  };
}

function detail(value: ArtifactReferenceResponse): RunArtifactDetailResponse {
  return { ok: true, run_id: 'run-1', artifact: value, issues: [] };
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
    outcome: { kind: 'indexed', count: 1 },
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
    renderBrowser({ outcome: { kind: 'indexed', count: 0 } });
    expect(await screen.findByText('No indexed artifacts')).toBeInTheDocument();
  });

  it('keeps nonzero empty responses honest', async () => {
    listArtifactsMock.mockResolvedValue(page([]));
    renderBrowser({ outcome: { kind: 'indexed', count: 1 } });
    expect(await screen.findByText('Artifact index is not ready')).toBeInTheDocument();
    expect(screen.queryByText('No indexed artifacts')).not.toBeInTheDocument();
  });
});

describe('ArtifactBrowser queries', () => {
  it('loads explicit pages of 50 and stops a repeated cursor', async () => {
    const user = userEvent.setup();
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], 'artifact-a'))
      .mockResolvedValueOnce(page([artifact('artifact-b')], 'artifact-a'));
    renderBrowser({ outcome: { kind: 'indexed', count: 2 } });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));
    expect(await screen.findAllByText('artifact-b')).not.toHaveLength(0);
    expect(listArtifactsMock).toHaveBeenNthCalledWith(
      1,
      'run-1',
      { after: undefined, limit: 50 },
    );
    expect(listArtifactsMock).toHaveBeenNthCalledWith(
      2,
      'run-1',
      { after: 'artifact-a', limit: 50 },
    );
    expect(screen.queryByRole('button', { name: 'Load more artifacts' })).not.toBeInTheDocument();
  });

  it('preserves existing rows when another page fails and offers retry', async () => {
    const user = userEvent.setup();
    listArtifactsMock
      .mockResolvedValueOnce(page([artifact('artifact-a')], 'artifact-a'))
      .mockRejectedValueOnce(new Error('/private/path'));
    renderBrowser({ outcome: { kind: 'indexed', count: 2 } });

    expect(await screen.findAllByText('artifact-a')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));
    expect(await screen.findByText(/Existing rows are preserved/)).toBeInTheDocument();
    expect(screen.getAllByText('artifact-a').length).toBeGreaterThan(0);
    expect(screen.queryByText('/private/path')).not.toBeInTheDocument();
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
    expect(getArtifactMock).toHaveBeenCalledWith('run-1', 'artifact-a');
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
      queryKey: ['run-artifact', 'run-1', 'artifact-a'],
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
