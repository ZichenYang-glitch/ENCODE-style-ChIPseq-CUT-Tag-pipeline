import { beforeEach, describe, expect, it, vi } from 'vitest';
import { act, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { appRoutes } from '../../app/router';
import { renderWithRouter } from '../../test/test-utils';
import { createStubRunApiClient } from '../../api/runClient';
import type { RunApiClient } from '../../api/runClient';
import type { ArtifactReferenceResponse } from '../../api/generated/models';
import {
  getRunArtifact,
  listRunArtifacts,
} from '../../api/generated/artifacts/artifacts';

vi.mock('../../api/generated/artifacts/artifacts', () => ({
  listRunArtifacts: vi.fn(),
  getRunArtifact: vi.fn(),
}));

const listArtifactsMock = vi.mocked(listRunArtifacts);
const getArtifactMock = vi.mocked(getRunArtifact);

const artifact: ArtifactReferenceResponse = {
  artifact_id: 'artifact-a',
  run_id: 'run-1',
  artifact_type: 'file',
  name: 'result_manifest.tsv',
  uri: 'run://runs/run-1/artifacts/artifact-a',
  mime_type: 'text/tab-separated-values',
  produced_at: '2026-07-12T12:01:00Z',
  revision: `artifactrev-${'a'.repeat(64)}`,
  relative_path: 'results/multiqc/result_manifest.tsv',
  output_type: 'result_manifest',
  size_bytes: 42,
  metadata: { scope: 'project' },
};

function succeededRunClient(): RunApiClient {
  const base = createStubRunApiClient();
  return {
    ...base,
    getRun: vi.fn(async () => ({
      ok: true,
      run: {
        run_id: 'run-1',
        workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
        inputs: {},
        status: 'succeeded',
        created_at: '2026-07-12T12:00:00Z',
        updated_at: '2026-07-12T12:01:00Z',
        started_at: '2026-07-12T12:00:10Z',
        ended_at: '2026-07-12T12:01:00Z',
        current_stage: 'execution',
        cancellation_reason: null,
        error: null,
        tags: {},
      },
      issues: [],
    })),
    listRunEvents: vi.fn(async () => ({
      ok: true,
      run_id: 'run-1',
      events: [
        {
          event_id: 'event-indexed',
          run_id: 'run-1',
          sequence: 1,
          event_type: 'artifacts_indexed',
          timestamp: '2026-07-12T12:01:00Z',
          status: 'succeeded',
          stage: 'artifact_extraction',
          message: 'Artifacts indexed.',
          context: {
            artifact_count: 1,
            artifact_generation: `artifactgen-${'a'.repeat(64)}`,
          },
          issue: null,
        },
      ],
      next_cursor: null,
      issues: [],
    })),
    listRunLogs: vi.fn(async (_runId, options = {}) => ({
      ok: true,
      run_id: 'run-1',
      stream_name: options.streamName ?? 'stdout',
      chunks: [],
      next_cursor: null,
      issues: [],
    })),
  };
}

beforeEach(() => {
  listArtifactsMock.mockReset();
  getArtifactMock.mockReset();
  listArtifactsMock.mockResolvedValue({
    ok: true,
    run_id: 'run-1',
    artifact_generation: `artifactgen-${'a'.repeat(64)}`,
    artifacts: [artifact],
    next_cursor: null,
    issues: [],
  });
  getArtifactMock.mockResolvedValue({
    ok: true,
    run_id: 'run-1',
    artifact_generation: `artifactgen-${'a'.repeat(64)}`,
    artifact,
    issues: [],
  });
});

describe('run artifact browser route state', () => {
  it('stores tabs and selection in browser history', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs/run-1'],
      clients: { runClient: succeededRunClient() },
    });

    expect(await screen.findByRole('tab', { name: 'Activity' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    expect(screen.getByText('run-1')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Refresh run progress' })).toBeInTheDocument();

    await user.click(screen.getByRole('tab', { name: 'Artifacts' }));
    await waitFor(() => expect(router.state.location.search).toBe('?view=artifacts'));
    const openButtons = await screen.findAllByRole('button', {
      name: 'Open artifact result_manifest.tsv',
    });
    await user.click(openButtons[0]);
    await waitFor(() =>
      expect(router.state.location.search).toBe(
        '?view=artifacts&artifact=artifact-a',
      ),
    );
    expect(await screen.findByText(artifact.uri)).toBeInTheDocument();

    await act(async () => router.navigate(-1));
    await waitFor(() => expect(router.state.location.search).toBe('?view=artifacts'));
    await act(async () => router.navigate(-1));
    await waitFor(() => expect(router.state.location.search).toBe(''));
    expect(screen.getByRole('tab', { name: 'Activity' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    await act(async () => router.navigate(1));
    await waitFor(() => expect(router.state.location.search).toBe('?view=artifacts'));
  });

  it('restores a detail deep link and lets artifact imply the view', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/runs/run-1?artifact=artifact-a'],
      clients: { runClient: succeededRunClient() },
    });

    expect(await screen.findByText(artifact.uri)).toBeInTheDocument();
    expect(screen.getByRole('tab', { name: 'Artifacts' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    expect(getArtifactMock).toHaveBeenCalledWith('run-1', 'artifact-a', {
      generation: `artifactgen-${'a'.repeat(64)}`,
    });
  });

  it('rejects an unsafe decoded artifact deep link before the generated call', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/runs/run-1?view=artifacts&artifact=..%2Fescape'],
      clients: { runClient: succeededRunClient() },
    });

    expect(
      await screen.findByText('This artifact link is invalid or unavailable.'),
    ).toBeInTheDocument();
    expect(getArtifactMock).not.toHaveBeenCalled();
  });
});
