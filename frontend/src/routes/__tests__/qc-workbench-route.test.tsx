import { act, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import type {
  ArtifactReferenceResponse,
  QcMetricResponse,
} from '../../api/generated/models';
import { listRunQcMetrics } from '../../api/generated/qc-metrics/qc-metrics';
import {
  getRunArtifact,
  listRunArtifacts,
} from '../../api/generated/artifacts/artifacts';
import { createStubRunApiClient, type RunApiClient } from '../../api/runClient';
import { appRoutes } from '../../app/router';
import { renderWithRouter } from '../../test/test-utils';

vi.mock('../../api/generated/qc-metrics/qc-metrics', () => ({
  listRunQcMetrics: vi.fn(),
}));

vi.mock('../../api/generated/artifacts/artifacts', () => ({
  listRunArtifacts: vi.fn(),
  getRunArtifact: vi.fn(),
}));

const listQcMetricsMock = vi.mocked(listRunQcMetrics);
const listArtifactsMock = vi.mocked(listRunArtifacts);
const getArtifactMock = vi.mocked(getRunArtifact);
const METRIC_ID = `qcmetric-${'a'.repeat(64)}`;

const qcMetric: QcMetricResponse = {
  metric_id: METRIC_ID,
  metric_key: 'alignment.mapped_reads',
  display_name: 'Mapped reads',
  value: '9007199254740993.125',
  unit: 'count',
  scope: 'sample',
  sample_id: 'sample-1',
  experiment_id: null,
  assay: 'chipseq',
  qc_flag: 'pass',
  source_artifact_id: 'artifact-qc-summary',
  produced_at: '2026-07-13T08:00:00Z',
};

const sourceArtifact: ArtifactReferenceResponse = {
  artifact_id: 'artifact-qc-summary',
  run_id: 'run-1',
  artifact_type: 'file',
  name: 'qc_summary.tsv',
  uri: 'run://runs/run-1/artifacts/artifact-qc-summary',
  mime_type: 'text/tab-separated-values',
  produced_at: '2026-07-13T08:00:00Z',
  relative_path: 'results/qc/qc_summary.tsv',
  output_type: 'qc_summary',
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
        created_at: '2026-07-13T07:00:00Z',
        updated_at: '2026-07-13T08:00:00Z',
        started_at: '2026-07-13T07:00:10Z',
        ended_at: '2026-07-13T08:00:00Z',
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
          event_id: 'event-artifacts',
          run_id: 'run-1',
          sequence: 1,
          event_type: 'artifacts_indexed',
          timestamp: '2026-07-13T08:00:00Z',
          status: 'succeeded',
          stage: 'artifact_extraction',
          message: 'Artifacts indexed.',
          context: { artifact_count: 1 },
          issue: null,
        },
        {
          event_id: 'event-qc',
          run_id: 'run-1',
          sequence: 2,
          event_type: 'qc_metrics_indexed',
          timestamp: '2026-07-13T08:00:01Z',
          status: 'succeeded',
          stage: 'qc_summary_indexing',
          message: 'QC metrics indexed.',
          context: { metric_count: 1 },
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
  listQcMetricsMock.mockReset();
  listArtifactsMock.mockReset();
  getArtifactMock.mockReset();
  listQcMetricsMock.mockResolvedValue({
    ok: true,
    run_id: 'run-1',
    qc_metrics: [qcMetric],
    next_cursor: null,
    issues: [],
  });
  listArtifactsMock.mockResolvedValue({
    ok: true,
    run_id: 'run-1',
    artifacts: [sourceArtifact],
    next_cursor: null,
    issues: [],
  });
  getArtifactMock.mockResolvedValue({
    ok: true,
    run_id: 'run-1',
    artifact: sourceArtifact,
    issues: [],
  });
});

describe('run QC workbench route state', () => {
  it('restores the QC deep link and browser history', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs/run-1?view=qc'],
      clients: { runClient: succeededRunClient() },
    });

    expect(await screen.findAllByText('Mapped reads')).not.toHaveLength(0);
    expect(screen.getByRole('tab', { name: 'QC' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    expect(screen.getByText('run-1')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Refresh run progress' })).toBeInTheDocument();

    await user.click(screen.getByRole('tab', { name: 'Activity' }));
    await waitFor(() => expect(router.state.location.search).toBe(''));
    await act(async () => router.navigate(-1));
    await waitFor(() => expect(router.state.location.search).toBe('?view=qc'));
    expect(screen.getByRole('tab', { name: 'QC' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    await act(async () => router.navigate(1));
    await waitFor(() => expect(router.state.location.search).toBe(''));
  });

  it('opens a metric source in the existing artifact deep link', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs/run-1?view=qc'],
      clients: { runClient: succeededRunClient() },
    });

    const sourceButtons = await screen.findAllByRole('button', {
      name: 'Open source artifact for Mapped reads',
    });
    await user.click(sourceButtons[0]);
    await waitFor(() =>
      expect(router.state.location.search).toBe(
        '?view=artifacts&artifact=artifact-qc-summary',
      ),
    );
    expect(screen.getByRole('tab', { name: 'Artifacts' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    expect(await screen.findByText(sourceArtifact.uri)).toBeInTheDocument();
  });

  it('clears artifact selection when entering QC and safely falls back for unknown views', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: [
        '/runs/run-1?view=artifacts&artifact=artifact-qc-summary',
      ],
      clients: { runClient: succeededRunClient() },
    });

    expect(await screen.findByText(sourceArtifact.uri)).toBeInTheDocument();
    await user.click(screen.getByRole('tab', { name: 'QC' }));
    await waitFor(() => expect(router.state.location.search).toBe('?view=qc'));

    await act(async () => router.navigate('/runs/run-1?view=unknown'));
    expect(screen.getByRole('tab', { name: 'Activity' })).toHaveAttribute(
      'aria-selected',
      'true',
    );
  });

  it('supports roving keyboard navigation across all workbench tabs', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs/run-1'],
      clients: { runClient: succeededRunClient() },
    });

    const activity = await screen.findByRole('tab', { name: 'Activity' });
    activity.focus();
    await user.keyboard('{ArrowRight}');
    await waitFor(() => expect(router.state.location.search).toBe('?view=artifacts'));
    expect(screen.getByRole('tab', { name: 'Artifacts' })).toHaveFocus();

    await user.keyboard('{ArrowRight}');
    await waitFor(() => expect(router.state.location.search).toBe('?view=qc'));
    expect(screen.getByRole('tab', { name: 'QC' })).toHaveFocus();

    await user.keyboard('{Home}');
    await waitFor(() => expect(router.state.location.search).toBe(''));
    expect(screen.getByRole('tab', { name: 'Activity' })).toHaveFocus();

    await user.keyboard('{End}');
    await waitFor(() => expect(router.state.location.search).toBe('?view=qc'));
    expect(screen.getByRole('tab', { name: 'QC' })).toHaveFocus();
  });
});
