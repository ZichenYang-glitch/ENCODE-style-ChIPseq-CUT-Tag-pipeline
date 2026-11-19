import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen } from '@testing-library/react';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import { MemoryRouter } from 'react-router-dom';
import type { PublishedResultContextResponse } from '../../api/generated/models';
import { readPublishedResultContext } from '../../api/generated/result-context/result-context';
import { PublicationResultContext } from './PublicationResultContext';

vi.mock('../../api/generated/result-context/result-context', () => ({
  readPublishedResultContext: vi.fn(),
}));

const readContext = vi.mocked(readPublishedResultContext);
const PROJECT_ID = `prj_${'a'.repeat(32)}`;
const REVISION_ID = `expr_${'b'.repeat(32)}`;
const DIGEST = 'c'.repeat(64);
const GENERATION = `artifactgen-${'d'.repeat(64)}`;

function renderContext() {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  return render(
    <QueryClientProvider client={queryClient}>
      <MemoryRouter>
        <PublicationResultContext
          projectId={PROJECT_ID}
          runId="run-1"
          artifactId="artifact-1"
          generation={GENERATION}
        />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function context(
  overrides: Partial<PublishedResultContextResponse> = {},
): PublishedResultContextResponse {
  return {
    ok: true,
    project_id: PROJECT_ID,
    run_id: 'run-1',
    artifact_id: 'artifact-1',
    artifact_generation: GENERATION,
    run_status: 'succeeded',
    reference_evidence: 'present',
    qc_evidence: 'present',
    preflight_failure_evidence: 'present',
    experiment_trace: {
      ok: true,
      project_id: PROJECT_ID,
      run_id: 'run-1',
      artifact_id: 'artifact-1',
      artifact_generation: GENERATION,
      trace_status: 'experiment_snapshot',
      experiment_revision_id: REVISION_ID,
      payload_digest: DIGEST,
      result: {} as PublishedResultContextResponse['experiment_trace']['result'],
    },
    reference_profile: {
      profile_id: `refp_${'e'.repeat(32)}`,
      revision_id: `refpr_${'f'.repeat(32)}`,
      revision_number: 1,
      display_name: 'hg38 tiny',
      organism: 'human',
      assembly: 'hg38',
      identity_sha256: 'a'.repeat(64),
    },
    qc_metrics: [
      {
        metric_id: `qcmetric-${'1'.repeat(64)}`,
        metric_key: 'alignment.mapped_fraction',
        display_name: 'Mapped fraction',
        value: '0.5',
        unit: 'fraction',
        scope: 'run',
        sample_id: null,
        experiment_id: null,
        assay: null,
        qc_flag: null,
        source_artifact_id: 'artifact-1',
        produced_at: '2026-11-19T00:00:00Z',
      },
    ],
    preflight_failures: [
      {
        event_id: 'event-1',
        event_type: 'preflight_failed',
        message: 'Reference profile could not be prepared.',
        stage: 'preflight',
        code: 'PREFLIGHT_FAILED',
      },
    ],
    ...overrides,
  };
}

beforeEach(() => {
  readContext.mockReset();
});

describe('PublicationResultContext', () => {
  it('shows the copied revision and the recorded evidence', async () => {
    readContext.mockResolvedValue(context());
    renderContext();

    expect(await screen.findByTitle(REVISION_ID)).toBeInTheDocument();
    expect(screen.getByTitle(DIGEST)).toBeInTheDocument();
    expect(screen.getByText(/Run status succeeded/)).toBeInTheDocument();
    expect(screen.getByText(/hg38 tiny/)).toBeInTheDocument();
    expect(screen.getByText(/Mapped fraction: 0.5 fraction/)).toBeInTheDocument();
    expect(screen.getByText(/PREFLIGHT_FAILED/)).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'Open run' })).toHaveAttribute(
      'href',
      '/runs/run-1',
    );
    expect(readContext).toHaveBeenCalledWith(PROJECT_ID, 'run-1', 'artifact-1', {
      generation: GENERATION,
    });
  });

  it('shows why an unbound result cannot name a revision', async () => {
    readContext.mockResolvedValue(
      context({
        reference_evidence: 'absent',
        reference_profile: undefined,
        qc_evidence: 'absent',
        qc_metrics: undefined,
        preflight_failure_evidence: 'absent',
        preflight_failures: undefined,
        experiment_trace: {
          ok: true,
          project_id: PROJECT_ID,
          run_id: 'run-1',
          artifact_id: 'artifact-1',
          artifact_generation: GENERATION,
          trace_status: 'unresolved',
          reason: 'experiment_snapshot_binding_absent',
          result: {} as PublishedResultContextResponse['experiment_trace']['result'],
        },
      }),
    );
    renderContext();

    expect(
      await screen.findByText('This run has no copied experiment revision.'),
    ).toBeInTheDocument();
    expect(
      screen.getByText('No reference profile was recorded for this run.'),
    ).toBeInTheDocument();
    expect(
      screen.getByText('No QC metrics were recorded for this run.'),
    ).toBeInTheDocument();
    expect(screen.getByText('No preflight failure was recorded.')).toBeInTheDocument();
    expect(screen.queryByText(REVISION_ID)).not.toBeInTheDocument();
  });
});
