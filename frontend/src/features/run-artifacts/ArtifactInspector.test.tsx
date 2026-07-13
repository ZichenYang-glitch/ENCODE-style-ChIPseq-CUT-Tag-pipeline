import { afterEach, describe, expect, it, vi } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { ArtifactReferenceResponse } from '../../api/generated/models';
import { ArtifactInspector } from './ArtifactInspector';

const artifact: ArtifactReferenceResponse = {
  artifact_id: 'artifact-a',
  run_id: 'run-1',
  artifact_type: 'file',
  name: 'result_manifest.tsv',
  uri: 'run://runs/run-1/artifacts/artifact-a',
  mime_type: 'text/tab-separated-values',
  produced_at: '2026-07-12T12:00:00Z',
  relative_path: 'results/multiqc/result_manifest.tsv',
  output_type: 'result_manifest',
  size_bytes: 42,
  metadata: {
    catalog_id: 'result_manifest',
    scope: 'project',
    sample_id: 'sample-1',
    experiment_id: 'experiment-1',
    assay: 'ChIP-seq',
    target: 'H3K27ac',
    genome: 'hg38',
    method: 'manifest',
    qc_flag: 'pass',
  },
};

const originalScrollIntoView = Object.getOwnPropertyDescriptor(
  HTMLElement.prototype,
  'scrollIntoView',
);

afterEach(() => {
  vi.unstubAllGlobals();
  if (originalScrollIntoView) {
    Object.defineProperty(
      HTMLElement.prototype,
      'scrollIntoView',
      originalScrollIntoView,
    );
  } else {
    Reflect.deleteProperty(HTMLElement.prototype, 'scrollIntoView');
  }
});

describe('ArtifactInspector', () => {
  it('shows persisted public detail and a controlled download action', async () => {
    const user = userEvent.setup();
    const onDownload = vi.fn();
    render(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError={false}
        invalidSelection={false}
        onRetry={() => undefined}
        onDownload={onDownload}
      />,
    );

    expect(screen.getAllByText('result_manifest')).toHaveLength(2);
    expect(screen.getByText('results/multiqc/result_manifest.tsv')).toBeInTheDocument();
    expect(screen.getByText('run://runs/run-1/artifacts/artifact-a')).toBeInTheDocument();
    expect(screen.getByText('sample-1')).toBeInTheDocument();
    expect(screen.getByText('experiment-1')).toBeInTheDocument();
    await user.click(
      screen.getByRole('button', { name: 'Download result_manifest.tsv' }),
    );
    expect(onDownload).toHaveBeenCalledWith(artifact);
  });

  it('shows honest download progress, success, and redacted failure', () => {
    const { rerender } = render(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError={false}
        invalidSelection={false}
        onRetry={() => undefined}
        onDownload={() => undefined}
        isDownloading
      />,
    );
    expect(screen.getByRole('button', { name: 'Download result_manifest.tsv' })).toBeDisabled();
    expect(screen.getByText('Downloading…')).toBeInTheDocument();

    rerender(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError={false}
        invalidSelection={false}
        onRetry={() => undefined}
        onDownload={() => undefined}
        downloadStatus="success"
      />,
    );
    expect(screen.getByText('Download prepared successfully.')).toHaveAttribute(
      'role',
      'status',
    );

    rerender(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError={false}
        invalidSelection={false}
        onRetry={() => undefined}
        onDownload={() => undefined}
        downloadStatus="error"
      />,
    );
    expect(screen.getByText(/Download could not be completed/)).toHaveAttribute(
      'role',
      'status',
    );
    expect(screen.queryByText(/private|exception|stderr/i)).not.toBeInTheDocument();
  });

  it('announces copy success and failure without changing the value', async () => {
    const user = userEvent.setup();
    const writeText = vi.fn().mockResolvedValueOnce(undefined).mockRejectedValueOnce(new Error('denied'));
    vi.stubGlobal('navigator', { ...navigator, clipboard: { writeText } });
    render(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError={false}
        invalidSelection={false}
        onRetry={() => undefined}
      />,
    );

    await user.click(screen.getByRole('button', { name: 'Copy relative path' }));
    expect(writeText).toHaveBeenCalledWith(artifact.relative_path);
    expect(screen.getByRole('status')).toHaveTextContent('Relative path copied.');
    await user.click(screen.getByRole('button', { name: 'Copy opaque uri' }));
    expect(screen.getByRole('status')).toHaveTextContent('Could not copy opaque uri.');
  });

  it('renders stable loading, invalid-link, and retry states', async () => {
    const user = userEvent.setup();
    const onRetry = vi.fn();
    const { rerender } = render(
      <ArtifactInspector
        artifact={null}
        selectedArtifactId="artifact-a"
        isLoading
        isError={false}
        invalidSelection={false}
        onRetry={onRetry}
      />,
    );
    expect(screen.getByRole('status', { name: 'Loading artifact details' })).toBeInTheDocument();

    rerender(
      <ArtifactInspector
        artifact={null}
        selectedArtifactId="artifact-a"
        isLoading={false}
        isError
        invalidSelection={false}
        onRetry={onRetry}
      />,
    );
    await user.click(screen.getByRole('button', { name: 'Retry details' }));
    expect(onRetry).toHaveBeenCalledTimes(1);

    rerender(
      <ArtifactInspector
        artifact={null}
        selectedArtifactId="../escape"
        isLoading={false}
        isError={false}
        invalidSelection
        onRetry={onRetry}
      />,
    );
    expect(screen.getByText('This artifact link is invalid or unavailable.')).toBeInTheDocument();
  });

  it('preserves confirmed detail when a refresh fails', async () => {
    const user = userEvent.setup();
    const onRetry = vi.fn();
    render(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError
        invalidSelection={false}
        onRetry={onRetry}
      />,
    );

    expect(screen.getByText(artifact.uri)).toBeInTheDocument();
    expect(
      screen.getByText('Artifact detail refresh failed. The last confirmed detail is preserved.'),
    ).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Retry details' }));
    expect(onRetry).toHaveBeenCalledTimes(1);
  });

  it('reveals a selected detail on narrow screens and returns to the list', async () => {
    const user = userEvent.setup();
    const scrollIntoView = vi.fn();
    const onBackToList = vi.fn();
    Object.defineProperty(HTMLElement.prototype, 'scrollIntoView', {
      configurable: true,
      value: scrollIntoView,
    });
    vi.stubGlobal('innerWidth', 390);

    render(
      <ArtifactInspector
        artifact={artifact}
        selectedArtifactId={artifact.artifact_id}
        isLoading={false}
        isError={false}
        invalidSelection={false}
        onRetry={() => undefined}
        onBackToList={onBackToList}
      />,
    );

    await waitFor(() =>
      expect(scrollIntoView).toHaveBeenCalledWith({ block: 'start' }),
    );
    await user.click(screen.getByRole('button', { name: 'Back to artifact list' }));
    expect(onBackToList).toHaveBeenCalledTimes(1);
  });
});
