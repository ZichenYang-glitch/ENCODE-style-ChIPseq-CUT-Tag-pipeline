import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { ArtifactReferenceResponse } from '../../api/generated/models';
import { ArtifactList } from './ArtifactList';

const longPath = `results/${'nested/'.repeat(15)}result_manifest.tsv`;
const artifact: ArtifactReferenceResponse = {
  artifact_id: 'result_manifest.abc',
  run_id: 'run-1',
  artifact_type: 'file',
  name: 'result_manifest.tsv',
  uri: 'run://runs/run-1/artifacts/result_manifest.abc',
  mime_type: 'text/tab-separated-values',
  produced_at: '2026-07-12T12:00:00Z',
  revision: `artifactrev-${'a'.repeat(64)}`,
  relative_path: longPath,
  output_type: 'result_manifest',
  size_bytes: 1536,
  metadata: {
    scope: 'project',
    sample_id: 'sample-1',
    experiment_id: 'experiment-1',
    assay: 'ChIP-seq',
  },
};

describe('ArtifactList', () => {
  it('renders required fields in desktop and mobile presentations', () => {
    render(
      <ArtifactList
        artifacts={[artifact]}
        selectedArtifactId={null}
        onSelect={() => undefined}
        hasNextPage={false}
        isFetchingNextPage={false}
        onLoadMore={() => undefined}
      />,
    );

    expect(screen.getAllByText('result_manifest')).toHaveLength(2);
    expect(screen.getAllByText('result_manifest.tsv')).toHaveLength(2);
    expect(screen.getAllByText(longPath)).toHaveLength(2);
    expect(screen.getAllByText('1.5 KiB')).toHaveLength(2);
    expect(screen.getAllByText('project')).toHaveLength(2);
    expect(screen.getAllByText('sample-1 · experiment-1')).toHaveLength(2);
    expect(screen.getAllByText('ChIP-seq')).toHaveLength(2);
    expect(screen.getByRole('table', { name: 'Indexed run artifacts' })).toBeInTheDocument();
  });

  it('selects an artifact accessibly and loads another page', async () => {
    const user = userEvent.setup();
    const onSelect = vi.fn();
    const onLoadMore = vi.fn();
    render(
      <ArtifactList
        artifacts={[artifact]}
        selectedArtifactId={artifact.artifact_id}
        onSelect={onSelect}
        hasNextPage
        isFetchingNextPage={false}
        onLoadMore={onLoadMore}
      />,
    );

    const openButtons = screen.getAllByRole('button', {
      name: 'Open artifact result_manifest.tsv',
    });
    expect(openButtons.every((button) => button.getAttribute('aria-pressed') === 'true')).toBe(true);
    await user.click(openButtons[0]);
    expect(onSelect).toHaveBeenCalledWith(artifact.artifact_id);
    await user.click(screen.getByRole('button', { name: 'Load more artifacts' }));
    expect(onLoadMore).toHaveBeenCalledTimes(1);
  });

  it('uses a stable pending Load more control', () => {
    render(
      <ArtifactList
        artifacts={[artifact]}
        selectedArtifactId={null}
        onSelect={() => undefined}
        hasNextPage
        isFetchingNextPage
        onLoadMore={() => undefined}
      />,
    );
    expect(screen.getByRole('button', { name: 'Load more artifacts' })).toBeDisabled();
    expect(screen.getByText('Loading more…')).toBeInTheDocument();
  });
});
