import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import type { QcMetricResponse } from '../../api/generated/models';
import { QcMetricList } from './QcMetricList';

const METRIC_ID = `qcmetric-${'a'.repeat(64)}`;

function metric(
  overrides: Partial<QcMetricResponse> = {},
): QcMetricResponse {
  return {
    metric_id: METRIC_ID,
    metric_key: 'alignment.mapped_reads',
    display_name: 'Mapped reads',
    value: '9007199254740993.125',
    unit: 'count',
    scope: 'sample',
    sample_id: '-S1',
    experiment_id: 'EXP-1',
    assay: 'chipseq',
    qc_flag: 'pass',
    source_artifact_id: 'artifact-qc-summary',
    produced_at: '2026-07-13T08:00:00Z',
    ...overrides,
  };
}

describe('QcMetricList', () => {
  it('renders exact persisted facts in semantic desktop and mobile layouts', () => {
    render(
      <QcMetricList
        metrics={[metric()]}
        hasNextPage={false}
        isFetchingNextPage={false}
        onLoadMore={vi.fn()}
        onOpenSourceArtifact={vi.fn()}
      />,
    );

    expect(screen.getByRole('table', { name: 'Indexed run QC metrics' })).toBeInTheDocument();
    expect(screen.getByRole('list', { name: 'Indexed run QC metrics' })).toBeInTheDocument();
    expect(screen.getAllByText('Mapped reads').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('alignment.mapped_reads').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('9007199254740993.125').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('count').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('sample').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('-S1 · EXP-1').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('chipseq').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('pass').length).toBeGreaterThanOrEqual(2);
    expect(screen.queryByText(/good|bad|threshold/i)).not.toBeInTheDocument();
    expect(document.querySelectorAll('time[datetime="2026-07-13T08:00:00Z"]')).toHaveLength(2);
  });

  it('navigates to the persisted source artifact with an accessible icon action', async () => {
    const user = userEvent.setup();
    const onOpenSourceArtifact = vi.fn();
    render(
      <QcMetricList
        metrics={[metric()]}
        hasNextPage={false}
        isFetchingNextPage={false}
        onLoadMore={vi.fn()}
        onOpenSourceArtifact={onOpenSourceArtifact}
      />,
    );

    const sourceButtons = screen.getAllByRole('button', {
      name: 'Open source artifact for Mapped reads',
    });
    expect(sourceButtons[0]).toHaveAttribute('title', 'Open source artifact');
    await user.click(sourceButtons[0]);
    expect(onOpenSourceArtifact).toHaveBeenCalledWith('artifact-qc-summary');
  });

  it('uses explicit null labels and protects long values from layout expansion', () => {
    const longKey = `alignment.${'metric'.repeat(30)}`;
    render(
      <QcMetricList
        metrics={[
          metric({
            metric_key: longKey,
            sample_id: null,
            experiment_id: null,
            assay: null,
            qc_flag: null,
          }),
        ]}
        hasNextPage={false}
        isFetchingNextPage={false}
        onLoadMore={vi.fn()}
        onOpenSourceArtifact={vi.fn()}
      />,
    );

    expect(screen.getAllByText('Not reported').length).toBeGreaterThanOrEqual(2);
    expect(screen.getAllByText('—').length).toBeGreaterThanOrEqual(4);
    expect(screen.getAllByTitle(longKey)).toHaveLength(2);
    expect(screen.getByTestId('qc-metric-list')).toHaveClass('min-w-0');
  });

  it('loads more explicitly and disables the control while pending', async () => {
    const user = userEvent.setup();
    const onLoadMore = vi.fn();
    const { rerender } = render(
      <QcMetricList
        metrics={[metric()]}
        hasNextPage
        isFetchingNextPage={false}
        onLoadMore={onLoadMore}
        onOpenSourceArtifact={vi.fn()}
      />,
    );

    await user.click(screen.getByRole('button', { name: 'Load more QC metrics' }));
    expect(onLoadMore).toHaveBeenCalledTimes(1);

    rerender(
      <QcMetricList
        metrics={[metric()]}
        hasNextPage
        isFetchingNextPage
        onLoadMore={onLoadMore}
        onOpenSourceArtifact={vi.fn()}
      />,
    );
    expect(screen.getByRole('button', { name: 'Load more QC metrics' })).toBeDisabled();
  });
});
