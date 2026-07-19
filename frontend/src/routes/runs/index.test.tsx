import { act, screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import type { RunHistoryResponse, RunSummaryResponse } from '../../api/generated/models';
import { listRuns } from '../../api/generated/runs/runs';
import { ApiError } from '../../api/fetcher';
import { appRoutes } from '../../app/router';
import { renderWithRouter } from '../../test/test-utils';

vi.mock('../../api/generated/runs/runs', async (importOriginal) => {
  const original = await importOriginal<typeof import('../../api/generated/runs/runs')>();
  return { ...original, listRuns: vi.fn() };
});

const listRunsMock = vi.mocked(listRuns);

beforeAll(async () => {
  await import('./index');
});

function run(
  runId: string,
  status: RunSummaryResponse['status'] = 'succeeded',
): RunSummaryResponse {
  return {
    run_id: runId,
    workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
    status,
    created_at: '2026-07-14T08:00:00Z',
    updated_at: '2026-07-14T08:01:00Z',
    started_at: status === 'created' ? null : '2026-07-14T08:00:01Z',
    ended_at: ['succeeded', 'failed', 'cancelled'].includes(status)
      ? '2026-07-14T08:01:00Z'
      : null,
    current_stage: status === 'created' ? null : 'execution',
  };
}

function page(
  runs: RunSummaryResponse[],
  nextCursor: string | null = null,
): RunHistoryResponse {
  return { ok: true, runs, next_cursor: nextCursor, issues: [] };
}

beforeEach(() => {
  listRunsMock.mockReset();
  listRunsMock.mockResolvedValue(page([run('run-1')]));
});

describe('run history route', () => {
  it('loads directly, marks Runs current, and exposes truthful result links', async () => {
    renderWithRouter(appRoutes, { initialEntries: ['/runs'] });

    expect(await screen.findAllByText('run-1')).not.toHaveLength(0);
    expect(screen.getByRole('link', { name: 'Runs' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    expect(screen.getAllByRole('link', { name: 'run-1' })[0]).toHaveAttribute(
      'href',
      '/runs/run-1',
    );
    expect(screen.getAllByRole('link', { name: 'Open QC for run run-1' })[0]).toHaveAttribute(
      'href',
      '/runs/run-1?view=qc',
    );
    expect(
      screen.getAllByRole('link', { name: 'Open artifacts for run run-1' })[0],
    ).toHaveAttribute('href', '/runs/run-1?view=artifacts');
  });

  it('shows stable loading, empty, and initial transport error states', async () => {
    let resolve!: (value: RunHistoryResponse) => void;
    listRunsMock.mockReturnValueOnce(new Promise((done) => { resolve = done; }));
    const loading = renderWithRouter(appRoutes, { initialEntries: ['/runs'] });
    expect(await screen.findByLabelText('Loading run history')).toBeInTheDocument();
    await act(async () => resolve(page([])));
    expect(await screen.findByRole('heading', { name: 'No runs yet' })).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'Browse workflows' })).toHaveAttribute(
      'href',
      '/workflows',
    );
    loading.unmount();

    listRunsMock.mockRejectedValue(new ApiError(503, 'API_ERROR', 'private'));
    renderWithRouter(appRoutes, { initialEntries: ['/runs'] });
    expect(
      await screen.findByRole('heading', { name: 'Run history could not be loaded' }),
    ).toBeInTheDocument();
    expect(screen.queryByText('private')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Retry run history' })).toBeInTheDocument();
  });

  it('uses URL filters, never carries rows across keys, and restores browser history', async () => {
    const user = userEvent.setup();
    listRunsMock.mockImplementation(async (params) => {
      if (params?.status === 'failed') return page([run('run-failed', 'failed')]);
      return page([run('run-success')]);
    });
    const { router } = renderWithRouter(appRoutes, { initialEntries: ['/runs'] });
    expect(await screen.findAllByText('run-success')).not.toHaveLength(0);

    await user.selectOptions(screen.getByLabelText('Filter runs by status'), 'failed');
    await waitFor(() => expect(router.state.location.search).toBe('?status=failed'));
    expect(screen.queryAllByText('run-success')).toHaveLength(0);
    expect(await screen.findAllByText('run-failed')).not.toHaveLength(0);
    expect(listRunsMock).toHaveBeenLastCalledWith(
      expect.objectContaining({ status: 'failed', limit: 50 }),
    );

    await act(async () => router.navigate(-1));
    await waitFor(() => expect(router.state.location.search).toBe(''));
    expect(await screen.findAllByText('run-success')).not.toHaveLength(0);
  });

  it('loads a second keyset page and keeps stable ordering', async () => {
    const user = userEvent.setup();
    listRunsMock
      .mockResolvedValueOnce(page([run('run-b')], 'runhist_next'))
      .mockResolvedValueOnce(page([run('run-a')], null));
    renderWithRouter(appRoutes, { initialEntries: ['/runs'] });

    expect(await screen.findAllByText('run-b')).not.toHaveLength(0);
    await user.click(screen.getByRole('button', { name: 'Load more runs' }));
    expect(await screen.findAllByText('run-a')).not.toHaveLength(0);
    expect(listRunsMock).toHaveBeenLastCalledWith(
      expect.objectContaining({ after: 'runhist_next', limit: 50 }),
    );
    const desktopRows = within(screen.getByRole('table', { name: 'Run history' }))
      .getAllByRole('row');
    expect(desktopRows[1]).toHaveTextContent('run-b');
    expect(desktopRows[2]).toHaveTextContent('run-a');
  });

  it('retains cached rows on refresh failure and offers transport retry', async () => {
    const user = userEvent.setup();
    listRunsMock.mockResolvedValueOnce(page([run('run-cached')]));
    renderWithRouter(appRoutes, { initialEntries: ['/runs'] });
    expect(await screen.findAllByText('run-cached')).not.toHaveLength(0);
    listRunsMock.mockRejectedValue(new ApiError(503, 'API_ERROR', 'private'));

    await user.click(screen.getByRole('button', { name: 'Refresh run history' }));

    expect(await screen.findByRole('alert')).toHaveTextContent(
      'Previously loaded canonical rows are retained',
    );
    expect(screen.getAllByText('run-cached')).not.toHaveLength(0);
    expect(screen.getByRole('button', { name: 'Retry' })).toBeInTheDocument();
  });

  it('blocks invalid or repeated URL filters without issuing a request', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs?status=failed&status=succeeded'],
    });

    expect(
      await screen.findByRole('heading', { name: 'Run filters could not be used' }),
    ).toBeInTheDocument();
    expect(listRunsMock).not.toHaveBeenCalled();
    await user.click(screen.getByRole('button', { name: 'Reset filters' }));
    await waitFor(() => expect(router.state.location.search).toBe(''));
    expect(await screen.findAllByText('run-1')).not.toHaveLength(0);
  });

  it('treats a repeated response cursor as protocol recovery, not empty history', async () => {
    const user = userEvent.setup();
    listRunsMock
      .mockResolvedValueOnce(page([run('run-first')], 'runhist_repeat'))
      .mockResolvedValue(page([run('run-second')], 'runhist_repeat'));
    renderWithRouter(appRoutes, { initialEntries: ['/runs'] });
    expect(await screen.findAllByText('run-first')).not.toHaveLength(0);

    await user.click(screen.getByRole('button', { name: 'Load more runs' }));

    expect(await screen.findByRole('alert')).toHaveTextContent(
      'Loaded rows are retained',
    );
    expect(screen.getAllByText('run-first')).not.toHaveLength(0);
    expect(screen.getByRole('button', { name: 'Reload from first page' })).toBeInTheDocument();
  });
});
