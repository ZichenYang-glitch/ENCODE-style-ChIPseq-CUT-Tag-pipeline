import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { RouteObject } from 'react-router-dom';
import { AppShell } from './AppShell';
import { renderWithRouter } from '../test/test-utils';

const routes: RouteObject[] = [
  {
    path: '/',
    element: <AppShell />,
    children: [
      { path: 'workflows', element: <p>Workflow catalog</p> },
      {
        path: 'workflows/:workflowId/new-run',
        element: <p>Input workbench</p>,
      },
      {
        path: 'workflows/:workflowId',
        element: <p>Workflow detail</p>,
      },
      { path: 'runs/:runId', element: <p>Run detail</p> },
    ],
  },
];

describe('AppShell navigation', () => {
  it('uses the brand as a link back to the workflow catalog', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(routes, {
      initialEntries: ['/runs/run-123'],
    });

    await user.click(
      screen.getByRole('link', { name: 'Workflow Platform' }),
    );

    expect(router.state.location.pathname).toBe('/workflows');
  });

  it('marks Workflows current on the catalog and exact workflow detail', () => {
    const { unmount } = renderWithRouter(routes, {
      initialEntries: ['/workflows'],
    });

    expect(
      screen.getByRole('navigation', { name: 'Primary' }),
    ).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'Workflows' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    unmount();

    renderWithRouter(routes, {
      initialEntries: ['/workflows/rna-seq'],
    });
    expect(screen.getByRole('link', { name: 'Workflows' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    expect(screen.getByRole('link', { name: 'New analysis' })).toHaveAttribute(
      'href',
      '/workflows/rna-seq/new-run',
    );
  });

  it('keeps Workflows current for the equivalent trailing-slash catalog URL', () => {
    renderWithRouter(routes, {
      initialEntries: ['/workflows/'],
    });

    expect(screen.getByRole('link', { name: 'Workflows' })).toHaveAttribute(
      'aria-current',
      'page',
    );
  });

  it('derives a workflow-neutral New analysis destination from the URL', () => {
    renderWithRouter(routes, {
      initialEntries: ['/workflows/rna-seq/new-run'],
    });

    const authoringLink = screen.getByRole('link', { name: 'New analysis' });
    expect(authoringLink).toHaveAttribute('href', '/workflows/rna-seq/new-run');
    expect(authoringLink).toHaveAttribute('aria-current', 'page');
    expect(screen.getByRole('link', { name: 'Workflows' })).not.toHaveAttribute(
      'aria-current',
    );
  });

  it('keeps Workflows available on run detail without inventing run navigation', () => {
    renderWithRouter(routes, {
      initialEntries: ['/runs/run-123'],
    });

    expect(screen.getByRole('link', { name: 'Workflows' })).toHaveAttribute(
      'href',
      '/workflows',
    );
    expect(screen.getByRole('link', { name: 'Workflows' })).not.toHaveAttribute(
      'aria-current',
    );
    expect(
      screen.queryByRole('link', { name: 'New analysis' }),
    ).not.toBeInTheDocument();
    expect(screen.queryByRole('link', { name: 'Runs' })).not.toBeInTheDocument();
  });
});
