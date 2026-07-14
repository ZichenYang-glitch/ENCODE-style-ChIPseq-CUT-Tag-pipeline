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
      { path: 'runs', element: <p>Run history</p> },
      { path: 'runs/:runId', element: <p>Run detail</p> },
      { path: '*', element: <p>Unknown route</p> },
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
      screen.getByRole('link', { name: 'HelixWeave' }),
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

  it('marks Runs current on history and run detail while keeping Workflows available', () => {
    const { unmount } = renderWithRouter(routes, {
      initialEntries: ['/runs'],
    });
    expect(screen.getByRole('link', { name: 'Runs' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    expect(screen.getByRole('link', { name: 'Workflows' })).not.toHaveAttribute(
      'aria-current',
    );
    unmount();

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
    expect(screen.getByRole('link', { name: 'Runs' })).toHaveAttribute(
      'href',
      '/runs',
    );
    expect(screen.getByRole('link', { name: 'Runs' })).toHaveAttribute(
      'aria-current',
      'page',
    );
  });

  it('does not mark Runs current on unknown run-shaped routes', () => {
    renderWithRouter(routes, {
      initialEntries: ['/runs/run-123/unknown'],
    });

    expect(screen.getByRole('link', { name: 'Runs' })).not.toHaveAttribute(
      'aria-current',
    );
  });
});
