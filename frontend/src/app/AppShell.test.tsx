import { describe, expect, it, vi } from 'vitest';
import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { RouteObject } from 'react-router-dom';
import { AppShell } from './AppShell';
import { renderWithRouter } from '../test/test-utils';
import { fetchSessionState } from '../api/authClient';
import {
  getTerminalEmailPreference,
  setTerminalEmailPreference,
} from '../api/generated/auth/auth';

vi.mock('../api/generated/auth/auth', () => ({
  getTerminalEmailPreference: vi.fn(),
  setTerminalEmailPreference: vi.fn(),
}));

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
      { path: 'artifacts', element: <p>Artifact publications</p> },
      {
        path: 'artifacts/:runId/:artifactId',
        element: <p>Artifact publication detail</p>,
      },
      { path: '*', element: <p>Unknown route</p> },
    ],
  },
];

describe('AppShell navigation', () => {
  it('lets a member opt out without exposing or editing the address', async () => {
    vi.mocked(fetchSessionState).mockResolvedValueOnce({
      ok: true,
      setup_required: false,
      authenticated: true,
      principal: {
        user_id: 'usr_11111111111111111111111111111111',
        username: 'lab-member',
        role: 'member',
      },
    });
    vi.mocked(getTerminalEmailPreference).mockResolvedValueOnce({
      terminal_email_enabled: true,
      address_configured: true,
    });
    vi.mocked(setTerminalEmailPreference).mockResolvedValueOnce({
      terminal_email_enabled: false,
      address_configured: true,
    });
    const user = userEvent.setup();

    renderWithRouter(routes, { initialEntries: ['/runs'] });

    await user.click(await screen.findByRole('button', { name: /lab-member/ }));
    const preference = await screen.findByRole('checkbox', {
      name: 'Email me when my runs finish',
    });
    expect(preference).toBeChecked();
    expect(screen.queryByText(/@/)).not.toBeInTheDocument();

    await user.click(preference);

    await waitFor(() => {
      expect(setTerminalEmailPreference).toHaveBeenCalledWith({
        terminal_email_enabled: false,
      });
      expect(preference).not.toBeChecked();
    });
  });

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
    expect(authoringLink).toHaveAttribute('title', 'New analysis');
    expect(authoringLink.querySelector('.sm\\:hidden')).toHaveTextContent('New');
    expect(authoringLink.querySelector('.sm\\:inline')).toHaveTextContent(
      'New analysis',
    );
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

  it('marks Artifacts current only on the exact list and two-segment detail', () => {
    const first = renderWithRouter(routes, {
      initialEntries: ['/artifacts'],
    });
    expect(screen.getByRole('link', { name: 'Artifacts' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    first.unmount();

    const detail = renderWithRouter(routes, {
      initialEntries: ['/artifacts/run-a/artifact-a'],
    });
    expect(screen.getByRole('link', { name: 'Artifacts' })).toHaveAttribute(
      'aria-current',
      'page',
    );
    detail.unmount();

    renderWithRouter(routes, {
      initialEntries: ['/artifacts/run-a/artifact-a/unknown'],
    });
    expect(screen.getByRole('link', { name: 'Artifacts' })).not.toHaveAttribute(
      'aria-current',
    );
  });
});
