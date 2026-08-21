import { describe, expect, it, vi, beforeEach } from 'vitest';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { createMemoryRouter, RouterProvider } from 'react-router-dom';
import { AuthProvider } from '../../app/auth';
import { appRoutes } from '../../app/router';
import { LoginPage } from '../login';
import * as authClient from '../../api/authClient';
import { ApiError } from '../../api/fetcher';
import { renderWithRouter } from '../../test/test-utils';

function renderLogin() {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  const router = createMemoryRouter(
    [{ path: '/login', element: <LoginPage /> }],
    { initialEntries: ['/login'] },
  );
  render(
    <QueryClientProvider client={queryClient}>
      <AuthProvider>
        <RouterProvider router={router} />
      </AuthProvider>
    </QueryClientProvider>,
  );
  return router;
}

const mockedFetchSessionState = vi.mocked(authClient.fetchSessionState);
const mockedLogin = vi.mocked(authClient.login);

beforeEach(() => {
  mockedFetchSessionState.mockResolvedValue({
    ok: true,
    setup_required: false,
    authenticated: true,
    principal: {
      user_id: 'usr_00000000000000000000000000000000',
      username: 'test-admin',
      role: 'administrator',
    },
  });
});

describe('session guard', () => {
  it('redirects anonymous visitors to /login with the origin path', async () => {
    mockedFetchSessionState.mockResolvedValue({
      ok: true,
      setup_required: false,
      authenticated: false,
      principal: null,
    });
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs'],
    });
    await screen.findByRole('heading', { name: 'HelixWeave' });
    expect(router.state.location.pathname).toBe('/login');
    expect(router.state.location.state).toEqual({ from: '/runs' });
  });

  it('lets an authenticated principal through to the requested page', async () => {
    renderWithRouter(appRoutes, { initialEntries: ['/workflows'] });
    expect(
      await screen.findByRole('heading', { name: 'HelixWeave' }),
    ).toBeInTheDocument();
    expect(
      await screen.findByTestId('current-username'),
    ).toHaveTextContent('test-admin');
  });
});

describe('login page', () => {
  it('shows the bootstrap instruction while setup is required', async () => {
    mockedFetchSessionState.mockResolvedValue({
      ok: true,
      setup_required: true,
      authenticated: false,
      principal: null,
    });
    renderLogin();
    expect(
      await screen.findByText(/helixweave admin account bootstrap/),
    ).toBeInTheDocument();
    expect(screen.queryByLabelText('Username')).not.toBeInTheDocument();
  });

  it('signs in and navigates back to the origin path', async () => {
    mockedFetchSessionState
      .mockResolvedValueOnce({
        ok: true,
        setup_required: false,
        authenticated: false,
        principal: null,
      })
      .mockResolvedValue({
        ok: true,
        setup_required: false,
        authenticated: true,
        principal: {
          user_id: 'usr_00000000000000000000000000000000',
          username: 'root-admin',
          role: 'administrator',
        },
      });
    const router = renderLogin();
    await userEvent.type(screen.getByLabelText('Username'), 'root-admin');
    await userEvent.type(
      screen.getByLabelText('Password'),
      'correct horse battery staple',
    );
    await userEvent.click(screen.getByRole('button', { name: 'Sign in' }));
    await waitFor(() => {
      expect(mockedLogin).toHaveBeenCalledWith(
        'root-admin',
        'correct horse battery staple',
      );
    });
    await waitFor(() => {
      expect(router.state.location.pathname).toBe('/workflows');
    });
  });

  it('shows the uniform rejection without echoing credentials', async () => {
    mockedFetchSessionState.mockResolvedValue({
      ok: true,
      setup_required: false,
      authenticated: false,
      principal: null,
    });
    mockedLogin.mockRejectedValueOnce(
      new ApiError(401, 'INVALID_CREDENTIALS', 'The username or password is invalid.'),
    );
    renderLogin();
    await userEvent.type(screen.getByLabelText('Username'), 'root-admin');
    await userEvent.type(screen.getByLabelText('Password'), 'wrong-password');
    await userEvent.click(screen.getByRole('button', { name: 'Sign in' }));
    expect(await screen.findByRole('alert')).toHaveTextContent(
      'The username or password is invalid.',
    );
  });
});
