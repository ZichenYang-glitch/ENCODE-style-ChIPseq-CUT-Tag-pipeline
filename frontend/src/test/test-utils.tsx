import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { createMemoryRouter, RouterProvider } from 'react-router-dom';
import { render } from '@testing-library/react';
import { ClientProvider, type ClientContextValue } from '../api/client-context';

export interface RenderWithRouterOptions {
  initialEntries?: string[];
  initialIndex?: number;
  clients?: Partial<ClientContextValue>;
}

export function renderWithRouter(
  routes: Parameters<typeof createMemoryRouter>[0],
  options: RenderWithRouterOptions = {},
) {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  const router = createMemoryRouter(routes, {
    initialEntries: options.initialEntries,
    initialIndex: options.initialIndex,
  });

  return {
    router,
    ...render(
      <QueryClientProvider client={queryClient}>
        <ClientProvider clients={options.clients}>
          <RouterProvider router={router} />
        </ClientProvider>
      </QueryClientProvider>,
    ),
  };
}
