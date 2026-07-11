import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { createMemoryRouter, RouterProvider } from 'react-router-dom';
import { render } from '@testing-library/react';
import { ClientProvider, type ClientContextValue } from '../api/client-context';
import { createStubWorkflowClient } from '../api/client';
import { createStubRunApiClient } from '../api/runClient';
import type { AgentApiClient } from '../api/agentClient';

export interface RenderWithRouterOptions {
  initialEntries?: string[];
  initialIndex?: number;
  clients?: Partial<ClientContextValue>;
}

export function renderWithRouter(
  routes: Parameters<typeof createMemoryRouter>[0],
  options: RenderWithRouterOptions = {},
) {
  const stubAgentClient: AgentApiClient = {
    async chat() {
      return {
        ok: true,
        session_id: null,
        message: 'Stub agent response.',
        suggestions: [],
        tool_calls: [],
        issues: [],
      };
    },
  };
  const clients: ClientContextValue = {
    workflowClient: createStubWorkflowClient(),
    runClient: createStubRunApiClient(),
    agentClient: stubAgentClient,
    ...options.clients,
  };
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
        <ClientProvider clients={clients}>
          <RouterProvider router={router} />
        </ClientProvider>
      </QueryClientProvider>,
    ),
  };
}
