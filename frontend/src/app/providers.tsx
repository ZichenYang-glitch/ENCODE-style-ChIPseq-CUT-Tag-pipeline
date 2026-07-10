import { useMemo } from 'react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { RouterProvider } from 'react-router-dom';
import {
  ClientProvider,
  type ClientProviderProps,
} from '../api/client-context';
import { createAppRouter } from './router';

export interface ProvidersProps {
  queryClient?: QueryClient;
  clients?: ClientProviderProps['clients'];
}

const defaultQueryClient = new QueryClient();

export function Providers({
  queryClient = defaultQueryClient,
  clients,
}: ProvidersProps) {
  const router = useMemo(() => createAppRouter(), []);

  return (
    <QueryClientProvider client={queryClient}>
      <ClientProvider clients={clients}>
        <RouterProvider router={router} />
      </ClientProvider>
    </QueryClientProvider>
  );
}
