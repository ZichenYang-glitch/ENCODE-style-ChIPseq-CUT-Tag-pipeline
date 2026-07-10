import { createContext, useContext, useMemo, type ReactNode } from 'react';
import type { AgentApiClient } from './agentClient';
import type { RunApiClient } from './runClient';
import type { WorkflowApiClient } from './client';
import {
  createGeneratedAgentClient,
  createGeneratedRunClient,
  createGeneratedWorkflowClient,
} from './generated-client-adapters';

export interface ClientContextValue {
  workflowClient: WorkflowApiClient;
  runClient: RunApiClient;
  agentClient: AgentApiClient;
}

const ClientContext = createContext<ClientContextValue | null>(null);

export interface ClientProviderProps {
  children: ReactNode;
  clients?: Partial<ClientContextValue>;
}

export function ClientProvider({ children, clients }: ClientProviderProps) {
  const value = useMemo<ClientContextValue>(
    () => ({
      workflowClient: clients?.workflowClient ?? createGeneratedWorkflowClient(),
      runClient: clients?.runClient ?? createGeneratedRunClient(),
      agentClient: clients?.agentClient ?? createGeneratedAgentClient(),
    }),
    [clients],
  );

  return (
    <ClientContext.Provider value={value}>{children}</ClientContext.Provider>
  );
}

export function useClients(): ClientContextValue {
  const context = useContext(ClientContext);
  if (!context) {
    throw new Error('useClients must be used within a ClientProvider');
  }
  return context;
}
