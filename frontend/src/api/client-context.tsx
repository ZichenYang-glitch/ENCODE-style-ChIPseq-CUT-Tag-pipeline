import { createContext, useContext, useMemo, type ReactNode } from 'react';
import { createAgentApiClient, type AgentApiClient } from './agentClient';
import { createStubRunApiClient, type RunApiClient } from './runClient';
import { createStubWorkflowClient, type WorkflowApiClient } from './client';

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
      workflowClient: clients?.workflowClient ?? createStubWorkflowClient(),
      runClient: clients?.runClient ?? createStubRunApiClient(),
      agentClient: clients?.agentClient ?? createAgentApiClient(),
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
