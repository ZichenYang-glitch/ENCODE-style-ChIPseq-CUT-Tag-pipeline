import { describe, it, expect } from 'vitest';
import { renderHook } from '@testing-library/react';
import { ClientProvider, useClients } from './client-context';

function wrapper({ children }: { children: React.ReactNode }) {
  return <ClientProvider>{children}</ClientProvider>;
}

describe('ClientProvider', () => {
  it('provides default generated clients', () => {
    const { result } = renderHook(() => useClients(), { wrapper });
    expect(result.current.workflowClient).toBeDefined();
    expect(result.current.runClient).toBeDefined();
    expect(result.current.agentClient).toBeDefined();
  });
});
