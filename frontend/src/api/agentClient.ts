import type { AgentApiError, AgentRequest, AgentResponse } from './types';

export interface AgentApiClient {
  chat(workflowId: string, request: AgentRequest): Promise<AgentResponse>;
}

function getBaseUrl(): string {
  const envBase = import.meta.env?.VITE_API_BASE_URL;
  return typeof envBase === 'string' ? envBase : '';
}

export function createAgentApiClient(baseUrl?: string): AgentApiClient {
  const root = baseUrl ?? getBaseUrl();

  return {
    async chat(workflowId: string, request: AgentRequest): Promise<AgentResponse> {
      const url = `${root}/api/v1/workflows/${encodeURIComponent(workflowId)}/agent/chat`;

      const response = await fetch(url, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(request),
        credentials: 'omit',
      });

      if (!response.ok) {
        const error = new Error(
          `Agent chat request failed: ${response.status} ${response.statusText}`,
        ) as AgentApiError;
        error.status = response.status;
        error.statusText = response.statusText;
        throw error;
      }

      return response.json() as Promise<AgentResponse>;
    },
  };
}
