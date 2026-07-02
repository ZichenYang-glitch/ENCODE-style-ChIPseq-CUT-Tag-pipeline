import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { createAgentApiClient } from './agentClient';
import type { AgentRequest } from './types';

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

function mockFetch(response: Partial<Response>) {
  return vi.fn().mockResolvedValue(response as Response);
}

const sampleRequest: AgentRequest = {
  session_id: null,
  message: 'Explain the issues.',
  context: {
    current_issues: [],
    current_config: {},
    current_schema: {},
  },
};

describe('createAgentApiClient', () => {
  let originalFetch: typeof globalThis.fetch;

  beforeEach(() => {
    originalFetch = globalThis.fetch;
  });

  afterEach(() => {
    globalThis.fetch = originalFetch;
    vi.restoreAllMocks();
  });

  it('POSTs to /api/v1/workflows/{id}/agent/chat with JSON body', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => ({
        ok: true,
        session_id: null,
        message: 'Deterministic mock explanation.',
        suggestions: [],
        tool_calls: [],
        issues: [],
      }),
    });
    globalThis.fetch = fetchMock;

    const client = createAgentApiClient('');
    const response = await client.chat(WORKFLOW_ID, sampleRequest);

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/workflows/${WORKFLOW_ID}/agent/chat`,
      expect.objectContaining({
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(sampleRequest),
        credentials: 'omit',
      }),
    );
    expect(response.message).toBe('Deterministic mock explanation.');
  });

  it('uses the provided baseUrl', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => ({
        ok: true,
        session_id: null,
        message: 'ok',
        suggestions: [],
        tool_calls: [],
        issues: [],
      }),
    });
    globalThis.fetch = fetchMock;

    const client = createAgentApiClient('http://localhost:8000');
    await client.chat(WORKFLOW_ID, sampleRequest);

    expect(fetchMock).toHaveBeenCalledWith(
      `http://localhost:8000/api/v1/workflows/${WORKFLOW_ID}/agent/chat`,
      expect.anything(),
    );
  });

  it('throws AgentApiError on non-2xx response', async () => {
    const fetchMock = mockFetch({
      ok: false,
      status: 400,
      statusText: 'Bad Request',
    });
    globalThis.fetch = fetchMock;

    const client = createAgentApiClient('');
    await expect(client.chat(WORKFLOW_ID, sampleRequest)).rejects.toMatchObject({
      status: 400,
      statusText: 'Bad Request',
    });
  });
});
