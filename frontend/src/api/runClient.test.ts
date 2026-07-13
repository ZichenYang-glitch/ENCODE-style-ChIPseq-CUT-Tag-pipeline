import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { createRunApiClient } from './runClient';
import type { RunCreateRequest } from './runTypes';

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';
const RUN_ID = 'run-123';

function mockFetch(response: Partial<Response>) {
  return vi.fn().mockResolvedValue(response as Response);
}

const sampleCreateRequest: RunCreateRequest = {
  snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
  tags: { env: 'test' },
};

const sampleRunRecord = {
  run_id: RUN_ID,
  workflow_id: WORKFLOW_ID,
  inputs: {
    config: { genome: 'hg38' },
    samples: [{ name: 'sample-1', fastq_r1: 's1_R1.fq.gz' }],
    options: { threads: 4 },
  },
  status: 'created',
  created_at: '2026-07-03T12:00:00Z',
  updated_at: '2026-07-03T12:00:00Z',
  started_at: null,
  ended_at: null,
  current_stage: null,
  cancellation_reason: null,
  error: null,
  tags: { env: 'test' },
};

const sampleRunResponse = {
  ok: true,
  run: sampleRunRecord,
  issues: [],
};

const sampleEventsResponse = {
  ok: true,
  run_id: RUN_ID,
  events: [
    {
      event_id: 'evt-1',
      run_id: RUN_ID,
      sequence: 1,
      event_type: 'run.created',
      timestamp: '2026-07-03T12:00:00Z',
      status: 'created',
      stage: null,
      message: 'Run created.',
      context: {},
      issue: null,
    },
  ],
  next_cursor: null,
  issues: [],
};

const sampleLogsResponse = {
  ok: true,
  run_id: RUN_ID,
  stream_name: 'stdout',
  chunks: [
    {
      chunk_id: 'chunk-1',
      run_id: RUN_ID,
      stream_name: 'stdout',
      sequence: 1,
      timestamp: '2026-07-03T12:00:00Z',
      lines: ['line 1'],
    },
  ],
  next_cursor: null,
  issues: [],
};

describe('createRunApiClient', () => {
  let originalFetch: typeof globalThis.fetch;

  beforeEach(() => {
    originalFetch = globalThis.fetch;
  });

  afterEach(() => {
    globalThis.fetch = originalFetch;
    vi.restoreAllMocks();
  });

  it('POSTs to /api/v1/workflows/{id}/runs with JSON body', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 201,
      statusText: 'Created',
      json: async () => sampleRunResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    const response = await client.createRun(WORKFLOW_ID, sampleCreateRequest);

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/workflows/${WORKFLOW_ID}/runs`,
      expect.objectContaining({
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(sampleCreateRequest),
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
    expect(response.run).toEqual(sampleRunRecord);
  });

  it('GETs /api/v1/runs/{id}', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => sampleRunResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    const response = await client.getRun(RUN_ID);

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/runs/${RUN_ID}`,
      expect.objectContaining({
        method: 'GET',
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
    expect(response.run?.run_id).toBe(RUN_ID);
  });

  it('POSTs to /api/v1/runs/{id}/preflight with empty body', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 202,
      statusText: 'Accepted',
      json: async () => sampleRunResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    const response = await client.preflightRun(RUN_ID);

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/runs/${RUN_ID}/preflight`,
      expect.objectContaining({
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({}),
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
  });

  it('POSTs to /api/v1/runs/{id}/start with empty body', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 202,
      statusText: 'Accepted',
      json: async () => sampleRunResponse,
    });
    globalThis.fetch = fetchMock;

    const response = await createRunApiClient('').startRun(RUN_ID);

    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/runs/${RUN_ID}/start`,
      expect.objectContaining({
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({}),
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
  });

  it('GETs /api/v1/runs/{id}/events with after and limit', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => sampleEventsResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    const response = await client.listRunEvents(RUN_ID, { after: 'evt-0', limit: 10 });

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/runs/${RUN_ID}/events?after=evt-0&limit=10`,
      expect.objectContaining({
        method: 'GET',
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
    expect(response.events).toHaveLength(1);
  });

  it('GETs /api/v1/runs/{id}/logs with stream_name, after, and limit', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => sampleLogsResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    const response = await client.listRunLogs(RUN_ID, {
      streamName: 'stderr',
      after: 'chunk-0',
      limit: 25,
    });

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/runs/${RUN_ID}/logs?stream_name=stderr&after=chunk-0&limit=25`,
      expect.objectContaining({
        method: 'GET',
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
    expect(response.chunks).toHaveLength(1);
  });

  it('POSTs to /api/v1/runs/{id}/cancel with empty body', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => sampleRunResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    const response = await client.cancelRun(RUN_ID);

    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/runs/${RUN_ID}/cancel`,
      expect.objectContaining({
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({}),
        credentials: 'omit',
      }),
    );
    expect(response.ok).toBe(true);
  });

  it('uses the provided baseUrl', async () => {
    const fetchMock = mockFetch({
      ok: true,
      status: 200,
      statusText: 'OK',
      json: async () => sampleRunResponse,
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('http://localhost:8000');
    await client.getRun(RUN_ID);

    expect(fetchMock).toHaveBeenCalledWith(
      `http://localhost:8000/api/v1/runs/${RUN_ID}`,
      expect.anything(),
    );
  });

  it('throws RunApiError on non-2xx response', async () => {
    const fetchMock = mockFetch({
      ok: false,
      status: 404,
      statusText: 'Not Found',
    });
    globalThis.fetch = fetchMock;

    const client = createRunApiClient('');
    await expect(client.getRun(RUN_ID)).rejects.toMatchObject({
      status: 404,
      statusText: 'Not Found',
    });
  });
});

import { createStubRunApiClient } from './runClient';

describe('createStubRunApiClient', () => {
  it('returns a synthetic created run', async () => {
    const client = createStubRunApiClient();
    const response = await client.createRun(WORKFLOW_ID, sampleCreateRequest);

    expect(response.ok).toBe(true);
    expect(response.run).not.toBeNull();
    expect(response.run?.status).toBe('created');
    expect(response.run?.workflow_id).toBe(WORKFLOW_ID);
    expect(response.run?.inputs).toEqual({
      validated_snapshot_id: sampleCreateRequest.snapshot_id,
    });
    expect(response.run?.tags).toEqual({ env: 'test' });
    expect(response.run?.run_id).toMatch(/^stub-run-\d+$/);
  });

  it('returns the created run from getRun and reports unknown runs', async () => {
    const client = createStubRunApiClient();
    const createResponse = await client.createRun(WORKFLOW_ID, sampleCreateRequest);
    const runId = createResponse.run!.run_id;

    const getResponse = await client.getRun(runId);
    expect(getResponse.ok).toBe(true);
    expect(getResponse.run?.run_id).toBe(runId);

    const unknownResponse = await client.getRun('no-such-run');
    expect(unknownResponse.ok).toBe(false);
    expect(unknownResponse.issues[0].code).toBe('RUN_NOT_FOUND');
  });

  it('transitions a created run to cancelled and preserves terminal runs', async () => {
    const client = createStubRunApiClient();
    const createResponse = await client.createRun(WORKFLOW_ID, sampleCreateRequest);
    const runId = createResponse.run!.run_id;

    const cancelResponse = await client.cancelRun(runId);
    expect(cancelResponse.ok).toBe(true);
    expect(cancelResponse.run?.status).toBe('cancelled');
    expect(cancelResponse.run?.cancellation_reason).toBe('User requested cancellation.');

    const secondCancelResponse = await client.cancelRun(runId);
    expect(secondCancelResponse.ok).toBe(true);
    expect(secondCancelResponse.run?.status).toBe('cancelled');
  });

  it('preflights a created run exactly once', async () => {
    const client = createStubRunApiClient();
    const createResponse = await client.createRun(WORKFLOW_ID, sampleCreateRequest);
    const runId = createResponse.run!.run_id;

    const preflightResponse = await client.preflightRun(runId);
    expect(preflightResponse.ok).toBe(true);
    expect(preflightResponse.run).toMatchObject({
      status: 'planned',
      current_stage: 'preflight',
    });

    const eventsResponse = await client.listRunEvents(runId);
    expect(eventsResponse.events.at(-1)).toMatchObject({
      event_type: 'preflight_completed',
      status: 'planned',
    });

    const repeatedResponse = await client.preflightRun(runId);
    expect(repeatedResponse.ok).toBe(false);
    expect(repeatedResponse.issues[0].code).toBe('PREFLIGHT_ALREADY_TRIGGERED');
  });

  it('returns the created event, and a cancelled event after cancellation', async () => {
    const client = createStubRunApiClient();
    const createResponse = await client.createRun(WORKFLOW_ID, sampleCreateRequest);
    const runId = createResponse.run!.run_id;

    const eventsBeforeCancel = await client.listRunEvents(runId);
    expect(eventsBeforeCancel.ok).toBe(true);
    expect(eventsBeforeCancel.events).toHaveLength(1);
    expect(eventsBeforeCancel.events[0]).toMatchObject({
      event_type: 'status_changed',
      status: 'created',
      context: { previous_status: null, new_status: 'created' },
    });

    await client.cancelRun(runId);

    const eventsAfterCancel = await client.listRunEvents(runId);
    expect(eventsAfterCancel.events).toHaveLength(2);
    expect(eventsAfterCancel.events[1]).toMatchObject({
      event_type: 'status_changed',
      status: 'cancelled',
      context: {
        previous_status: 'created',
        new_status: 'cancelled',
        cancellation_reason: 'User requested cancellation.',
      },
    });
  });

  it('returns an empty chunk list by default', async () => {
    const client = createStubRunApiClient();
    const createResponse = await client.createRun(WORKFLOW_ID, sampleCreateRequest);
    const runId = createResponse.run!.run_id;

    const logsResponse = await client.listRunLogs(runId);
    expect(logsResponse.ok).toBe(true);
    expect(logsResponse.chunks).toEqual([]);
    expect(logsResponse.stream_name).toBe('stdout');
  });
});
