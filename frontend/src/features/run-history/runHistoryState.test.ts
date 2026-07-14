import { describe, expect, it } from 'vitest';
import type { RunHistoryResponse } from '../../api/generated/models';
import {
  flattenRunHistoryPages,
  parseRunHistoryFilters,
  RunHistoryProtocolError,
  runHistoryPaginationAnomaly,
  validateRunHistoryResponse,
} from './runHistoryState';

function response(
  runId = 'run-a',
  nextCursor: string | null = null,
): RunHistoryResponse {
  return {
    ok: true,
    runs: [
      {
        run_id: runId,
        workflow_id: 'workflow-a',
        status: 'succeeded',
        created_at: '2026-07-14T08:00:00Z',
        updated_at: '2026-07-14T08:01:00Z',
        started_at: '2026-07-14T08:00:01Z',
        ended_at: '2026-07-14T08:01:00Z',
        current_stage: 'execution',
      },
    ],
    next_cursor: nextCursor,
    issues: [],
  };
}

describe('run history URL and protocol state', () => {
  it('parses exact workflow and status filters', () => {
    expect(
      parseRunHistoryFilters(
        new URLSearchParams('workflow_id=workflow-a&status=succeeded'),
      ),
    ).toEqual({
      ok: true,
      filters: { workflowId: 'workflow-a', status: 'succeeded' },
    });
  });

  it.each([
    'status=unknown',
    'status=succeeded&status=failed',
    'workflow_id=workflow-a&workflow_id=workflow-b',
    'workflow_id=%2Fprivate',
    'workflow_id=workflow-a%0Asecret',
    `workflow_id=${'x'.repeat(256)}`,
  ])('rejects unsafe or ambiguous filters without returning their value: %s', (query) => {
    expect(parseRunHistoryFilters(new URLSearchParams(query))).toEqual({
      ok: false,
      filters: null,
    });
  });

  it('projects a strict response and rejects unsafe cursors or timestamps', () => {
    expect(validateRunHistoryResponse(response(), undefined).runs[0]?.run_id).toBe(
      'run-a',
    );
    expect(() =>
      validateRunHistoryResponse(response('run-a', 'cursor with spaces'), undefined),
    ).toThrow(RunHistoryProtocolError);
    const invalid = response();
    invalid.runs[0]!.updated_at = 'not-a-time';
    expect(() => validateRunHistoryResponse(invalid, undefined)).toThrow(
      RunHistoryProtocolError,
    );
  });

  it('detects duplicate cursors and run identities while retaining first rows', () => {
    const first = response('run-a', 'runhist_next');
    const second = response('run-a', null);
    expect(runHistoryPaginationAnomaly([first, second], [undefined, 'runhist_next'])).toBe(
      true,
    );
    expect(flattenRunHistoryPages([first, second]).map((run) => run.run_id)).toEqual([
      'run-a',
    ]);
  });

  it('accepts the normal previous-next to following-page cursor handoff', () => {
    const first = response('run-b', 'runhist_next');
    const second = response('run-a', null);
    expect(
      runHistoryPaginationAnomaly([first, second], [undefined, 'runhist_next']),
    ).toBe(false);
  });
});
