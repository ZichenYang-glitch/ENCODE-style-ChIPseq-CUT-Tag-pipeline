import { describe, expect, it } from 'vitest';
import { ApiError } from '../../api/fetcher';
import type { AuthoringTurnResponse, AuthoringValidationResponse } from '../../api/generated/models';
import {
  canonicalDraftSnapshot,
  capabilityView,
  classifyAuthoringError,
  connectionView,
  draftSaveDisposition,
  draftWriteRequest,
  explanationView,
  isValidationBusinessFailure,
  proposalView,
  questionRequest,
  reasonCopy,
  receiptListParams,
  turnResponseIsCurrent,
  validationFromResponse,
  validationView,
  withMissingWorkflowSampleContext,
  workflowSampleRowsAbsent,
} from './authoringContract';
import { createCapabilityFixture } from './testFixtures';

const scope = {
  userId: 'user-1',
  projectId: 'prj_' + 'a'.repeat(32),
  workflowId: 'demo',
  sid: 'aas_' + 'b'.repeat(32),
  schemaDigest: 'c'.repeat(64),
  projectionDigest: 'd'.repeat(64),
};

describe('authoring contract mapping', () => {
  it('omits an empty receipt cursor instead of sending the string null', () => {
    expect(receiptListParams(null)).toEqual({ limit: 20 });
    expect(receiptListParams('cursor-2', 20)).toEqual({ cursor: 'cursor-2', limit: 20 });
  });

  it('maps known reason codes and hides unrecognized codes', () => {
    expect(reasonCopy('BRIDGE_UNAVAILABLE')).toMatch(/unavailable/i);
    expect(reasonCopy('HOST_PATH_/tmp/secret')).toBe('The request was not accepted.');
    expect(reasonCopy('HOST_PATH_/tmp/secret')).not.toMatch(/tmp|secret/);
  });

  it('does not treat a blocked bridge or a click as a ready connection', () => {
    expect(
      connectionView(
        { state: 'blocked', reason_code: 'BRIDGE_UNAVAILABLE' },
        { draftRestored: true, connectPending: false },
      ).state,
    ).toBe('bridge-unavailable');
    expect(
      connectionView({ state: 'closed' }, { draftRestored: true, connectPending: false }),
    ).toMatchObject({ state: 'disconnected', draftRestored: true });
    expect(
      connectionView({ state: 'closed' }, { draftRestored: false, connectPending: true }).state,
    ).toBe('connecting');
    expect(
      connectionView({ state: 'ready' }, { draftRestored: false, connectPending: false }).state,
    ).toBe('available');
  });

  it('sends only integer preferences and omits the explanation issue code', () => {
    const capability = createCapabilityFixture();
    const suggest = questionRequest(
      {
        intent: 'suggest',
        fieldTokens: ['resource.limit', 'resource.enabled'],
        preferences: [
          { fieldToken: 'resource.limit', value: 0 },
          { fieldToken: 'resource.enabled', value: false },
        ],
        issueCode: null,
      },
      capability,
    );
    expect(suggest.ok).toBe(true);
    if (!suggest.ok) return;
    expect(suggest.value.question).toEqual({
      intent: 'suggest',
      field_tokens: ['resource.limit', 'resource.enabled'],
      preferences: [{ field_token: 'resource.limit', value: 0 }],
    });
    expect(suggest.value.omittedNonIntegerPreferences).toBe(true);
    expect(JSON.stringify(suggest.value.question)).not.toMatch(/issue|false/);

    const explanation = questionRequest(
      {
        intent: 'explain_issues',
        fieldTokens: ['resource.limit'],
        preferences: [],
        issueCode: 'EXAMPLE_ISSUE',
      },
      capability,
    );
    expect(explanation.ok).toBe(true);
    if (!explanation.ok) return;
    expect(explanation.value.issueCodeNotSent).toBe(true);
    expect(explanation.value.question).toEqual({
      intent: 'explain_issues',
      field_tokens: [],
      preferences: [],
    });
  });

  it('keeps the server digest and refuses to adopt an explanation', () => {
    const capability = capabilityView({
      available: true,
      projection_version: 'hw-a1/1',
      allowed_intents: ['suggest', 'explain_issues'],
      allowed_issue_codes: [],
      fields: [
        {
          token: 'resource.limit',
          label: 'Example limit',
          type: 'integer',
          constraints: { type: 'integer' },
        },
      ],
    });
    const digest = 'ab'.repeat(32);
    const proposal = proposalView(turn({
      result_kind: 'proposal',
      adoptable: true,
      proposal_digest: digest,
      proposal_state: 'pending',
      proposal: {
        changes: [
          {
            field_token: 'resource.limit',
            value: 2,
            before: { present: false, token: 'resource.limit', value: null },
            after: { present: true, token: 'resource.limit', value: 2 },
          },
          {
            field_token: 'not-in-catalog',
            value: 1,
            before: { present: false, value: null },
            after: { present: true, value: 1 },
          },
        ],
        reason_codes: ['MODEL_OUTPUT_REJECTED'],
      },
    }), capability);
    expect(proposal?.digestAvailable).toBe(true);
    expect(proposal?.changes[0]).toMatchObject({
      label: 'Example limit',
      before: { presence: 'missing' },
      after: { presence: 'present', value: 2 },
      recognized: true,
    });
    expect(proposal?.changes[1]?.recognized).toBe(false);
    expect(JSON.stringify(proposal)).not.toContain(digest);

    const explanation = explanationView(
      turn({
        result_kind: 'explanation',
        adoptable: false,
        explanation: { reason_codes: ['SAMPLE_SELECTION_REQUIRED', 'not-a-public-code'] },
      }),
      {
        ...capability!,
        issues: [
          {
            code: 'SAMPLE_SELECTION_REQUIRED',
            label: 'Select at least one sample revision before validation.',
          },
        ],
      },
    );
    expect(proposalView(turn({
      result_kind: 'explanation',
      adoptable: false,
      explanation: { reason_codes: ['SAMPLE_SELECTION_REQUIRED'] },
    }), capability)).toBeNull();
    expect(explanation?.labels[0]).toMatch(/sample revision/i);
    expect(explanation?.labels[1]).toBe('This explanation is not in the current catalog.');
    expect(explanation?.labels.join(' ')).not.toMatch(/not-a-public-code/);
  });

  it('drops a turn response from another draft or catalog', () => {
    expect(turnResponseIsCurrent(scope, scope)).toBe(true);
    expect(turnResponseIsCurrent(scope, { ...scope, sid: 'aas_' + 'e'.repeat(32) })).toBe(false);
    expect(turnResponseIsCurrent(scope, { ...scope, projectionDigest: 'f'.repeat(64) })).toBe(false);
    expect(turnResponseIsCurrent(scope, { ...scope, userId: 'user-2' })).toBe(false);
  });

  it('labels an approved validation explanation without echoing an unknown code', () => {
    const view = explanationView(
      turn({
        result_kind: 'explanation',
        adoptable: false,
        explanation: {
          reason_codes: ['ENCODE_OPTIONS_INVALID', 'REFERENCE_PROFILE_REQUIRED', 'Traceback: secret'],
        },
      }),
      {
        available: true,
        unavailableReason: null,
        intents: ['explain_issues'],
        fields: [],
        issues: [
          { code: 'ENCODE_OPTIONS_INVALID', label: 'Options do not match the workflow.' },
          { code: 'REFERENCE_PROFILE_REQUIRED', label: 'A reference profile revision is required.' },
        ],
      },
    );
    expect(view?.labels).toEqual([
      'Options do not match the workflow.',
      'A reference profile revision is required.',
      'This explanation is not in the current catalog.',
    ]);
    expect(view?.labels.join(' ')).not.toMatch(/Traceback|secret/);
  });

  it('binds validation to the revision and schema that were checked', () => {
    const stored = validationFromResponse(
      {
        draft_revision: 3,
        ok: true,
        snapshot: null,
      },
      scope.schemaDigest,
    );
    expect(validationView({
      stored,
      currentRevision: 3,
      schemaDigest: scope.schemaDigest,
      dirty: false,
    }).status).toBe('no-snapshot');
    expect(validationView({
      stored,
      currentRevision: 4,
      schemaDigest: scope.schemaDigest,
      dirty: false,
    }).status).toBe('stale');
    expect(validationView({
      stored,
      currentRevision: 3,
      schemaDigest: scope.schemaDigest,
      dirty: true,
    }).status).toBe('stale');
  });

  it('keeps the real snapshot id, revision, and expiry from a successful validation', () => {
    const snapshotId = `vsnap_${'a'.repeat(32)}`;
    const stored = validationFromResponse(
      {
        draft_revision: 4,
        ok: true,
        snapshot: {
          snapshot_id: snapshotId,
          expires_at: '2026-11-17T00:00:00Z',
        } as NonNullable<AuthoringValidationResponse['snapshot']>,
      },
      scope.schemaDigest,
    );
    expect(stored.view).toMatchObject({
      status: 'snapshot',
      snapshotId,
      revisionLabel: 'Revision 4',
    });
    expect(stored.view.expiresLabel).toMatch(/2026/);
    expect(stored.view.snapshotId).not.toMatch(/generated|placeholder/i);
  });

  it('separates validation business failures from transport failures', () => {
    const business = classifyAuthoringError(
      new ApiError(409, 'SAMPLE_SELECTION_REQUIRED', 'server text', []),
    );
    const stale = classifyAuthoringError(new ApiError(409, 'STALE_DRAFT', 'server text', []));
    const lost = classifyAuthoringError(new TypeError('socket hang up'));
    expect(isValidationBusinessFailure(business)).toBe(true);
    expect(isValidationBusinessFailure(stale)).toBe(false);
    expect(lost.kind).toBe('network');
    expect(lost.message).not.toMatch(/socket|hang/);
    expect(stale.message).not.toBe('server text');
  });

  it('preserves server-managed samples without accepting a replacement path', () => {
    const body = draftWriteRequest({
      expectedRevision: 2,
      config: { limit: 1 },
      options: {},
      sampleRevisionIds: ['smpr_' + '1'.repeat(32)],
      referenceProfileRevisionId: 'refpr_' + '2'.repeat(32),
      preservedSamples: '/srv/private/samples.tsv',
      preservedInputSelections: [
        {
          input_use_key: 'fastq',
          occurrence: 0,
          input_file_revision_ids: ['ifr_' + '3'.repeat(32)],
        },
      ],
    });
    expect(body).toMatchObject({
      expected_revision: 2,
      sample_revision_ids: ['smpr_' + '1'.repeat(32)],
      reference_profile_revision_id: 'refpr_' + '2'.repeat(32),
      samples: '/srv/private/samples.tsv',
    });
    expect(body?.input_selections).toEqual([
      {
        input_use_key: 'fastq',
        occurrence: 0,
        input_file_revision_ids: ['ifr_' + '3'.repeat(32)],
      },
    ]);
  });

  it('compares the full normalized draft, including selections', () => {
    const sampleA = 'smpr_' + '1'.repeat(32);
    const sampleB = 'smpr_' + '2'.repeat(32);
    const referenceA = 'refpr_' + 'a'.repeat(32);
    const referenceB = 'refpr_' + 'b'.repeat(32);
    const selection = {
      input_use_key: 'reads',
      occurrence: 0,
      input_file_revision_ids: ['ifr_' + '3'.repeat(32)],
    };
    const sent = {
      config: { b: 1, a: 2 },
      options: { cores: 1 },
      samples: null,
      sample_revision_ids: [sampleA],
      reference_profile_revision_id: referenceA,
      input_selections: [selection],
    };
    expect(canonicalDraftSnapshot(sent)).toBe(
      canonicalDraftSnapshot({
        config: { a: 2, b: 1 },
        options: { cores: 1 },
        sample_revision_ids: [sampleA],
        reference_profile_revision_id: referenceA,
        input_selections: [selection],
      }),
    );
    expect(canonicalDraftSnapshot(sent)).not.toBe(
      canonicalDraftSnapshot({ ...sent, sample_revision_ids: [sampleB] }),
    );
    expect(canonicalDraftSnapshot(sent)).not.toBe(
      canonicalDraftSnapshot({ ...sent, sample_revision_ids: [sampleA, sampleB] }),
    );
    expect(canonicalDraftSnapshot(sent)).not.toBe(
      canonicalDraftSnapshot({ ...sent, reference_profile_revision_id: referenceB }),
    );
    expect(canonicalDraftSnapshot(sent)).not.toBe(
      canonicalDraftSnapshot({ ...sent, input_selections: [] }),
    );
    expect(canonicalDraftSnapshot(sent)).not.toBe(
      canonicalDraftSnapshot({ ...sent, samples: [{ sample: 'A' }] }),
    );
    expect(workflowSampleRowsAbsent(null)).toBe(true);
    expect(workflowSampleRowsAbsent([])).toBe(true);
    expect(workflowSampleRowsAbsent('/srv/private/samples.tsv')).toBe(false);
    expect(
      withMissingWorkflowSampleContext(
        'ENCODE reference binding is unavailable.',
        [{ code: 'ENCODE_REFERENCE_BINDING_CONFLICT' }],
        null,
      ),
    ).toMatch(/no workflow sample rows/i);
    expect(
      withMissingWorkflowSampleContext(
        'ENCODE reference binding is unavailable.',
        [{ code: 'ENCODE_REFERENCE_BINDING_CONFLICT' }],
        '/srv/private/samples.tsv',
      ),
    ).toBe('ENCODE reference binding is unavailable.');
  });

  it('confirms a matching draft without treating a readback as exactly-once', () => {
    const sent = {
      config: { limit: 5 },
      options: {},
      samples: null,
      sample_revision_ids: [] as string[],
      reference_profile_revision_id: null,
      input_selections: [],
    };
    expect(
      draftSaveDisposition({
        sent,
        observed: sent,
        previousRevision: 1,
        observedRevision: 2,
        acknowledged: false,
      }),
    ).toBe('confirmed');
    expect(
      draftSaveDisposition({
        sent,
        observed: sent,
        previousRevision: 1,
        observedRevision: 1,
        acknowledged: false,
      }),
    ).toBe('unconfirmed');
    expect(
      draftSaveDisposition({
        sent,
        observed: { ...sent, sample_revision_ids: ['smpr_' + '2'.repeat(32)] },
        previousRevision: 1,
        observedRevision: 2,
        acknowledged: false,
      }),
    ).toBe('conflict');
    expect(
      draftSaveDisposition({
        sent,
        observed: { ...sent, reference_profile_revision_id: 'refpr_' + 'b'.repeat(32) },
        previousRevision: 1,
        observedRevision: 2,
        acknowledged: true,
      }),
    ).toBe('conflict');
  });
});

function turn(overrides: Partial<AuthoringTurnResponse>): AuthoringTurnResponse {
  return {
    ok: true,
    turn_id: 'atr_' + '1'.repeat(32),
    client_turn_id: 'client1',
    status: 'completed',
    adoptable: false,
    ...overrides,
  };
}
