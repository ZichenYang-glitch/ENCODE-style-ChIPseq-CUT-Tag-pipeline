import { describe, expect, it } from 'vitest';
import { adoptBlockReason, connectBlockReason } from './actionGate';
import { createPendingProposal } from './testFixtures';

describe('authoring action gates', () => {
  it('keeps an unavailable assistant from looking connected', () => {
    expect(connectBlockReason(true, 'bridge-unavailable')).toBe(
      'The assistant is unavailable. You can keep editing by hand.',
    );
    expect(connectBlockReason(false, 'disconnected')).toMatch(/reviewed configuration-draft interface/);
  });

  it('does not treat an unconfirmed or expired suggestion as adoptable', () => {
    expect(
      adoptBlockReason({
        writesEnabled: true,
        dirty: false,
        proposal: createPendingProposal({ state: 'unconfirmed' }),
      }),
    ).toMatch(/not confirmed/);
    expect(
      adoptBlockReason({
        writesEnabled: true,
        dirty: false,
        proposal: createPendingProposal({ state: 'expired' }),
      }),
    ).toMatch(/expired/);
    expect(
      adoptBlockReason({
        writesEnabled: true,
        dirty: false,
        proposal: createPendingProposal({ state: 'conflict' }),
      }),
    ).toMatch(/not adopted/);
    expect(
      adoptBlockReason({
        writesEnabled: true,
        dirty: true,
        proposal: createPendingProposal(),
      }),
    ).toMatch(/unsaved edits/);
    expect(
      adoptBlockReason({
        writesEnabled: true,
        dirty: false,
        proposal: createPendingProposal({ digestAvailable: false }),
      }),
    ).toMatch(/review identity/);
  });
});
