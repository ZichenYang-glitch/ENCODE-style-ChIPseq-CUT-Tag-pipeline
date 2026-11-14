import {
  AUTHORING_WRITES_UNAVAILABLE,
  type AssistantConnectionState,
  type ProposalView,
} from './model';

export const UNSAVED_EDIT_BLOCK =
  'Save or discard unsaved edits before this action.';

export function writeBlockReason(
  writesEnabled: boolean,
  dirty: boolean,
  unavailable = AUTHORING_WRITES_UNAVAILABLE,
): string | null {
  if (!writesEnabled) return unavailable;
  if (dirty) return UNSAVED_EDIT_BLOCK;
  return null;
}

export function connectBlockReason(
  writesEnabled: boolean,
  state: AssistantConnectionState,
  unavailable = AUTHORING_WRITES_UNAVAILABLE,
): string | null {
  if (state === 'connecting') return 'The assistant is still connecting.';
  if (state === 'available') return 'The assistant is already available.';
  if (!writesEnabled) return unavailable;
  return null;
}

export function adoptBlockReason(input: {
  writesEnabled: boolean;
  dirty: boolean;
  proposal: ProposalView | null;
  unavailable?: string;
}): string | null {
  const blocked = writeBlockReason(input.writesEnabled, input.dirty, input.unavailable);
  if (blocked) return blocked;
  return proposalDecisionBlock(input.proposal, 'adopt');
}

export function rejectBlockReason(input: {
  writesEnabled: boolean;
  dirty: boolean;
  proposal: ProposalView | null;
  unavailable?: string;
}): string | null {
  const blocked = writeBlockReason(input.writesEnabled, input.dirty, input.unavailable);
  if (blocked) return blocked;
  return proposalDecisionBlock(input.proposal, 'reject');
}

function proposalDecisionBlock(
  proposal: ProposalView | null,
  action: 'adopt' | 'reject',
): string | null {
  if (!proposal) return 'There is no suggestion to review.';
  if (!proposal.digestAvailable) {
    return 'This suggestion cannot be reviewed until the server provides its review identity.';
  }
  switch (proposal.state) {
    case 'pending':
      if (action === 'adopt' && !proposal.adoptable) {
        return 'This suggestion cannot be adopted.';
      }
      return null;
    case 'adopting':
      return 'Adoption is already in progress.';
    case 'adopted':
      return 'This suggestion is already adopted.';
    case 'rejected':
      return 'This suggestion was rejected. The draft was not changed.';
    case 'expired':
      return 'This suggestion has expired and was not adopted.';
    case 'conflict':
      return 'The draft changed before this action completed. The suggestion was not adopted.';
    case 'invalidated':
      return 'This suggestion is no longer valid and was not adopted.';
    case 'unconfirmed':
      return 'The last response was not confirmed. Refresh the server status before trying again.';
    default: {
      const unreachable: never = proposal.state;
      return unreachable;
    }
  }
}

export function validateBlockReason(input: {
  writesEnabled: boolean;
  dirty: boolean;
  draftExpired: boolean;
  unavailable?: string;
}): string | null {
  if (input.draftExpired) {
    return 'This draft has expired, so it cannot be validated.';
  }
  return writeBlockReason(input.writesEnabled, input.dirty, input.unavailable);
}
