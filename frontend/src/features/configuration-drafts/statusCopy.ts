import type {
  AssistantConnectionState,
  DraftValidationView,
  ProposalDecisionState,
} from './model';

export function connectionHeading(state: AssistantConnectionState): string {
  switch (state) {
    case 'disconnected':
      return 'Assistant not connected';
    case 'connecting':
      return 'Connecting the assistant';
    case 'available':
      return 'Assistant available';
    case 'ended':
      return 'Assistant connection ended';
    case 'bridge-unavailable':
      return 'Assistant unavailable';
    default: {
      const unreachable: never = state;
      return unreachable;
    }
  }
}

export function connectionDetail(input: {
  state: AssistantConnectionState;
  draftRestored: boolean;
}): string {
  if (input.draftRestored && input.state === 'bridge-unavailable') {
    return 'Your draft is still here after signing in. The assistant is unavailable, so it was not reconnected. You can keep editing by hand.';
  }
  if (
    input.draftRestored &&
    input.state !== 'available' &&
    input.state !== 'connecting'
  ) {
    return 'Your draft is still here after signing in. The assistant is not connected. Connect it again before requesting a suggestion.';
  }
  switch (input.state) {
    case 'disconnected':
      return 'Connect the assistant before requesting a suggestion. Opening this page does not connect it.';
    case 'connecting':
      return 'The assistant connection is still in progress.';
    case 'available':
      return 'The assistant can answer a structured suggestion request.';
    case 'ended':
      return 'The assistant connection has ended. Connect it again before requesting a suggestion.';
    case 'bridge-unavailable':
      return 'The assistant is unavailable. You can keep editing by hand.';
    default: {
      const unreachable: never = input.state;
      return unreachable;
    }
  }
}

export function proposalStatusCopy(state: ProposalDecisionState): string {
  switch (state) {
    case 'pending':
      return 'Suggestion waiting for your review. It is not applied to the draft.';
    case 'adopting':
      return 'Adoption is in progress. This page is waiting for the server result.';
    case 'adopted':
      return 'Suggestion adopted.';
    case 'rejected':
      return 'Suggestion rejected. The draft was not changed.';
    case 'expired':
      return 'This suggestion has expired and was not adopted.';
    case 'conflict':
      return 'The draft changed before this action completed. The suggestion was not adopted.';
    case 'invalidated':
      return 'This suggestion is no longer valid and was not adopted.';
    case 'unconfirmed':
      return 'The last response was not confirmed. Refresh the server status before trying again. This page will not retry the action.';
    default: {
      const unreachable: never = state;
      return unreachable;
    }
  }
}

export function validationLead(view: DraftValidationView): string {
  switch (view.status) {
    case 'unavailable':
      return 'Draft validation is not connected yet.';
    case 'unconfirmed':
      return 'The validation response was not confirmed. This page has no validation result.';
    case 'issues':
      return 'Validation reported problems.';
    case 'no-snapshot':
      return 'The configuration check passed, and there is no execution snapshot.';
    case 'snapshot':
      return 'Validation finished. This does not start a run.';
    default: {
      const unreachable: never = view.status;
      return unreachable;
    }
  }
}
