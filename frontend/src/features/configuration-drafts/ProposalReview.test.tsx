import { render, screen, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { ExplanationNotice } from './ExplanationNotice';
import { ProposalReview } from './ProposalReview';
import { createPendingProposal } from './testFixtures';

describe('proposal review', () => {
  it('shows missing, null, zero, false, and empty values without editing them', () => {
    render(
      <ProposalReview
        proposal={createPendingProposal({
          changes: [
            {
              label: 'Missing value',
              before: { presence: 'missing' },
              after: { presence: 'present', value: 0 },
              recognized: true,
            },
            {
              label: 'Null value',
              before: { presence: 'null' },
              after: { presence: 'present', value: false },
              recognized: true,
            },
            {
              label: 'Empty value',
              before: { presence: 'present', value: '' },
              after: { presence: 'present', value: 'kept' },
              recognized: true,
            },
            {
              label: 'Outside field',
              before: { presence: 'present', value: 1 },
              after: { presence: 'present', value: 2 },
              recognized: false,
            },
          ],
        })}
        writesEnabled={false}
        dirty={false}
      />,
    );

    const table = screen.getByRole('table', { name: 'Suggested changes' });
    expect(table.querySelector('input, textarea, select')).toBeNull();
    expect(within(table).getByText('Not set')).toBeInTheDocument();
    expect(within(table).getByText('0')).toBeInTheDocument();
    expect(within(table).getByText('Null')).toBeInTheDocument();
    expect(within(table).getByText('False')).toBeInTheDocument();
    expect(within(table).getByText('Empty')).toBeInTheDocument();
    expect(within(table).getByText(/outside the current suggestion list/i)).toBeInTheDocument();
  });

  it('does not mark a suggestion adopted when the action is only requested', async () => {
    const onAdopt = vi.fn();
    const onReject = vi.fn();
    render(
      <ProposalReview
        proposal={createPendingProposal()}
        writesEnabled
        dirty={false}
        onAdopt={onAdopt}
        onReject={onReject}
      />,
    );
    const user = userEvent.setup();
    await user.click(screen.getByRole('button', { name: 'Adopt suggestion' }));
    await user.click(screen.getByRole('button', { name: 'Reject suggestion' }));
    expect(onAdopt).toHaveBeenCalledTimes(1);
    expect(onReject).toHaveBeenCalledTimes(1);
    expect(screen.getByRole('status')).toHaveTextContent(
      'Suggestion waiting for your review. It is not applied to the draft.',
    );
    expect(screen.queryByText('Suggestion adopted.')).not.toBeInTheDocument();
    expect(screen.queryByText('Suggestion rejected. The draft was not changed.')).not.toBeInTheDocument();
  });

  it.each([
    ['expired', 'This suggestion has expired and was not adopted.'],
    ['conflict', 'The draft changed before this action completed. The suggestion was not adopted.'],
    ['invalidated', 'This suggestion is no longer valid and was not adopted.'],
    ['unconfirmed', 'The last response was not confirmed. Refresh the server status before trying again. This page will not retry the action.'],
    ['adopting', 'Adoption is in progress. This page is waiting for the server result.'],
  ] as const)('keeps the %s state unsuccessful', async (state, copy) => {
    const onAdopt = vi.fn();
    const onRefresh = vi.fn();
    render(
      <ProposalReview
        proposal={createPendingProposal({ state })}
        writesEnabled
        dirty={false}
        onAdopt={onAdopt}
        onRefresh={onRefresh}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(copy);
    expect(screen.queryByText('Suggestion adopted.')).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Adopt suggestion' })).toBeDisabled();
    await userEvent.setup().click(
      screen.getByRole('button', { name: 'Adopt suggestion' }),
    );
    expect(onAdopt).not.toHaveBeenCalled();
    if (state === 'unconfirmed') {
      await userEvent.setup().click(
        screen.getByRole('button', { name: 'Refresh server status' }),
      );
      expect(onRefresh).toHaveBeenCalledTimes(1);
      expect(onAdopt).not.toHaveBeenCalled();
    }
  });

  it('points both actions at the rendered shared explanation', () => {
    render(
      <ProposalReview
        proposal={createPendingProposal()}
        writesEnabled={false}
        dirty={false}
      />,
    );
    const adopt = screen.getByRole('button', { name: 'Adopt suggestion' });
    const reject = screen.getByRole('button', { name: 'Reject suggestion' });
    const descriptionId = adopt.getAttribute('aria-describedby');
    expect(descriptionId).toBeTruthy();
    expect(reject).toHaveAttribute('aria-describedby', descriptionId);
    expect(document.getElementById(descriptionId ?? '')).toHaveTextContent(
      /reviewed configuration-draft interface/i,
    );
  });

  it('shows an adopted suggestion only when the server state says so', () => {
    render(
      <ProposalReview
        proposal={createPendingProposal({ state: 'adopted' })}
        writesEnabled
        dirty={false}
      />,
    );
    expect(screen.getByText('Suggestion adopted.')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Adopt suggestion' })).toBeDisabled();
  });
});

describe('explanation notice', () => {
  it('shows explanation text without adopt or reject actions', () => {
    const { rerender } = render(
      <ExplanationNotice
        explanation={{
          contentAvailable: true,
          labels: ['Example needs a sample selection.'],
        }}
      />,
    );
    expect(screen.getByText('Example needs a sample selection.')).toBeInTheDocument();
    expect(screen.getByText(/does not change the draft/i)).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: /adopt|reject/i })).not.toBeInTheDocument();

    rerender(<ExplanationNotice explanation={{ contentAvailable: false, labels: [] }} />);
    expect(screen.getByRole('status')).toHaveTextContent(
      'The explanation was not included in the response.',
    );
    expect(screen.queryByRole('button', { name: /adopt|reject/i })).not.toBeInTheDocument();
    expect(screen.queryByText('Suggestion adopted.')).not.toBeInTheDocument();
  });
});
