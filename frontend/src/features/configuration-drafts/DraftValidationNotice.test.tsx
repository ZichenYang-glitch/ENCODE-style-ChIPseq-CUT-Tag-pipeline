import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { DraftValidationNotice } from './DraftValidationNotice';

describe('draft validation notice', () => {
  it.each([
    ['unavailable', 'Draft validation is not connected yet.'],
    ['unconfirmed', 'The validation response was not confirmed. This page has no validation result.'],
    ['issues', 'Validation reported problems.'],
  ] as const)('does not present %s as a finished validation', (status, copy) => {
    render(
      <DraftValidationNotice
        validation={{ status, summary: 'The server reported a problem.' }}
        writesEnabled
        dirty={false}
        draftExpired={false}
        onValidate={vi.fn()}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(copy);
    expect(screen.queryByText(/validation finished/i)).not.toBeInTheDocument();
    expect(screen.queryByText(/no execution snapshot/i)).not.toBeInTheDocument();
  });

  it('describes a completed check without an execution snapshot', () => {
    render(
      <DraftValidationNotice
        validation={{ status: 'no-snapshot' }}
        writesEnabled={false}
        dirty={false}
        draftExpired={false}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(
      'The configuration check passed, and there is no execution snapshot.',
    );
    expect(screen.queryByText(/validation finished/i)).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Validate draft' })).toBeDisabled();
  });

  it('shows a finished validation without offering to start a run', () => {
    render(
      <DraftValidationNotice
        validation={{
          status: 'snapshot',
          revisionLabel: 'Revision 2',
          expiresLabel: '17 Nov 2026',
        }}
        writesEnabled={false}
        dirty={false}
        draftExpired={false}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(
      'Validation finished. This does not start a run.',
    );
    expect(screen.getByText('Revision 2')).toBeInTheDocument();
    expect(screen.getByText('17 Nov 2026')).toBeInTheDocument();
    expect(screen.queryByRole('link', { name: /run|preflight|start/i })).not.toBeInTheDocument();
  });

  it('does not record success when validate is clicked', async () => {
    const onValidate = vi.fn();
    render(
      <DraftValidationNotice
        validation={{ status: 'unavailable' }}
        writesEnabled
        dirty={false}
        draftExpired={false}
        onValidate={onValidate}
      />,
    );
    await userEvent.setup().click(
      screen.getByRole('button', { name: 'Validate draft' }),
    );
    expect(onValidate).toHaveBeenCalledTimes(1);
    expect(screen.getByRole('status')).toHaveTextContent(
      'Draft validation is not connected yet.',
    );
    expect(screen.queryByText(/validation finished/i)).not.toBeInTheDocument();
  });
});
