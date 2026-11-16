import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { ConfigurationDraftList } from './ConfigurationDraftList';
import { createDraftSummary } from './testFixtures';

describe('configuration draft list', () => {
  it('shows loading, empty, unavailable, and error without a success message', async () => {
    const retry = vi.fn();
    const { rerender } = render(
      <ConfigurationDraftList
        model={{
          status: 'loading',
          drafts: [],
          pageLimit: 20,
          hasPrevious: false,
          hasNext: false,
        }}
        selectedSid={null}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(
      'Loading your configuration drafts.',
    );

    rerender(
      <ConfigurationDraftList
        model={{
          status: 'empty',
          drafts: [],
          pageLimit: 20,
          hasPrevious: false,
          hasNext: false,
        }}
        selectedSid={null}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(
      'No configuration drafts for this workflow.',
    );

    rerender(
      <ConfigurationDraftList
        model={{
          status: 'unavailable',
          drafts: [],
          pageLimit: 20,
          hasPrevious: false,
          hasNext: false,
        }}
        selectedSid={null}
      />,
    );
    expect(screen.getByRole('status')).toHaveTextContent(
      'The draft list is not connected yet. Nothing was changed.',
    );

    rerender(
      <ConfigurationDraftList
        model={{
          status: 'error',
          drafts: [],
          pageLimit: 20,
          hasPrevious: false,
          hasNext: false,
          errorMessage: 'The draft list could not be loaded. Nothing was changed.',
        }}
        selectedSid={null}
        onRetry={retry}
      />,
    );
    expect(screen.getByRole('alert')).toHaveTextContent(
      'Nothing was changed.',
    );
    expect(screen.queryByText(/saved|adopted|validated/i)).not.toBeInTheDocument();
    await userEvent.setup().click(
      screen.getByRole('button', { name: 'Retry draft list' }),
    );
    expect(retry).toHaveBeenCalledTimes(1);
  });

  it('bounds the page and moves between drafts with the keyboard', async () => {
    const user = userEvent.setup();
    const onSelect = vi.fn();
    const onNext = vi.fn();
    render(
      <ConfigurationDraftList
        model={{
          status: 'ready',
          drafts: [
            createDraftSummary('one', 'Alpha project'),
            createDraftSummary('two', 'Beta project', 'expired'),
            createDraftSummary('three', 'Hidden project'),
          ],
          pageLimit: 2,
          hasPrevious: false,
          hasNext: true,
        }}
        selectedSid="one"
        onSelect={onSelect}
        onNext={onNext}
      />,
    );

    expect(screen.getByText('Up to 2 drafts per page.')).toBeInTheDocument();
    expect(screen.queryByText('Hidden project')).not.toBeInTheDocument();
    expect(screen.getByText(/Adoption receipts can still be reviewed/)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Alpha project/ })).toHaveAttribute(
      'aria-pressed',
      'true',
    );

    await user.tab();
    expect(screen.getByRole('button', { name: /Alpha project/ })).toHaveFocus();
    await user.keyboard('{ArrowDown}');
    expect(screen.getByRole('button', { name: /Beta project/ })).toHaveFocus();
    await user.keyboard('{Enter}');
    expect(onSelect).toHaveBeenCalledTimes(1);
    expect(onSelect).toHaveBeenCalledWith('two');

    await user.click(screen.getByRole('button', { name: 'Next page' }));
    expect(onNext).toHaveBeenCalledTimes(1);
    expect(screen.queryByText('Hidden project')).not.toBeInTheDocument();
  });

  it('caps the advertised page size at the product maximum', () => {
    render(
      <ConfigurationDraftList
        model={{
          status: 'ready',
          drafts: [createDraftSummary('one', 'Alpha project')],
          pageLimit: 80,
          hasPrevious: true,
          hasNext: false,
        }}
        selectedSid={null}
      />,
    );
    expect(screen.getByText('Up to 50 drafts per page.')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Previous page' })).toBeDisabled();
    expect(screen.getByRole('button', { name: 'Next page' })).toBeDisabled();
  });
});
