import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { SampleIdentity } from './SampleIdentity';

describe('SampleIdentity', () => {
  it('retains the unmodified value on copy failure without revealing the exception', async () => {
    const user = userEvent.setup();
    vi.spyOn(navigator.clipboard, 'writeText').mockRejectedValueOnce(new Error('/private/clipboard-path'));
    render(<SampleIdentity sampleId="a  b " />);
    await user.click(screen.getByRole('button', { name: 'Copy sample ID' }));
    expect(screen.getByRole('status')).toHaveTextContent('Could not copy sample ID.');
    expect(document.querySelector('[data-sample-label]')).toHaveTextContent('"a  b "', { normalizeWhitespace: false });
    expect(document.body).not.toHaveTextContent('/private/clipboard-path');
  });

  it('clears stale copied feedback when the exact sample identity changes', async () => {
    const user = userEvent.setup();
    const { rerender } = render(<SampleIdentity sampleId="a b" />);
    await user.click(screen.getByRole('button', { name: 'Copy sample ID' }));
    expect(await navigator.clipboard.readText()).toBe('a b');
    expect(screen.getByRole('status')).toHaveTextContent('Sample ID copied.');
    rerender(<SampleIdentity sampleId="a b " />);
    expect(screen.queryByRole('status')).not.toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Copy sample ID' }));
    expect(await navigator.clipboard.readText()).toBe('a b ');
  });
});
