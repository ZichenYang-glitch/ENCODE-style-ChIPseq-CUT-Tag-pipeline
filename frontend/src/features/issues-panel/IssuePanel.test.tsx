import { describe, it, expect, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { IssuePanel } from './IssuePanel';
import type { Issue } from '../../api/types';

const sampleIssue: Issue = {
  code: 'ENCODE_SAMPLES_INVALID',
  message: 'Sample sheet is invalid.',
  severity: 'error',
  path: 'samples',
  source: 'samples',
  technical_message: 'Missing samples field.',
  hint: 'Provide a path to a TSV sample sheet.',
  context: {},
};

describe('IssuePanel', () => {
  it('calls onAskAgent with the issue when Ask Agent is clicked', async () => {
    const onAskAgent = vi.fn();
    const user = userEvent.setup();

    render(
      <IssuePanel issues={[sampleIssue]} onAskAgent={onAskAgent} />,
    );

    const askButton = screen.getByRole('button', {
      name: /Ask Agent about ENCODE_SAMPLES_INVALID/i,
    });
    await user.click(askButton);

    expect(onAskAgent).toHaveBeenCalledTimes(1);
    expect(onAskAgent).toHaveBeenCalledWith(sampleIssue);
  });

  it('does not render Ask Agent button when onAskAgent is omitted', () => {
    render(<IssuePanel issues={[sampleIssue]} />);

    expect(
      screen.queryByRole('button', { name: /Ask Agent/i }),
    ).not.toBeInTheDocument();
  });

  it('groups issues by source', () => {
    const issues: Issue[] = [
      sampleIssue,
      {
        code: 'UI_ERROR',
        message: 'UI issue.',
        severity: 'warning',
        path: null,
        source: 'ui',
        technical_message: null,
        hint: null,
        context: {},
      },
    ];

    render(<IssuePanel issues={issues} />);

    expect(screen.getByRole('heading', { name: 'samples' })).toBeInTheDocument();
    expect(screen.getByRole('heading', { name: 'ui' })).toBeInTheDocument();
  });
});
