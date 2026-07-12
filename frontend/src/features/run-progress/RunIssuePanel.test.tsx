import { render, screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';
import { RunIssuePanel } from './RunIssuePanel';

describe('RunIssuePanel', () => {
  it('shows public issue guidance without technical details or context', () => {
    render(
      <RunIssuePanel
        issues={[
          {
            code: 'RUN_FAILED',
            message: 'The run failed.',
            severity: 'error',
            path: 'execution',
            source: 'worker',
            technical_message: 'sqlite at /private/runtime.sqlite raised an exception',
            hint: 'Review stderr and retry.',
            context: { private_path: '/private/runtime.sqlite' },
          },
        ]}
      />,
    );

    expect(screen.getByRole('alert')).toHaveTextContent('RUN_FAILED: The run failed.');
    expect(screen.getByRole('alert')).toHaveTextContent('Review stderr and retry.');
    expect(screen.queryByText(/private\/runtime/)).not.toBeInTheDocument();
    expect(screen.queryByText(/sqlite at/)).not.toBeInTheDocument();
  });
});
