import { render, screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';
import { ViewerHost } from './ViewerHost';

describe('ViewerHost', () => {
  it('embeds the same-origin viewer and opens the full page without a wildcard target', () => {
    render(<ViewerHost runId="run-1" />);
    const frame = screen.getByTitle('Genome browser');
    expect(frame).toHaveAttribute('src', expect.stringMatching(/^\/viewer\.html#instance=[0-9a-f-]{36}$/i));
    expect(frame.getAttribute('src')).not.toContain('*');
    const fullPage = screen.getByRole('link', { name: 'Open genome browser full page' });
    expect(fullPage).toHaveAttribute('href', '/viewer.html?run=run-1');
    expect(fullPage).toHaveAttribute('rel', 'noopener noreferrer');
  });

  it('ignores a ready notice from a retired instance', () => {
    render(<ViewerHost runId="run-1" />);
    const frame = screen.getByTitle('Genome browser') as HTMLIFrameElement;
    const retired = '22222222-2222-4222-8222-222222222222';
    window.dispatchEvent(
      new MessageEvent('message', {
        data: {
          version: 1,
          type: 'helixweave.viewer.ready',
          instanceId: retired,
          payload: {
            runId: 'run-1',
            referenceStatus: 'bound',
            artifactIds: ['artifact-peak'],
          },
        },
        origin: window.location.origin,
        source: frame.contentWindow,
      }),
    );
    expect(screen.getByRole('status')).toHaveTextContent('Opening genome browser.');
  });
});
