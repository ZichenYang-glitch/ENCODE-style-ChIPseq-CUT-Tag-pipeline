import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { AssistantConnectionPanel } from './AssistantConnectionPanel';

describe('assistant connection panel', () => {
  it.each([
    ['disconnected', 'Assistant not connected'],
    ['connecting', 'Connecting the assistant'],
    ['available', 'Assistant available'],
    ['ended', 'Assistant connection ended'],
    ['bridge-unavailable', 'Assistant unavailable'],
  ] as const)('labels the %s state', (state, heading) => {
    render(
      <AssistantConnectionPanel
        connection={{ state, draftRestored: false }}
        writesEnabled={false}
      />,
    );
    expect(screen.getByText(heading)).toBeInTheDocument();
  });

  it('keeps a restored draft separate from an assistant connection', () => {
    render(
      <AssistantConnectionPanel
        connection={{ state: 'disconnected', draftRestored: true }}
        writesEnabled={false}
      />,
    );
    expect(screen.getByText('Assistant not connected')).toBeInTheDocument();
    expect(screen.getByText(/draft is still here after signing in/i)).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
    expect(
      screen.queryByRole('button', { name: /resume|fork|permission|provider/i }),
    ).not.toBeInTheDocument();
  });

  it('does not report a successful connection when the assistant is unavailable', async () => {
    const onConnect = vi.fn();
    render(
      <AssistantConnectionPanel
        connection={{ state: 'bridge-unavailable', draftRestored: true }}
        writesEnabled
        onConnect={onConnect}
      />,
    );
    expect(screen.getByText('Assistant unavailable')).toBeInTheDocument();
    expect(screen.getByText(/was not reconnected/i)).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
    const connect = screen.getByRole('button', { name: 'Connect assistant' });
    await userEvent.setup().click(connect);
    expect(onConnect).toHaveBeenCalledTimes(1);
    expect(screen.getByText('Assistant unavailable')).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
  });

  it('does not change the visible state when connect is requested', async () => {
    const onConnect = vi.fn();
    render(
      <AssistantConnectionPanel
        connection={{ state: 'disconnected', draftRestored: false }}
        writesEnabled
        onConnect={onConnect}
      />,
    );
    await userEvent.setup().click(
      screen.getByRole('button', { name: 'Connect assistant' }),
    );
    expect(onConnect).toHaveBeenCalledTimes(1);
    expect(screen.getByText('Assistant not connected')).toBeInTheDocument();
    expect(screen.queryByText('Assistant available')).not.toBeInTheDocument();
  });
});
