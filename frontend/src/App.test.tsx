import { describe, it, expect } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';

describe('App shell', () => {
  it('renders the workflow platform heading', async () => {
    render(<App />);
    expect(
      await screen.findByRole('heading', { name: /workflow platform/i }),
    ).toBeInTheDocument();
  });

  it('shows the stub workflow catalog', async () => {
    render(<App />);
    expect(
      await screen.findByText(/ENCODE-style ChIP-seq/i),
    ).toBeInTheDocument();
  });

  it('loads schema hints when a workflow is selected', async () => {
    const user = userEvent.setup();
    render(<App />);
    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);
    expect(await screen.findByText(/Config schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Sample schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Options schema/i)).toBeInTheDocument();
  });

  it('renders stub validation issues when Validate is clicked', async () => {
    const user = userEvent.setup();
    render(<App />);
    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);
    const validateButton = await screen.findByTestId('validate-button');
    await user.click(validateButton);
    expect(await screen.findByText(/Sample sheet is invalid/i)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Validate/i })).toBeInTheDocument();
  });

  it('renders a frontend parse error for invalid JSON and does not require a backend', async () => {
    const user = userEvent.setup();
    render(<App />);
    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);

    const configInput = screen.getByLabelText(/Config \(JSON\)/i);
    await user.clear(configInput);
    await user.type(configInput, 'not json');

    const validateButton = screen.getByTestId('validate-button');
    await user.click(validateButton);

    expect(await screen.findByText(/FRONTEND_INPUT_INVALID/i)).toBeInTheDocument();
    expect(screen.getByText(/Invalid JSON for config/i)).toBeInTheDocument();
  });

  it('renders the read-only agent sidebar label', async () => {
    render(<App />);
    expect(
      await screen.findByText(/Validation Assistant — Read Only/i),
    ).toBeInTheDocument();
  });
});
