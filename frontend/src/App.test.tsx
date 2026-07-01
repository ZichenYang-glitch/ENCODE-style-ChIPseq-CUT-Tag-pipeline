import { describe, it, expect } from 'vitest';
import { render, screen } from '@testing-library/react';
import { App } from './App';

describe('App shell', () => {
  it('renders the workflow platform heading', async () => {
    render(<App />);
    expect(
      await screen.findByRole('heading', { name: /workflow platform/i }),
    ).toBeInTheDocument();
  });

  it('loads and shows the stub workflow', async () => {
    render(<App />);
    expect(
      await screen.findByText(/ENCODE-style ChIP-seq/i),
    ).toBeInTheDocument();
  });
});
