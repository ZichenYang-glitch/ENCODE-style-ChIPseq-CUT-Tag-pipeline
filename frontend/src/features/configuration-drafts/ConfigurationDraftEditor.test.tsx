import { useCallback, useState } from 'react';
import { fireEvent, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { ConfigurationDraftEditor } from './ConfigurationDraftEditor';
import { createDraftEditorSchema } from './testFixtures';

vi.mock('@uiw/react-codemirror', () => ({
  default: ({
    value,
    onChange,
  }: {
    value: string;
    onChange: (text: string) => void;
  }) => (
    <textarea
      aria-label="Advanced config YAML"
      value={value}
      onChange={(event) => onChange(event.target.value)}
    />
  ),
}));

const schema = createDraftEditorSchema();

describe('configuration draft editor', () => {
  it('does not save or mark unsaved edits while the form initializes', async () => {
    const onSave = vi.fn();
    const onDiscard = vi.fn();
    render(
      <ConfigurationDraftEditor
        schema={schema}
        schemaIdentity="schema-a"
        discardEpoch={0}
        writesEnabled
        onDirtyChange={() => undefined}
        onSave={onSave}
        onDiscard={onDiscard}
      />,
    );

    expect(await screen.findByLabelText('Local limit')).toHaveValue(1);
    await waitFor(() => {
      expect(screen.queryByText(/unsaved edits/i)).not.toBeInTheDocument();
    });
    expect(onSave).not.toHaveBeenCalled();
    expect(onDiscard).not.toHaveBeenCalled();
  });

  it('protects a user edit without treating save as confirmed', async () => {
    const onSave = vi.fn();
    function Harness() {
      const [epoch, setEpoch] = useState(0);
      const [identity, setIdentity] = useState('schema-a');
      const onDirtyChange = useCallback(() => undefined, []);
      return (
        <>
          <button type="button" onClick={() => setIdentity('schema-b')}>
            Refresh schema
          </button>
          <ConfigurationDraftEditor
            schema={schema}
            schemaIdentity={identity}
            discardEpoch={epoch}
            writesEnabled
            onDirtyChange={onDirtyChange}
            onSave={onSave}
            onDiscard={() => setEpoch((value) => value + 1)}
          />
        </>
      );
    }
    const user = userEvent.setup();
    render(<Harness />);
    const limit = await screen.findByLabelText('Local limit');
    await user.clear(limit);
    await user.type(limit, '0');
    expect(await screen.findByText(/you have unsaved edits/i)).toBeInTheDocument();
    expect(onSave).not.toHaveBeenCalled();

    await user.click(screen.getByRole('button', { name: 'Save edits' }));
    expect(onSave).toHaveBeenCalledTimes(1);
    expect(screen.getByText(/you have unsaved edits/i)).toBeInTheDocument();
    expect(screen.queryByText(/edits saved|suggestion adopted|validation finished/i)).not.toBeInTheDocument();

    await user.click(screen.getByRole('button', { name: 'Discard edits' }));
    await waitFor(() => {
      expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
    });
    expect(onSave).toHaveBeenCalledTimes(1);

    await user.clear(await screen.findByLabelText('Local limit'));
    await user.type(screen.getByLabelText('Local limit'), '4');
    expect(await screen.findByText(/you have unsaved edits/i)).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Refresh schema' }));
    expect(await screen.findByLabelText('Local limit')).toHaveValue(4);
    expect(screen.getByText(/you have unsaved edits/i)).toBeInTheDocument();
    expect(
      screen.getByText(/these edits stay here until you discard them/i),
    ).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Save edits' })).toBeDisabled();

    await user.click(screen.getByRole('button', { name: 'Keep these edits' }));
    expect(screen.getByLabelText('Local limit')).toHaveValue(4);
    expect(screen.getByText(/remain on the previous contract/i)).toBeInTheDocument();

    await user.click(screen.getByRole('button', { name: 'Discard edits and reload' }));
    await waitFor(() => {
      expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
    });
    expect(
      screen.queryByText(/these edits stay here until you discard them/i),
    ).not.toBeInTheDocument();
    expect(
      screen.queryByText(/discard these edits and reload before saving/i),
    ).not.toBeInTheDocument();
    expect(screen.getByLabelText('Local limit')).toHaveValue(1);
    expect(onSave).toHaveBeenCalledTimes(1);
  });

  it('reloads a clean editor when the schema identity changes', async () => {
    function Harness() {
      const [identity, setIdentity] = useState('schema-a');
      const onDirtyChange = useCallback(() => undefined, []);
      return (
        <>
          <button type="button" onClick={() => setIdentity('schema-b')}>
            Refresh schema
          </button>
          <ConfigurationDraftEditor
            schema={schema}
            schemaIdentity={identity}
            discardEpoch={0}
            writesEnabled
            onDirtyChange={onDirtyChange}
            onDiscard={() => undefined}
          />
        </>
      );
    }
    const user = userEvent.setup();
    render(<Harness />);
    expect(await screen.findByLabelText('Local limit')).toHaveValue(1);
    await user.click(screen.getByRole('button', { name: 'Refresh schema' }));
    await waitFor(() => {
      expect(screen.queryByText(/editor contract changed/i)).not.toBeInTheDocument();
    });
    expect(screen.getByLabelText('Local limit')).toHaveValue(1);
    expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
  });

  it('does not treat a mode switch as an edit', async () => {
    const onDirtyChange = vi.fn();
    render(
      <ConfigurationDraftEditor
        schema={schema}
        schemaIdentity="schema-a"
        discardEpoch={0}
        writesEnabled
        onDirtyChange={onDirtyChange}
        onDiscard={() => undefined}
      />,
    );
    const user = userEvent.setup();
    await screen.findByLabelText('Local limit');
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    expect(screen.getByLabelText('Advanced config YAML')).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Form mode' }));
    expect(screen.getByLabelText('Local limit')).toHaveValue(1);
    expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
    expect(onDirtyChange).not.toHaveBeenCalledWith(true);
  });

  it('protects invalid YAML until it is fixed or discarded', async () => {
    const onDiscard = vi.fn();
    function Harness() {
      const [epoch, setEpoch] = useState(0);
      return (
        <ConfigurationDraftEditor
          schema={schema}
          schemaIdentity="schema-a"
          discardEpoch={epoch}
          writesEnabled
          onDirtyChange={() => undefined}
          onSave={vi.fn()}
          onDiscard={() => {
            onDiscard();
            setEpoch((value) => value + 1);
          }}
        />
      );
    }
    const user = userEvent.setup();
    render(<Harness />);
    await screen.findByLabelText('Local limit');
    await user.click(screen.getByRole('button', { name: 'YAML mode' }));
    const input = screen.getByLabelText('Advanced config YAML');
    fireEvent.pointerDown(input);
    fireEvent.change(input, { target: { value: 'limit: [' } });
    expect(await screen.findByRole('alert')).toBeInTheDocument();
    expect(input).toHaveValue('limit: [');
    expect(screen.getByText(/you have unsaved edits/i)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Save edits' })).toBeDisabled();
    expect(screen.getByText(/fix the editor errors before saving/i)).toBeInTheDocument();

    fireEvent.change(input, { target: { value: 'limit: 1\n' } });
    await waitFor(() => {
      expect(screen.queryByRole('alert')).not.toBeInTheDocument();
    });
    expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();

    fireEvent.change(input, { target: { value: 'limit: [' } });
    expect(await screen.findByText(/you have unsaved edits/i)).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Discard edits' }));
    await waitFor(() => {
      expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
    });
    expect(onDiscard).toHaveBeenCalledTimes(1);
    expect(screen.queryByDisplayValue('limit: [')).not.toBeInTheDocument();
  });
});
