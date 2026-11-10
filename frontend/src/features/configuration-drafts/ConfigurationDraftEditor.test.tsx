import { useCallback, useState } from 'react';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import { ConfigurationDraftEditor } from './ConfigurationDraftEditor';
import { createDraftEditorSchema } from './testFixtures';

vi.mock('@uiw/react-codemirror', () => ({
  default: () => <textarea aria-label="Advanced config YAML" />,
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
    await waitFor(() => {
      expect(screen.queryByText(/you have unsaved edits/i)).not.toBeInTheDocument();
    });
    expect(onSave).toHaveBeenCalledTimes(1);
    expect(screen.getByLabelText('Local limit')).toHaveValue(1);
  });
});
