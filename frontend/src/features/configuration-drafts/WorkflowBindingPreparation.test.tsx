import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { describe, expect, it, vi } from 'vitest';
import type { SampleColumn } from '../input-workbench/schemaContract';
import {
  acceptCatalogPage,
  WorkflowBindingPreparation,
} from './WorkflowBindingPreparation';

const columns: SampleColumn[] = [
  {
    key: 'layout',
    label: 'layout',
    description: null,
    required: true,
    enumValues: ['PE', 'SE'],
    defaultValue: '',
    readOnly: false,
  },
  {
    key: 'sample',
    label: 'sample',
    description: 'Sample identifier.',
    required: true,
    enumValues: null,
    defaultValue: '',
    readOnly: false,
  },
];

describe('workflow binding preparation', () => {
  it('does not offer a path field when the catalog is not connected', () => {
    render(
      <WorkflowBindingPreparation
        {...baseProps()}
        catalogStatus="unavailable"
      />,
    );
    expect(screen.getByText(/not connected yet/i)).toBeInTheDocument();
    expect(screen.getByLabelText('layout')).toBeInTheDocument();
    expect(screen.queryByRole('textbox')).not.toBeInTheDocument();
    expect(document.querySelector('input[type="file"]')).toBeNull();
    const sampleField = screen.getByText(/^sample\./).closest('p');
    expect(sampleField).toHaveTextContent('This field is not connected to save yet.');
    expect(sampleField).not.toHaveTextContent(/input revision/i);
  });

  it('reports an empty catalog and a failed catalog separately', () => {
    const { rerender } = render(
      <WorkflowBindingPreparation {...baseProps()} catalogStatus="empty" />,
    );
    expect(screen.getByText(/no admitted input revisions/i)).toBeInTheDocument();
    rerender(
      <WorkflowBindingPreparation
        {...baseProps()}
        catalogStatus="error"
        catalogError="The catalog could not be loaded."
      />,
    );
    expect(screen.getByRole('alert')).toHaveTextContent('The catalog could not be loaded.');
  });

  it('keeps a control or input revision visible when it leaves the list', async () => {
    const user = userEvent.setup();
    render(
      <WorkflowBindingPreparation
        {...baseProps()}
        catalogStatus="ready"
        catalogItems={[{ revisionId: 'listed', label: 'Listed read' }]}
        selectedInputRevisionIds={['listed', 'gone']}
        unlistedInputCount={1}
        controlSampleRevisionId="smpr_missing"
      />,
    );
    expect(screen.getByText(/control sample revision is no longer listed/i)).toBeInTheDocument();
    expect(screen.getByText(/input revision is no longer listed/i)).toBeInTheDocument();
    await user.click(screen.getByRole('checkbox', { name: 'Listed read' }));
  });

  it('pages the catalog and ignores a page from an older request', async () => {
    const user = userEvent.setup();
    const onNext = vi.fn();
    render(
      <WorkflowBindingPreparation
        {...baseProps()}
        catalogStatus="ready"
        catalogItems={[{ revisionId: 'one', label: 'Read one' }]}
        hasNext
        onNext={onNext}
      />,
    );
    await user.click(screen.getByRole('button', { name: 'Next input revisions' }));
    expect(onNext).toHaveBeenCalledTimes(1);
    const first = { label: 'Read one' };
    const late = { label: 'Read from the previous page' };
    expect(acceptCatalogPage(1, 2, late)).toBeNull();
    expect(acceptCatalogPage(2, 2, first)).toEqual(first);
  });
});

function baseProps() {
  return {
    columns,
    rowChoices: {},
    onChangeChoice: vi.fn(),
    samples: [],
    controlSampleRevisionId: null,
    onSelectControl: vi.fn(),
    catalogStatus: 'unavailable' as const,
    catalogError: null,
    catalogItems: [],
    selectedInputRevisionIds: [],
    unlistedInputCount: 0,
    onToggleInput: vi.fn(),
    hasPrevious: false,
    hasNext: false,
  };
}
