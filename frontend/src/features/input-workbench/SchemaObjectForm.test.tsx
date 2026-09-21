import { useState } from 'react';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { RJSFSchema } from '@rjsf/utils';
import { describe, expect, it, vi } from 'vitest';
import type { ValidationRequestConfig } from '../../api/generated/models';
import { SchemaObjectForm } from './SchemaObjectForm';
import { createDefaultObject } from './schemaContract';

const flagsSchema: RJSFSchema = {
  type: 'object',
  properties: {
    name: { type: 'string', title: 'Draft name' },
    settings: {
      type: 'object',
      properties: {
        checks: {
          type: 'object',
          default: { enabled: true, first: true, second: true, optional: false },
          properties: {
            enabled: { type: 'boolean', title: 'Checks enabled', default: true },
            first: { type: 'boolean', title: 'First check', default: true },
            second: { type: 'boolean', title: 'Second check', default: true },
            optional: { type: 'boolean', title: 'Optional check', default: false },
          },
          additionalProperties: false,
        },
      },
    },
  },
};

function renderForm(schema: RJSFSchema, initial = createDefaultObject(schema)) {
  const changed = vi.fn();
  function ControlledForm() {
    const [value, setValue] = useState(initial);
    return (
      <SchemaObjectForm
        schema={schema}
        value={value}
        resetRevision={0}
        onChange={(next) => {
          changed(next);
          setValue(next as ValidationRequestConfig);
        }}
        ariaLabel="Test config form"
      />
    );
  }
  render(<ControlledForm />);
  return changed;
}

describe('schema-owned section switches', () => {
  it('clears nested boolean options on disable and keeps them off on re-enable', async () => {
    const changed = renderForm(flagsSchema);
    const user = userEvent.setup();
    const master = screen.getByRole('checkbox', { name: 'Checks enabled' });
    expect(screen.getByRole('checkbox', { name: 'First check' })).toBeChecked();
    await user.click(master);
    expect(changed).toHaveBeenLastCalledWith({
      settings: { checks: { enabled: false, first: false, second: false, optional: false } },
    });
    for (const name of ['First check', 'Second check', 'Optional check']) {
      expect(screen.getByRole('checkbox', { name })).not.toBeChecked();
      expect(screen.getByRole('checkbox', { name })).toBeDisabled();
    }
    await user.click(master);
    expect(changed).toHaveBeenLastCalledWith({
      settings: { checks: { enabled: true, first: false, second: false, optional: false } },
    });
    await user.click(screen.getByRole('checkbox', { name: 'First check' }));
    expect(changed).toHaveBeenLastCalledWith({
      settings: { checks: { enabled: true, first: true, second: false, optional: false } },
    });
  });

  it('restores a declared disabled-only default, clearing optional settings and nested resources', async () => {
    const schema: RJSFSchema = {
      type: 'object',
      properties: {
        processing: {
          type: 'object',
          default: { enabled: false },
          properties: {
            enabled: { type: 'boolean', title: 'Processing enabled', default: false },
            method: { type: 'string', title: 'Method', enum: ['first', 'second'] },
            keep_stats: { type: 'boolean', title: 'Keep stats' },
            resource: {
              type: 'object',
              properties: { path: { type: 'string', title: 'Resource path' } },
            },
          },
          required: ['enabled'],
          allOf: [{ if: { properties: { enabled: { const: true } } }, then: { required: ['method'] } }],
          additionalProperties: false,
        },
      },
    };
    const changed = renderForm(schema, {
      processing: { enabled: true, method: 'second', keep_stats: true, resource: { path: '/reference' } },
    });
    const user = userEvent.setup();
    const master = screen.getByRole('checkbox', { name: 'Processing enabled' });
    await user.click(master);
    expect(changed).toHaveBeenLastCalledWith({ processing: { enabled: false } });
    expect(screen.getByRole('combobox', { name: 'Method' })).toBeDisabled();
    expect(screen.getByRole('checkbox', { name: 'Keep stats' })).toBeDisabled();
    expect(screen.getByRole('textbox', { name: 'Resource path' })).toBeDisabled();
    await user.click(master);
    expect(changed).toHaveBeenLastCalledWith({ processing: { enabled: true } });
    expect(screen.getByRole('combobox', { name: /Method/ })).toBeEnabled();
  });

  it('preserves required non-boolean settings when their master is disabled', async () => {
    const changed = renderForm({
      type: 'object',
      properties: {
        trimming: {
          type: 'object',
          default: { enabled: true, tool: 'fastp' },
          properties: {
            enabled: { type: 'boolean', title: 'Trimming enabled' },
            tool: { type: 'string', title: 'Tool', enum: ['fastp', 'trimgalore'] },
          },
          required: ['enabled', 'tool'],
        },
      },
    });
    await userEvent.setup().click(screen.getByRole('checkbox', { name: 'Trimming enabled' }));
    expect(changed).toHaveBeenLastCalledWith({ trimming: { enabled: false, tool: 'fastp' } });
  });

  it('leaves imported conflicts intact on mount and unrelated edits for normal validation', async () => {
    const changed = renderForm(flagsSchema, {
      settings: { checks: { enabled: false, first: true, second: false, optional: false } },
    });
    expect(changed).not.toHaveBeenCalled();
    await userEvent.setup().type(screen.getByRole('textbox', { name: 'Draft name' }), 'x');
    expect(changed).toHaveBeenLastCalledWith({
      name: 'x',
      settings: { checks: { enabled: false, first: true, second: false, optional: false } },
    });
  });
});
