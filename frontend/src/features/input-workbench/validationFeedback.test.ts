import type { RJSFSchema } from '@rjsf/utils';
import { afterEach, describe, expect, it } from 'vitest';
import { configIssueErrors, configIssueLocation, focusConfigIssue } from './validationFeedback';

const schema: RJSFSchema = {
  type: 'object',
  properties: { processing: { type: 'object', properties: { method: { type: 'string' } } } },
};

afterEach(() => { document.body.innerHTML = ''; });

describe('schema-based validation location', () => {
  it('falls back to known sections for unknown fields without interpreting hint text', () => {
    expect(configIssueLocation(schema, 'config.processing.missing')).toEqual({ path: 'config.processing', segments: ['processing'] });
    expect(configIssueLocation(schema, 'samples[0].sample')).toBeNull();
    expect(configIssueErrors(schema, [{ code: 'ANY', path: 'config.processing.method', message: 'Check value.', hint: 'config.other.field' }]))
      .toEqual({ processing: { method: { __errors: ['Check value. config.other.field'] } } });
    expect(configIssueErrors(schema, [{ code: 'NOTE', severity: 'warning', path: 'config.processing.method', message: 'Optional.' }])).toEqual({});
    const reserved = JSON.parse('{"properties":{"__proto__":{"type":"object"},"__errors":{"type":"string"}}}');
    expect(configIssueLocation(reserved, 'config.__proto__.polluted')?.segments).toEqual([]);
    expect(configIssueLocation(reserved, 'config.__errors')?.segments).toEqual([]);
  });

  it.each(['disabled', 'hidden'])('focuses a visible containing legend for a %s control', (attribute) => {
    document.body.innerHTML = `<section tabindex="-1"><fieldset id="root_processing"><legend>Processing</legend><input id="root_processing_method" ${attribute}></fieldset></section>`;
    focusConfigIssue(document.querySelector('section')!, schema, 'config.processing.method');
    expect(document.activeElement).toBe(document.querySelector('legend'));
    expect(document.querySelector('input')).toHaveAttribute(attribute);
  });

  it('opens a containing disclosure and scopes lookup to the config editor', () => {
    document.body.innerHTML = '<input id="root_processing_method"><section><details><summary>Processing</summary><input id="root_processing_method"></details></section>';
    const container = document.querySelector('section')!;
    focusConfigIssue(container, schema, 'config.processing.method');
    expect(document.querySelector('details')?.open).toBe(true);
    expect(document.activeElement).toBe(container.querySelector('input'));
  });
});
