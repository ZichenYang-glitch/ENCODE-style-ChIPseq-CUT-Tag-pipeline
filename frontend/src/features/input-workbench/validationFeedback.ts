import type { ErrorSchema, RJSFSchema } from '@rjsf/utils';
import { isPlainObject } from './jsonSafety';

export interface SafeIssue {
  code: string;
  message: string;
  severity?: 'error' | 'warning' | 'info';
  path?: string | null;
  hint?: string | null;
}

export interface ValidationFeedback {
  revision: number;
  issues: SafeIssue[];
  notice?: string;
  changedInFlight?: boolean;
  snapshotId?: string;
}

export interface ConfigFocusRequest {
  path: string;
  sequence: number;
}

// Paths come from the served schema, never from hint text or an adapter code table.
export function configIssueLocation(schema: RJSFSchema, path?: string | null) {
  if (path !== 'config' && !path?.startsWith('config.')) return null;
  const segments: string[] = [];
  let node = schema;
  let matched = 'config';
  while (path !== matched) {
    const matches = Object.entries(node.properties ?? {}).filter(([key, candidate]) =>
      !['__proto__', 'prototype', 'constructor', '__errors'].includes(key) &&
      isPlainObject(candidate) &&
      (path === `${matched}.${key}` || path?.startsWith(`${matched}.${key}.`) ||
        path?.startsWith(`${matched}.${key}[`)),
    );
    // An ambiguous or conditional path falls back to its known containing section.
    if (matches.length !== 1) break;
    const [key, candidate] = matches[0];
    segments.push(key);
    matched += `.${key}`;
    node = candidate as RJSFSchema;
  }
  return { path: matched, segments };
}

export function configIssueErrors(schema: RJSFSchema, issues: SafeIssue[]): ErrorSchema {
  const result: ErrorSchema = {};
  for (const issue of issues) {
    if (issue.severity === 'warning' || issue.severity === 'info') continue;
    const location = configIssueLocation(schema, issue.path);
    if (!location) continue;
    let node = result;
    for (const key of location.segments) {
      const child = Object.prototype.hasOwnProperty.call(node, key) ? node[key] : undefined;
      node[key] = isPlainObject(child) ? child : {};
      node = node[key] as ErrorSchema;
    }
    node.__errors = [...(node.__errors ?? []), [issue.message, issue.hint].filter(Boolean).join(' ')];
  }
  return result;
}

export function focusConfigIssue(container: HTMLElement, schema: RJSFSchema, path: string) {
  const location = configIssueLocation(schema, path);
  const segments = location?.segments.slice() ?? [];
  let target: HTMLElement | undefined;
  for (;;) {
    const id = ['root', ...segments].join('_');
    const matches = Array.from(container.querySelectorAll<HTMLElement>('[id]'))
      .filter((element) => element.id === id);
    if (matches.length === 1) {
      const element = matches[0];
      const hidden = (node: HTMLElement): boolean => {
        for (let item: HTMLElement | null = node; item && item !== container; item = item.parentElement) {
          if (item.hidden || item.getAttribute('aria-hidden') === 'true' ||
            getComputedStyle(item).display === 'none' || getComputedStyle(item).visibility === 'hidden') return true;
        }
        return false;
      };
      if (!element.matches(':disabled') && !hidden(element)) {
        target = element.tagName === 'FIELDSET'
          ? element.querySelector<HTMLElement>(':scope > legend') ?? element
          : element;
        break;
      }
    }
    if (segments.length === 0) break;
    segments.pop();
  }
  target ??= container;
  for (let parent = target.parentElement; parent; parent = parent.parentElement) {
    if (parent instanceof HTMLDetailsElement) parent.open = true;
  }
  if (!target.matches('input, select, textarea, button, a[href], [contenteditable="true"]')) {
    target.tabIndex = -1;
  }
  target.focus();
  target.scrollIntoView?.({ block: 'center' });
}
