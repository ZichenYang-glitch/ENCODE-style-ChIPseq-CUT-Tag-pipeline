import type { RJSFSchema } from '@rjsf/utils';
import { isPlainObject } from './jsonSafety';

export function sectionDisableMode(schema: RJSFSchema): 'reset' | 'flags' | null {
  const properties = schema.properties ?? {};
  if (
    schema.type !== 'object' ||
    !isPlainObject(properties.enabled) ||
    properties.enabled.type !== 'boolean' ||
    !Object.keys(properties).some((key) => key !== 'enabled')
  ) {
    return null;
  }
  // A disabled-only object default declares the section's neutral state.
  if (
    isPlainObject(schema.default) &&
    schema.default.enabled === false &&
    Object.keys(schema.default).length === 1 &&
    (schema.required ?? []).every((key) => key === 'enabled')
  ) {
    return 'reset';
  }
  return Object.entries(properties).some(
    ([key, property]) =>
      key !== 'enabled' && isPlainObject(property) && property.type === 'boolean',
  )
    ? 'flags'
    : null;
}

export function cascadeSectionSwitches(
  schema: RJSFSchema,
  previous: unknown,
  next: unknown,
): unknown {
  if (!isPlainObject(next)) return next;
  const result = { ...next };
  const properties = schema.properties ?? {};
  const mode = sectionDisableMode(schema);
  if (
    mode &&
    isPlainObject(previous) &&
    previous.enabled === true &&
    next.enabled === false
  ) {
    for (const [key, property] of Object.entries(properties)) {
      if (key === 'enabled') continue;
      if (mode === 'reset') {
        delete result[key];
      } else if (isPlainObject(property) && property.type === 'boolean') {
        result[key] = false;
      }
    }
    // Unknown keys remain available to the existing validator.
    return result;
  }
  for (const [key, property] of Object.entries(properties)) {
    if (
      !isPlainObject(property) ||
      !Object.prototype.hasOwnProperty.call(next, key)
    ) {
      continue;
    }
    result[key] = cascadeSectionSwitches(
      property as RJSFSchema,
      isPlainObject(previous) ? previous[key] : undefined,
      next[key],
    );
  }
  return result;
}
