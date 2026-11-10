import { useId, useState, type FormEvent } from 'react';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { writeBlockReason } from './actionGate';
import type {
  AuthoringCapabilityView,
  CapabilityField,
  StructuredQuestion,
  SuggestionIntent,
  SuggestionRequestStatus,
} from './model';

interface StructuredRequestPanelProps {
  capability: AuthoringCapabilityView | null;
  connectionAvailable: boolean;
  draftExpired: boolean;
  writesEnabled: boolean;
  dirty: boolean;
  requestStatus: SuggestionRequestStatus;
  onRequest?: (question: StructuredQuestion) => void;
}

export function StructuredRequestPanel({
  capability,
  connectionAvailable,
  draftExpired,
  writesEnabled,
  dirty,
  requestStatus,
  onRequest,
}: StructuredRequestPanelProps) {
  const baseId = useId();
  const [intent, setIntent] = useState<SuggestionIntent | null>(null);
  const [selected, setSelected] = useState<ReadonlySet<string>>(new Set());
  const [numbers, setNumbers] = useState<Record<string, string>>({});
  const [texts, setTexts] = useState<Record<string, string>>({});
  const [flags, setFlags] = useState<Record<string, boolean>>({});
  const [issueCode, setIssueCode] = useState<string | null>(null);
  const [formError, setFormError] = useState<string | null>(null);

  if (!capability || !capability.available) {
    return (
      <Panel title="Suggestion request">
        <p role="status">
          {!capability
            ? 'Suggestion fields are not connected yet. You can still edit by hand.'
            : capability.unavailableReason ??
              'This workflow has no suggestion fields. You can still edit by hand.'}
        </p>
        <p className="mt-2 text-xs text-[var(--color-text-muted)]">
          Free-text questions, file uploads, and pasted logs are not part of
          this review.
        </p>
      </Panel>
    );
  }

  const blockReason = requestBlockReason({
    writesEnabled,
    dirty,
    connectionAvailable,
    draftExpired,
    requestStatus,
    onRequest: onRequest !== undefined,
  });

  function toggleField(token: string) {
    setSelected((current) => {
      const next = new Set(current);
      if (next.has(token)) next.delete(token);
      else next.add(token);
      return next;
    });
  }

  function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (blockReason !== null || !intent || !onRequest) return;
    const question = buildQuestion({
      intent,
      fields: capability?.fields ?? [],
      selected,
      numbers,
      texts,
      flags,
      issueCode,
    });
    if (!question.ok) {
      setFormError(question.message);
      return;
    }
    setFormError(null);
    onRequest(question.value);
  }

  return (
    <Panel title="Suggestion request">
      <form className="space-y-4" onSubmit={submit}>
        <p className="text-xs text-[var(--color-text-muted)]">
          Choose a declared intent and field. Free-text questions, file
          uploads, and pasted logs are not part of this review.
        </p>
        <fieldset className="space-y-2">
          <legend className="text-sm font-medium">Request type</legend>
          {capability.intents.map((candidate) => (
            <label key={candidate} className="flex min-h-11 items-center gap-2 text-sm">
              <input
                type="radio"
                name={`${baseId}-intent`}
                checked={intent === candidate}
                onChange={() => setIntent(candidate)}
              />
              {candidate === 'suggest' ? 'Suggest values' : 'Explain a validation issue'}
            </label>
          ))}
        </fieldset>
        {intent === 'suggest' && (
          <fieldset className="space-y-3">
            <legend className="text-sm font-medium">Fields</legend>
            {capability.fields.length === 0 && (
              <p role="status">No suggestion fields were declared.</p>
            )}
            {capability.fields.map((field) => (
              <FieldPreference
                key={field.token}
                field={field}
                baseId={baseId}
                selected={selected.has(field.token)}
                numberValue={numbers[field.token] ?? ''}
                textValue={texts[field.token] ?? ''}
                flagValue={flags[field.token] ?? false}
                onToggle={() => toggleField(field.token)}
                onNumber={(value) =>
                  setNumbers((current) => ({ ...current, [field.token]: value }))
                }
                onText={(value) =>
                  setTexts((current) => ({ ...current, [field.token]: value }))
                }
                onFlag={(value) =>
                  setFlags((current) => ({ ...current, [field.token]: value }))
                }
              />
            ))}
          </fieldset>
        )}
        {intent === 'explain_issues' && (
          <fieldset className="space-y-2">
            <legend className="text-sm font-medium">Validation issue</legend>
            {capability.issues.length === 0 && (
              <p role="status">No explainable validation issues were declared.</p>
            )}
            {capability.issues.map((issue) => (
              <label key={issue.code} className="flex min-h-11 items-center gap-2 text-sm">
                <input
                  type="radio"
                  name={`${baseId}-issue`}
                  checked={issueCode === issue.code}
                  onChange={() => setIssueCode(issue.code)}
                />
                {issue.label}
              </label>
            ))}
          </fieldset>
        )}
        {requestStatus === 'unconfirmed' && (
          <p role="status">
            The suggestion response was not confirmed. This page will not send
            it again on its own.
          </p>
        )}
        {requestStatus === 'failed' && (
          <p role="alert">
            The suggestion request did not succeed. Nothing was adopted or
            validated.
          </p>
        )}
        {formError && <p role="alert">{formError}</p>}
        <Button
          type="submit"
          variant="primary"
          disabled={blockReason !== null || intent === null}
          aria-describedby={
            blockReason !== null ? `${baseId}-request-reason` : undefined
          }
        >
          Request suggestion
        </Button>
        {blockReason !== null && (
          <p id={`${baseId}-request-reason`} className="text-xs text-[var(--color-text-muted)]">
            {blockReason}
          </p>
        )}
      </form>
    </Panel>
  );
}

function FieldPreference({
  field,
  baseId,
  selected,
  numberValue,
  textValue,
  flagValue,
  onToggle,
  onNumber,
  onText,
  onFlag,
}: {
  field: CapabilityField;
  baseId: string;
  selected: boolean;
  numberValue: string;
  textValue: string;
  flagValue: boolean;
  onToggle: () => void;
  onNumber: (value: string) => void;
  onText: (value: string) => void;
  onFlag: (value: boolean) => void;
}) {
  const inputId = `${baseId}-${field.token}-preference`;
  return (
    <div className="space-y-2 rounded-[4px] border border-[var(--color-border)] p-3">
      <label className="flex min-h-11 items-center gap-2 text-sm">
        <input type="checkbox" checked={selected} onChange={onToggle} />
        {field.label}
      </label>
      {selected && field.type === 'boolean' && (
        <label className="flex min-h-11 items-center gap-2 text-sm" htmlFor={inputId}>
          <input
            id={inputId}
            type="checkbox"
            checked={flagValue}
            onChange={(event) => onFlag(event.target.checked)}
          />
          Preferred value
        </label>
      )}
      {selected && (field.type === 'integer' || field.type === 'number') && (
        <label className="block text-sm" htmlFor={inputId}>
          Preferred value
          <input
            id={inputId}
            className="mt-1 block w-full rounded-[4px] border border-[var(--color-border)] px-2 py-2"
            type="text"
            inputMode={field.type === 'integer' ? 'numeric' : 'decimal'}
            value={numberValue}
            onChange={(event) => onNumber(event.target.value)}
          />
          <ConstraintHint field={field} />
        </label>
      )}
      {selected && field.type === 'string' && (
        <label className="block text-sm" htmlFor={inputId}>
          Preferred value
          <input
            id={inputId}
            className="mt-1 block w-full rounded-[4px] border border-[var(--color-border)] px-2 py-2"
            type="text"
            value={textValue}
            onChange={(event) => onText(event.target.value)}
          />
        </label>
      )}
    </div>
  );
}

function ConstraintHint({ field }: { field: CapabilityField }) {
  const parts: string[] = [];
  if (typeof field.minimum === 'number') parts.push(`Minimum ${field.minimum}`);
  if (typeof field.maximum === 'number') parts.push(`Maximum ${field.maximum}`);
  if (parts.length === 0) return null;
  return (
    <span className="mt-1 block text-xs text-[var(--color-text-muted)]">
      {parts.join('. ')}.
    </span>
  );
}

function requestBlockReason(input: {
  writesEnabled: boolean;
  dirty: boolean;
  connectionAvailable: boolean;
  draftExpired: boolean;
  requestStatus: SuggestionRequestStatus;
  onRequest: boolean;
}): string | null {
  if (input.draftExpired) {
    return 'This draft has expired. You can review adoption receipts, not request a new suggestion.';
  }
  if (input.requestStatus === 'unconfirmed') {
    return 'The suggestion response was not confirmed. This page will not send it again on its own.';
  }
  const blocked = writeBlockReason(input.writesEnabled, input.dirty);
  if (blocked) return blocked;
  if (!input.connectionAvailable) {
    return 'Connect the assistant before requesting a suggestion.';
  }
  if (!input.onRequest) return 'The suggestion request is not connected yet.';
  return null;
}

function buildQuestion(input: {
  intent: SuggestionIntent;
  fields: readonly CapabilityField[];
  selected: ReadonlySet<string>;
  numbers: Record<string, string>;
  texts: Record<string, string>;
  flags: Record<string, boolean>;
  issueCode: string | null;
}): { ok: true; value: StructuredQuestion } | { ok: false; message: string } {
  if (input.intent === 'explain_issues') {
    if (!input.issueCode) {
      return { ok: false, message: 'Choose a declared validation issue.' };
    }
    return {
      ok: true,
      value: {
        intent: input.intent,
        fieldTokens: [],
        preferences: [],
        issueCode: input.issueCode,
      },
    };
  }
  const chosen = input.fields.filter((field) => input.selected.has(field.token));
  if (chosen.length === 0) {
    return { ok: false, message: 'Select a declared suggestion field.' };
  }
  const preferences: { fieldToken: string; value: string | number | boolean }[] = [];
  for (const field of chosen) {
    const parsed = readPreference(field, input);
    if (parsed === 'invalid') {
      return {
        ok: false,
        message: `Enter a valid preferred value for ${field.label}.`,
      };
    }
    if (parsed !== 'omit') {
      preferences.push({ fieldToken: field.token, value: parsed });
    }
  }
  return {
    ok: true,
    value: {
      intent: 'suggest',
      fieldTokens: chosen.map((field) => field.token),
      preferences,
      issueCode: null,
    },
  };
}

function readPreference(
  field: CapabilityField,
  input: {
    numbers: Record<string, string>;
    texts: Record<string, string>;
    flags: Record<string, boolean>;
  },
): string | number | boolean | 'omit' | 'invalid' {
  if (field.type === 'boolean') return input.flags[field.token] ?? false;
  if (field.type === 'string') {
    const text = input.texts[field.token] ?? '';
    return text === '' ? 'omit' : text;
  }
  const raw = (input.numbers[field.token] ?? '').trim();
  if (raw === '') return 'omit';
  const numeric = Number(raw);
  if (!Number.isFinite(numeric)) return 'invalid';
  if (field.type === 'integer' && !Number.isInteger(numeric)) return 'invalid';
  if (typeof field.minimum === 'number' && numeric < field.minimum) return 'invalid';
  if (typeof field.maximum === 'number' && numeric > field.maximum) return 'invalid';
  return numeric;
}
