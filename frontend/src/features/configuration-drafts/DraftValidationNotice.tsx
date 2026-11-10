import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { validateBlockReason } from './actionGate';
import type { DraftValidationView } from './model';
import { validationLead } from './statusCopy';

interface DraftValidationNoticeProps {
  validation: DraftValidationView;
  writesEnabled: boolean;
  unavailableMessage?: string;
  dirty: boolean;
  draftExpired: boolean;
  onValidate?: () => void;
}

export function DraftValidationNotice({
  validation,
  writesEnabled,
  unavailableMessage,
  dirty,
  draftExpired,
  onValidate,
}: DraftValidationNoticeProps) {
  const reason = validateBlockReason({
    writesEnabled,
    dirty,
    draftExpired,
    unavailable: unavailableMessage,
  });
  return (
    <Panel title="Draft validation">
      <div className="space-y-3">
        <p role="status" aria-live="polite">
          {validationLead(validation)}
        </p>
        {validation.summary && <p className="text-sm">{validation.summary}</p>}
        {validation.status === 'snapshot' && (
          <dl className="grid grid-cols-1 gap-2 text-sm sm:grid-cols-2">
            <div>
              <dt className="text-xs text-[var(--color-text-muted)]">Draft revision</dt>
              <dd>{validation.revisionLabel ?? 'Not reported'}</dd>
            </div>
            <div>
              <dt className="text-xs text-[var(--color-text-muted)]">Valid until</dt>
              <dd>{validation.expiresLabel ?? 'Not reported'}</dd>
            </div>
          </dl>
        )}
        <Button
          type="button"
          variant="secondary"
          disabled={reason !== null || !onValidate}
          aria-describedby={reason !== null ? 'draft-validation-reason' : undefined}
          onClick={() => {
            if (reason !== null) return;
            onValidate?.();
          }}
        >
          Validate draft
        </Button>
        {reason !== null && (
          <p id="draft-validation-reason" className="text-xs text-[var(--color-text-muted)]">
            {reason}
          </p>
        )}
      </div>
    </Panel>
  );
}
