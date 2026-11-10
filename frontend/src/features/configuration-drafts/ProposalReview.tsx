import { useId } from 'react';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { adoptBlockReason, rejectBlockReason } from './actionGate';
import type { ProposalView } from './model';
import { formatScalarValue } from './scalarDisplay';
import { proposalStatusCopy } from './statusCopy';

interface ProposalReviewProps {
  proposal: ProposalView | null;
  writesEnabled: boolean;
  dirty: boolean;
  onAdopt?: () => void;
  onReject?: () => void;
  onRefresh?: () => void;
}

export function ProposalReview({
  proposal,
  writesEnabled,
  dirty,
  onAdopt,
  onReject,
  onRefresh,
}: ProposalReviewProps) {
  if (!proposal) {
    return (
      <Panel title="Suggested changes">
        <p role="status">No suggestion is waiting for review.</p>
      </Panel>
    );
  }
  const adoptReason = adoptBlockReason({ writesEnabled, dirty, proposal });
  const rejectReason = rejectBlockReason({ writesEnabled, dirty, proposal });
  const refreshDisabled = proposal.state !== 'unconfirmed' || !onRefresh;
  const adoptDescriptionId = useId();
  const rejectDescriptionId = useId();
  const sharedReason = adoptReason !== null && adoptReason === rejectReason;
  return (
    <Panel title="Suggested changes">
      <div className="space-y-3">
        <p role="status" aria-live="polite">
          {proposalStatusCopy(proposal.state)}
        </p>
        {proposal.changes.length === 0 ? (
          <p>This suggestion does not list any field changes.</p>
        ) : (
          <div className="overflow-x-auto">
            <table className="w-full min-w-[18rem] border-collapse text-left text-sm">
              <caption className="sr-only">Suggested changes</caption>
              <thead>
                <tr className="border-b border-[var(--color-border)] text-xs text-[var(--color-text-muted)]">
                  <th scope="col" className="px-2 py-2 font-medium">Field</th>
                  <th scope="col" className="px-2 py-2 font-medium">Before</th>
                  <th scope="col" className="px-2 py-2 font-medium">After</th>
                </tr>
              </thead>
              <tbody>
                {proposal.changes.map((change, index) => (
                  <tr
                    key={`${change.label}-${index}`}
                    className="border-b border-[var(--color-border)] align-top"
                  >
                    <th scope="row" className="px-2 py-2 font-medium">
                      {change.label || 'Unlabeled field'}
                      {!change.recognized && (
                        <span className="mt-1 block text-xs font-normal text-[var(--color-text-muted)]">
                          This field is outside the current suggestion list, so
                          it cannot be edited or requested.
                        </span>
                      )}
                    </th>
                    <td className="px-2 py-2">{formatScalarValue(change.before)}</td>
                    <td className="px-2 py-2">{formatScalarValue(change.after)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
        <div className="flex flex-col gap-2 sm:flex-row">
          <Button
            type="button"
            variant="primary"
            disabled={adoptReason !== null || !onAdopt}
            aria-describedby={adoptReason !== null ? adoptDescriptionId : undefined}
            onClick={() => {
              if (adoptReason !== null) return;
              onAdopt?.();
            }}
          >
            Adopt suggestion
          </Button>
          <Button
            type="button"
            variant="danger"
            disabled={rejectReason !== null || !onReject}
            aria-describedby={
              rejectReason === null
                ? undefined
                : sharedReason
                  ? adoptDescriptionId
                  : rejectDescriptionId
            }
            onClick={() => {
              if (rejectReason !== null) return;
              onReject?.();
            }}
          >
            Reject suggestion
          </Button>
          <Button
            type="button"
            variant="secondary"
            disabled={refreshDisabled}
            onClick={() => onRefresh?.()}
          >
            Refresh server status
          </Button>
        </div>
        {adoptReason !== null && (
          <p id={adoptDescriptionId} className="text-xs text-[var(--color-text-muted)]">
            {adoptReason}
          </p>
        )}
        {rejectReason !== null && !sharedReason && (
          <p id={rejectDescriptionId} className="text-xs text-[var(--color-text-muted)]">
            {rejectReason}
          </p>
        )}
      </div>
    </Panel>
  );
}
