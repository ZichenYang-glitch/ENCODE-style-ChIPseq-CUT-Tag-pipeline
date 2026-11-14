import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { formatScalarValue } from './scalarDisplay';
import type { ProposalView } from './model';

export type ReceiptHistoryStatus = 'idle' | 'loading' | 'ready' | 'error' | 'not-found';

export interface AdoptedReceiptRecord {
  turnId: string;
  adoptedLabel: string;
  revisionLabel: string;
  changedLabel: string;
  changes: ProposalView['changes'] | null;
  differenceMessage: string | null;
}

interface AdoptedReceiptHistoryProps {
  status: ReceiptHistoryStatus;
  message: string | null;
  expired: boolean;
  receipts: readonly AdoptedReceiptRecord[];
  hasPrevious: boolean;
  hasNext: boolean;
  onPrevious?: () => void;
  onNext?: () => void;
  onRetry?: () => void;
}

export function AdoptedReceiptHistory({
  status,
  message,
  expired,
  receipts,
  hasPrevious,
  hasNext,
  onPrevious,
  onNext,
  onRetry,
}: AdoptedReceiptHistoryProps) {
  return (
    <Panel title="Adopted changes">
      <div className="space-y-3" aria-label="Adopted changes">
        <p>
          These records were already adopted. This list is not a suggestion you
          can adopt again.
        </p>
        {expired && status === 'ready' && (
          <p role="status">
            This draft has expired. The adopted records stay here until the
            draft is removed.
          </p>
        )}
        {status === 'loading' && <p role="status">Loading adopted changes.</p>}
        {status === 'not-found' && (
          <p role="status">{message ?? 'This draft is no longer available.'}</p>
        )}
        {status === 'error' && (
          <div className="space-y-2">
            <p role="alert">{message ?? 'The adopted changes could not be loaded.'}</p>
            <Button type="button" variant="secondary" onClick={() => onRetry?.()}>
              Try adopted changes again
            </Button>
          </div>
        )}
        {status === 'ready' && receipts.length === 0 && (
          <p role="status">No adopted changes are recorded for this draft.</p>
        )}
        {status === 'ready' &&
          receipts.map((receipt) => (
            <article key={receipt.turnId} className="space-y-2" aria-label="Adopted record">
              <p className="text-sm">{receipt.adoptedLabel}</p>
              <p className="text-sm text-[var(--color-text-muted)]">
                {receipt.revisionLabel} {receipt.changedLabel}
              </p>
              {receipt.differenceMessage && <p role="status">{receipt.differenceMessage}</p>}
              {receipt.changes && receipt.changes.length === 0 && (
                <p>This adopted record does not list any field changes.</p>
              )}
              {receipt.changes && receipt.changes.length > 0 && (
                <div className="overflow-x-auto">
                  <table className="w-full min-w-[18rem] border-collapse text-left text-sm">
                    <caption className="sr-only">Adopted field changes</caption>
                    <thead>
                      <tr className="border-b border-[var(--color-border)] text-xs text-[var(--color-text-muted)]">
                        <th scope="col" className="px-2 py-2 font-medium">Field</th>
                        <th scope="col" className="px-2 py-2 font-medium">Before</th>
                        <th scope="col" className="px-2 py-2 font-medium">After</th>
                      </tr>
                    </thead>
                    <tbody>
                      {receipt.changes.map((change, index) => (
                        <tr
                          key={`${change.label}-${index}`}
                          className="border-b border-[var(--color-border)] align-top"
                        >
                          <th scope="row" className="px-2 py-2 font-medium">
                            {change.label || 'Unlabeled field'}
                          </th>
                          <td className="px-2 py-2">{formatScalarValue(change.before)}</td>
                          <td className="px-2 py-2">{formatScalarValue(change.after)}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </article>
          ))}
        {status === 'ready' && (hasPrevious || hasNext) && (
          <div className="flex flex-col gap-2 sm:flex-row">
            <Button
              type="button"
              variant="secondary"
              disabled={!hasPrevious}
              onClick={() => onPrevious?.()}
            >
              Previous adopted changes
            </Button>
            <Button
              type="button"
              variant="secondary"
              disabled={!hasNext}
              onClick={() => onNext?.()}
            >
              Next adopted changes
            </Button>
          </div>
        )}
      </div>
    </Panel>
  );
}
