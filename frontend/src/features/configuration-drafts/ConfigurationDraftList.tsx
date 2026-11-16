import { useRef } from 'react';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { boundedPageLimit, type ConfigurationDraftSummary, type DraftListModel } from './model';

interface ConfigurationDraftListProps {
  model: DraftListModel;
  selectedSid: string | null;
  onSelect?: (sid: string) => void;
  onPrevious?: () => void;
  onNext?: () => void;
  onRetry?: () => void;
}

export function ConfigurationDraftList({
  model,
  selectedSid,
  onSelect,
  onPrevious,
  onNext,
  onRetry,
}: ConfigurationDraftListProps) {
  const limit = boundedPageLimit(model.pageLimit);
  const visible = model.drafts.slice(0, limit);
  const listRef = useRef<HTMLUListElement>(null);

  function move(current: HTMLButtonElement, direction: 1 | -1) {
    const buttons = Array.from(
      listRef.current?.querySelectorAll<HTMLButtonElement>(
        'button[data-draft-option]',
      ) ?? [],
    );
    const index = buttons.indexOf(current);
    buttons[index + direction]?.focus();
  }

  return (
    <Panel title="My configuration drafts">
      <div className="space-y-3">
        {model.status === 'loading' && (
          <p role="status">Loading your configuration drafts.</p>
        )}
        {model.status === 'empty' && (
          <p role="status">No configuration drafts for this workflow.</p>
        )}
        {model.status === 'unavailable' && (
          <p role="status">
            The draft list is not connected yet. Nothing was changed.
          </p>
        )}
        {model.status === 'error' && (
          <div className="space-y-2">
            <p role="alert">
              {model.errorMessage ??
                'The draft list could not be loaded. Nothing was changed.'}
            </p>
            <Button
              type="button"
              variant="secondary"
              disabled={!onRetry}
              onClick={() => onRetry?.()}
            >
              Retry draft list
            </Button>
          </div>
        )}
        {model.status === 'ready' && visible.length === 0 && (
          <p role="status">No configuration drafts for this workflow.</p>
        )}
        {model.status === 'ready' && visible.length > 0 && (
          <ul
            ref={listRef}
            aria-label="My configuration drafts"
            className="space-y-2"
            data-draft-list
          >
            {visible.map((draft) => (
              <li key={draft.sid}>
                <DraftOption
                  draft={draft}
                  selected={draft.sid === selectedSid}
                  onSelect={onSelect}
                  onMove={move}
                />
              </li>
            ))}
          </ul>
        )}
        {model.status === 'ready' && (
          <div className="flex flex-col gap-2 sm:flex-row sm:items-center sm:justify-between">
            <p className="text-xs text-[var(--color-text-muted)]">
              Up to {limit} drafts per page.
            </p>
            <div className="flex gap-2">
              <Button
                type="button"
                variant="secondary"
                disabled={!model.hasPrevious || !onPrevious}
                onClick={() => onPrevious?.()}
              >
                Previous page
              </Button>
              <Button
                type="button"
                variant="secondary"
                disabled={!model.hasNext || !onNext}
                onClick={() => onNext?.()}
              >
                Next page
              </Button>
            </div>
          </div>
        )}
      </div>
    </Panel>
  );
}

function DraftOption({
  draft,
  selected,
  onSelect,
  onMove,
}: {
  draft: ConfigurationDraftSummary;
  selected: boolean;
  onSelect?: (sid: string) => void;
  onMove: (current: HTMLButtonElement, direction: 1 | -1) => void;
}) {
  const expired = draft.state === 'expired';
  return (
    <button
      type="button"
      data-draft-option
      aria-pressed={selected}
      className="min-h-11 w-full rounded-[4px] border border-[var(--color-border)] px-3 py-2 text-left focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-[var(--color-focus)]"
      onClick={() => onSelect?.(draft.sid)}
      onKeyDown={(event) => {
        if (event.key !== 'ArrowDown' && event.key !== 'ArrowUp') return;
        event.preventDefault();
        onMove(event.currentTarget, event.key === 'ArrowDown' ? 1 : -1);
      }}
    >
      <span className="block text-sm font-medium">{draft.projectName}</span>
      <span className="mt-1 block text-xs text-[var(--color-text-muted)]">
        {draft.workflowName}
      </span>
      <span className="mt-1 block text-xs text-[var(--color-text-muted)]">
        Updated {draft.updatedLabel}. Expires {draft.expiresLabel}.
      </span>
      {expired && (
        <span className="mt-1 block text-xs text-[var(--color-error)]">
          Expired. The configuration stays hidden. Adoption receipts can still
          be reviewed.
        </span>
      )}
    </button>
  );
}
