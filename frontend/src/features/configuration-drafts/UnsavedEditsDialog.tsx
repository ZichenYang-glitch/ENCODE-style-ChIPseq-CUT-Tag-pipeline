import { useEffect, useRef, type KeyboardEvent } from 'react';
import { Button } from '../../components/Button';

interface UnsavedEditsDialogProps {
  open: boolean;
  writesEnabled: boolean;
  saveAttempted: boolean;
  onSave: () => void;
  onDiscard: () => void;
  onCancel: () => void;
}

export function UnsavedEditsDialog({
  open,
  writesEnabled,
  saveAttempted,
  onSave,
  onDiscard,
  onCancel,
}: UnsavedEditsDialogProps) {
  const dialogRef = useRef<HTMLDivElement>(null);
  const cancelRef = useRef<HTMLButtonElement>(null);

  useEffect(() => {
    if (!open) return;
    cancelRef.current?.focus();
  }, [open]);

  if (!open) return null;

  function trapTab(event: KeyboardEvent<HTMLDivElement>) {
    if (event.key === 'Escape') {
      event.preventDefault();
      onCancel();
      return;
    }
    if (event.key !== 'Tab' || !dialogRef.current) return;
    const items = Array.from(
      dialogRef.current.querySelectorAll<HTMLButtonElement>(
        'button:not(:disabled)',
      ),
    );
    if (items.length === 0) return;
    const first = items[0];
    const last = items[items.length - 1];
    if (event.shiftKey && document.activeElement === first) {
      event.preventDefault();
      last.focus();
    } else if (!event.shiftKey && document.activeElement === last) {
      event.preventDefault();
      first.focus();
    }
  }

  return (
    <div className="fixed inset-0 z-40 flex items-end justify-center bg-black/40 p-3 sm:items-center">
      <div
        ref={dialogRef}
        role="dialog"
        aria-modal="true"
        aria-labelledby="unsaved-edits-title"
        className="w-full max-w-lg rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4 shadow-sm"
        onKeyDown={trapTab}
      >
        <h2 id="unsaved-edits-title" className="text-base font-semibold">
          Unsaved edits
        </h2>
        <p className="mt-2 text-sm text-[var(--color-text-muted)]">
          These edits stay on this page until you save them or discard them.
          Cancel keeps the edits and stays here.
        </p>
        {saveAttempted && (
          <p className="mt-2 text-sm" role="status">
            The save response has not been confirmed. These edits are still
            unsaved.
          </p>
        )}
        {!writesEnabled && (
          <p className="mt-2 text-sm" id="unsaved-save-reason">
            Saving waits for the reviewed configuration-draft interface.
          </p>
        )}
        <div className="mt-4 flex flex-col gap-2 sm:flex-row sm:justify-end">
          <Button
            type="button"
            variant="primary"
            disabled={!writesEnabled}
            aria-describedby={writesEnabled ? undefined : 'unsaved-save-reason'}
            onClick={onSave}
          >
            Save edits
          </Button>
          <Button type="button" variant="danger" onClick={onDiscard}>
            Discard edits
          </Button>
          <Button ref={cancelRef} type="button" variant="secondary" onClick={onCancel}>
            Cancel
          </Button>
        </div>
      </div>
    </div>
  );
}
