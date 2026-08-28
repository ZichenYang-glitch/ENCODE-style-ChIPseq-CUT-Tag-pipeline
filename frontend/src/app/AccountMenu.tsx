import { useEffect, useRef, useState } from 'react';
import { ChevronDown, LogOut, UserRound } from 'lucide-react';
import type {
  PrincipalResponse,
  TerminalEmailPreferenceResponse,
} from '../api/generated/models';
import { Button } from '../components/Button';

interface AccountMenuProps {
  principal: PrincipalResponse;
  terminalEmailPreference?: TerminalEmailPreferenceResponse;
  preferencePending: boolean;
  onPreferenceChange: (enabled: boolean) => void;
  onSignOut: () => void;
}

export function AccountMenu({
  principal,
  terminalEmailPreference,
  preferencePending,
  onPreferenceChange,
  onSignOut,
}: AccountMenuProps) {
  const [open, setOpen] = useState(false);
  const rootRef = useRef<HTMLDivElement>(null);
  const triggerRef = useRef<HTMLButtonElement>(null);
  const panelRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!open) return;
    panelRef.current
      ?.querySelector<HTMLElement>('input:not(:disabled), button:not(:disabled)')
      ?.focus();
    function handlePointerDown(event: MouseEvent) {
      if (!rootRef.current?.contains(event.target as Node)) setOpen(false);
    }
    function handleKeyDown(event: KeyboardEvent) {
      if (event.key !== 'Escape') return;
      setOpen(false);
      triggerRef.current?.focus();
    }
    document.addEventListener('mousedown', handlePointerDown);
    document.addEventListener('keydown', handleKeyDown);
    return () => {
      document.removeEventListener('mousedown', handlePointerDown);
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [open]);

  return (
    <div ref={rootRef} className="relative min-w-0">
      <Button
        ref={triggerRef}
        type="button"
        variant="quiet"
        className="max-w-[11rem] gap-1.5 px-2"
        aria-expanded={open}
        aria-controls="account-menu-panel"
        onClick={() => setOpen((current) => !current)}
      >
        <UserRound aria-hidden="true" size={16} />
        <span className="truncate" data-testid="current-username">
          {principal.username}
        </span>
        <ChevronDown aria-hidden="true" size={14} />
      </Button>
      {open && (
        <div
          ref={panelRef}
          id="account-menu-panel"
          aria-label="Account"
          className="absolute right-0 z-20 mt-1 w-[min(20rem,calc(100vw-1.5rem))] rounded-[6px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] p-3 shadow-md"
        >
          <div className="border-b border-[var(--color-border)] pb-3">
            <p className="truncate text-sm font-semibold text-[var(--color-text)]">
              {principal.username}
            </p>
            <p className="mt-0.5 text-xs capitalize text-[var(--color-text-muted)]">
              {principal.role}
            </p>
          </div>
          {principal.role === 'member' && terminalEmailPreference && (
            <div className="border-b border-[var(--color-border)] py-3">
              <label className="flex min-h-11 items-start gap-2 text-sm text-[var(--color-text)]">
                <input
                  className="mt-1 h-4 w-4"
                  aria-label="Email me when my runs finish"
                  checked={terminalEmailPreference.terminal_email_enabled}
                  disabled={
                    !terminalEmailPreference.address_configured ||
                    preferencePending
                  }
                  type="checkbox"
                  onChange={(event) =>
                    onPreferenceChange(event.currentTarget.checked)
                  }
                />
                <span>
                  Terminal emails
                  <span className="mt-0.5 block text-xs text-[var(--color-text-muted)]">
                    {terminalEmailPreference.address_configured
                      ? 'Notify me when my runs reach a terminal state.'
                      : 'Address not configured. Ask the local operator to configure it.'}
                  </span>
                </span>
              </label>
            </div>
          )}
          {principal.role === 'administrator' && (
            <p className="border-b border-[var(--color-border)] py-3 text-xs text-[var(--color-text-muted)]">
              Terminal notifications are operator-managed for administrators.
            </p>
          )}
          <Button
            className="mt-2 w-full justify-start gap-2"
            variant="quiet"
            onClick={onSignOut}
          >
            <LogOut aria-hidden="true" size={16} />
            Sign out
          </Button>
        </div>
      )}
    </div>
  );
}
