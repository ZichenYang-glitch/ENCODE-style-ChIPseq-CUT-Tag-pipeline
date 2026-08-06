import { RefreshCw } from 'lucide-react';
import type { ReferenceProfileSummary } from '../../api/runTypes';
import { Button } from '../../components/Button';

interface ReferenceProfileSelectorProps {
  profiles: readonly ReferenceProfileSummary[] | null;
  selectedRevisionId: string | null;
  loading: boolean;
  refreshing: boolean;
  errorMessage: string | null;
  selectionMessage: string | null;
  onSelect: (revisionId: string | null) => void;
  onRefresh: () => void;
}

export function ReferenceProfileSelector({
  profiles,
  selectedRevisionId,
  loading,
  refreshing,
  errorMessage,
  selectionMessage,
  onSelect,
  onRefresh,
}: ReferenceProfileSelectorProps) {
  const selected =
    profiles?.find((profile) => profile.revision_id === selectedRevisionId) ??
    null;
  const empty = profiles !== null && profiles.length === 0;

  return (
    <div className="border-b border-[var(--color-border)] py-3">
      <div className="flex min-w-0 flex-col gap-2 sm:flex-row sm:items-end sm:justify-between">
        <label className="min-w-0 flex-1 text-xs font-medium text-[var(--color-text-muted)]">
          Reference profile
          <select
            className="mt-1 block w-full max-w-xl rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2 py-2 text-sm text-[var(--color-text)] outline-none focus:ring-2 focus:ring-[var(--color-accent)] disabled:cursor-not-allowed disabled:opacity-60"
            value={selectedRevisionId ?? ''}
            disabled={loading || profiles === null || empty}
            onChange={(event) => onSelect(event.target.value || null)}
            aria-describedby="reference-profile-status"
          >
            <option value="">Select an enabled reference</option>
            {(profiles ?? []).map((profile) => (
              <option key={profile.revision_id} value={profile.revision_id}>
                {profile.display_name} — {profile.assembly} — {profile.organism}
              </option>
            ))}
          </select>
        </label>
        <Button
          type="button"
          variant="secondary"
          className="gap-1.5 self-start sm:self-auto"
          disabled={refreshing}
          onClick={onRefresh}
          aria-label="Refresh references"
        >
          <RefreshCw aria-hidden="true" size={14} />
          {refreshing ? 'Refreshing…' : 'Refresh'}
        </Button>
      </div>
      <div
        id="reference-profile-status"
        className="mt-1 min-h-5 text-xs text-[var(--color-text-muted)]"
        role={errorMessage !== null ? 'alert' : 'status'}
      >
        {errorMessage ??
          selectionMessage ??
          (loading
            ? 'Loading enabled references…'
            : empty
              ? 'No enabled reference profile is available for this workflow.'
              : selected
                ? `${selected.organism} · ${selected.assembly} · revision ${selected.revision_number}`
                : 'Select the exact reference revision to validate this draft.')}
      </div>
    </div>
  );
}
