import { useEffect, useState } from 'react';
import { Check, Copy } from 'lucide-react';
import { Button } from './Button';

interface SampleIdentityProps {
  sampleId: string | null | undefined;
  experimentId?: string | null;
}

/** Keep whitespace-bearing identifiers distinguishable without changing their value. */
export function SampleIdentity({ sampleId, experimentId }: SampleIdentityProps) {
  const [copyStatus, setCopyStatus] = useState<'idle' | 'copied' | 'failed'>('idle');
  useEffect(() => setCopyStatus('idle'), [sampleId]);
  const sampleLabel = sampleId?.includes(' ') ? JSON.stringify(sampleId) : sampleId;
  const visible = [sampleLabel, experimentId].filter(Boolean).join(' · ') || '—';

  async function copySample() {
    if (!sampleId) return;
    try {
      await navigator.clipboard.writeText(sampleId);
      setCopyStatus('copied');
    } catch {
      setCopyStatus('failed');
    }
  }

  return (
    <span
      className="inline-flex min-w-0 max-w-full items-start gap-1"
      data-sample-identity={sampleId ?? undefined}
    >
      <span
        className="min-w-0 font-mono whitespace-pre-wrap break-words [overflow-wrap:anywhere]"
        data-sample-label
      >
        {visible}
      </span>
      {sampleId && (
        <Button
          type="button"
          variant="quiet"
          size="icon"
          className="shrink-0"
          onClick={() => void copySample()}
          aria-label="Copy sample ID"
          title={copyStatus === 'copied' ? 'Sample ID copied' : 'Copy sample ID exactly'}
        >
          {copyStatus === 'copied' ? (
            <Check size={15} aria-hidden="true" />
          ) : (
            <Copy size={15} aria-hidden="true" />
          )}
        </Button>
      )}
      <span
        className="sr-only"
        role={copyStatus === 'idle' ? undefined : 'status'}
        aria-live="polite"
      >
        {copyStatus === 'copied'
          ? 'Sample ID copied.'
          : copyStatus === 'failed'
            ? 'Could not copy sample ID.'
            : ''}
      </span>
    </span>
  );
}
