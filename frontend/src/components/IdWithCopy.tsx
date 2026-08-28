import { useEffect, useState } from 'react';
import { Check, Copy, ExternalLink } from 'lucide-react';
import { Link } from 'react-router-dom';
import { Button } from './Button';

interface IdWithCopyProps {
  value: string;
  label?: string;
  to?: string;
  shortLength?: number;
  className?: string;
}

export function IdWithCopy({
  value,
  label = 'ID',
  to,
  shortLength = 8,
  className = '',
}: IdWithCopyProps) {
  const [copied, setCopied] = useState(false);

  useEffect(() => {
    if (!copied) return;
    const timer = window.setTimeout(() => setCopied(false), 1500);
    return () => window.clearTimeout(timer);
  }, [copied]);

  const visible =
    value.length > shortLength ? `${value.slice(0, shortLength)}…` : value;
  const identity = (
    <code
      className="block min-w-0 truncate font-mono text-xs tabular-nums"
      data-technical-value
      title={value}
    >
      {visible}
    </code>
  );

  async function handleCopy() {
    try {
      await navigator.clipboard.writeText(value);
      setCopied(true);
    } catch {
      setCopied(false);
    }
  }

  return (
    <span className={`inline-flex min-w-0 max-w-full items-center gap-1 ${className}`}>
      {to ? (
        <Link
          className="inline-flex min-w-0 items-center gap-1 text-[var(--color-link)] hover:underline"
          to={to}
          aria-label={`Open ${label} ${value}`}
          title={value}
        >
          {identity}
          <ExternalLink aria-hidden="true" className="shrink-0" size={13} />
        </Link>
      ) : (
        identity
      )}
      <Button
        type="button"
        variant="quiet"
        size="icon"
        className="shrink-0"
        onClick={() => void handleCopy()}
        aria-label={`Copy full ${label} ${value}`}
        title={copied ? `${label} copied` : `Copy full ${label}`}
      >
        {copied ? (
          <Check aria-hidden="true" size={15} />
        ) : (
          <Copy aria-hidden="true" size={15} />
        )}
      </Button>
      <span className="sr-only" role="status" aria-live="polite">
        {copied ? `${label} copied.` : ''}
      </span>
    </span>
  );
}
