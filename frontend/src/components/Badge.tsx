import type { HTMLAttributes, ReactNode } from 'react';
import type { Severity } from '../api/types';

export type BadgeTone = 'success' | 'warning' | 'error' | 'info' | 'neutral';

interface BadgeProps extends HTMLAttributes<HTMLSpanElement> {
  severity?: Severity;
  tone?: BadgeTone;
  children?: ReactNode;
}

const severityLabels: Record<Severity, string> = {
  error: 'Error',
  warning: 'Warning',
  info: 'Info',
};

const toneClasses: Record<BadgeTone, string> = {
  success:
    'border-[var(--color-success-border)] bg-[var(--color-success-bg)] text-[var(--color-success)]',
  error: 'border-[var(--color-error-border)] bg-[var(--color-error-bg)] text-[var(--color-error)]',
  warning:
    'border-[var(--color-warning-border)] bg-[var(--color-warning-bg)] text-[var(--color-warning)]',
  info: 'border-[var(--color-info-border)] bg-[var(--color-info-bg)] text-[var(--color-info)]',
  neutral:
    'border-[var(--color-neutral-border)] bg-[var(--color-neutral-bg)] text-[var(--color-neutral)]',
};

export function Badge({
  severity,
  tone,
  children,
  className = '',
  ...props
}: BadgeProps) {
  const resolvedTone = tone ?? severity ?? 'neutral';
  return (
    <span
      className={`inline-flex items-center rounded-[4px] border px-1.5 py-0.5 text-xs font-semibold ${toneClasses[resolvedTone]} ${className}`}
      {...props}
    >
      {children ?? (severity ? severityLabels[severity] : 'Neutral')}
    </span>
  );
}
