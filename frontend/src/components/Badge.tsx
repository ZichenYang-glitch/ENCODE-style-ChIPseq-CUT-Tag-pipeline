import type { Severity } from '../api/types';

interface BadgeProps {
  severity: Severity;
}

const severityLabels: Record<Severity, string> = {
  error: 'Error',
  warning: 'Warning',
  info: 'Info',
};

const severityClasses: Record<Severity, string> = {
  error: 'bg-[var(--color-error-bg)] text-[var(--color-error)] border-[var(--color-error)]',
  warning:
    'bg-[var(--color-warning-bg)] text-[var(--color-warning)] border-[var(--color-warning)]',
  info: 'bg-[var(--color-info-bg)] text-[var(--color-info)] border-[var(--color-info)]',
};

export function Badge({ severity }: BadgeProps) {
  return (
    <span
      className={`inline-flex items-center rounded border px-1.5 py-0.5 text-xs font-semibold ${severityClasses[severity]}`}
    >
      {severityLabels[severity]}
    </span>
  );
}
