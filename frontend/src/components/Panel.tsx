import type { ReactNode } from 'react';

interface PanelProps {
  title: string;
  children: ReactNode;
  className?: string;
}

export function Panel({ title, children, className = '' }: PanelProps) {
  return (
    <section
      className={`rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] shadow-sm ${className}`}
    >
      <header className="border-b border-[var(--color-border)] bg-[var(--color-bg)] px-3 py-2">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]">
          {title}
        </h2>
      </header>
      <div className="p-3">{children}</div>
    </section>
  );
}
