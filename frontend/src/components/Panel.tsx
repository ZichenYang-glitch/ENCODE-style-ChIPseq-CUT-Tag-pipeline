import type { ReactNode } from 'react';

interface PanelProps {
  title: string;
  children: ReactNode;
  className?: string;
}

export function Panel({ title, children, className = '' }: PanelProps) {
  return (
    <section
      className={`rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] ${className}`}
    >
      <header className="border-b border-[var(--color-border)] bg-[var(--color-surface-subtle)] px-4 py-2">
        <h2 className="text-xs font-medium text-[var(--color-text-muted)]">
          {title}
        </h2>
      </header>
      <div className="p-4">{children}</div>
    </section>
  );
}
