import type { ReactNode } from 'react';

interface PanelProps {
  title: string;
  children: ReactNode;
  className?: string;
}

export function Panel({ title, children, className = '' }: PanelProps) {
  return (
    <section
      className={`rounded border border-[var(--color-border)] bg-[var(--color-surface)] ${className}`}
    >
      <header className="border-b border-[var(--color-border)] px-3 py-2 text-sm font-medium text-[var(--color-text)]">
        {title}
      </header>
      <div className="p-3">{children}</div>
    </section>
  );
}
