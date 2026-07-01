import type { ButtonHTMLAttributes, ReactNode } from 'react';

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  children: ReactNode;
  variant?: 'primary' | 'secondary';
}

export function Button({
  children,
  variant = 'secondary',
  className = '',
  ...props
}: ButtonProps) {
  const base =
    'inline-flex items-center justify-center rounded px-3 py-1.5 text-sm font-medium transition-colors';
  const variants =
    variant === 'primary'
      ? 'bg-[var(--color-accent)] text-white hover:opacity-90'
      : 'border border-[var(--color-border)] bg-[var(--color-bg)] text-[var(--color-text)] hover:bg-[var(--color-surface)]';
  return (
    <button className={`${base} ${variants} ${className}`} {...props}>
      {children}
    </button>
  );
}
