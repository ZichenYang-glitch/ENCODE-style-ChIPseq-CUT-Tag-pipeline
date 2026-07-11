import { Outlet } from 'react-router-dom';

export function AppShell() {
  return (
    <div className="flex min-h-screen flex-col bg-[var(--color-bg)] text-[var(--color-text)]">
      <header className="border-b border-[var(--color-border)] bg-[var(--color-surface)] px-4 py-3">
        <h1 className="text-base font-semibold tracking-wide text-[var(--color-accent)]">
          Workflow Platform
        </h1>
      </header>
      <main className="mx-auto flex w-full max-w-screen-2xl flex-1 flex-col gap-3 p-3 lg:flex-row">
        <Outlet />
      </main>
    </div>
  );
}
