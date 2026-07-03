import type { RunEventResponse } from '../../api/runTypes';

interface RunEventFeedProps {
  events: RunEventResponse[];
}

export function RunEventFeed({ events }: RunEventFeedProps) {
  if (events.length === 0) {
    return (
      <p
        className="text-sm text-[var(--color-text-muted)]"
        data-testid="run-event-feed-empty"
      >
        No events recorded yet.
      </p>
    );
  }

  return (
    <ul className="space-y-2" data-testid="run-event-feed">
      {events.map((event) => (
        <li
          key={event.event_id}
          className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-2 text-sm"
        >
          <div className="flex flex-wrap items-center gap-2">
            <span className="font-medium text-[var(--color-text)]">
              {event.event_type}
            </span>
            {event.status && (
              <span className="text-xs text-[var(--color-text-muted)]">
                {event.status}
              </span>
            )}
            <span className="ml-auto text-xs text-[var(--color-text-muted)]">
              {new Date(event.timestamp).toLocaleString()}
            </span>
          </div>
          <p className="mt-1 text-[var(--color-text)]">{event.message}</p>
          {Object.keys(event.context).length > 0 && (
            <details className="mt-1">
              <summary className="cursor-pointer text-xs text-[var(--color-accent)] hover:underline">
                Context
              </summary>
              <pre className="mt-1 overflow-auto rounded bg-[var(--color-surface)] p-1.5 text-xs text-[var(--color-text)]">
                {JSON.stringify(event.context, null, 2)}
              </pre>
            </details>
          )}
        </li>
      ))}
    </ul>
  );
}
