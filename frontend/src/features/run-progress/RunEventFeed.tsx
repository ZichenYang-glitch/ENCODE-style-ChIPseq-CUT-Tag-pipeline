import type { RunEventResponse } from '../../api/runTypes';

interface RunEventFeedProps {
  events: RunEventResponse[];
}

const publicContextFields = new Set([
  'previous_status',
  'new_status',
  'reason_code',
  'cancellation_reason',
]);

function publicContext(context: Record<string, unknown>) {
  return Object.fromEntries(
    Object.entries(context).filter(
      ([key, value]) => publicContextFields.has(key) && typeof value === 'string',
    ),
  );
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
      {events.map((event) => {
        const context = publicContext(event.context);
        return (
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
          {event.issue && (
            <div className="mt-1 text-xs text-[var(--color-error)]">
              <p>
                {event.issue.code}: {event.issue.message}
              </p>
              {event.issue.hint && <p className="mt-0.5">{event.issue.hint}</p>}
            </div>
          )}
          {Object.keys(context).length > 0 && (
            <details className="mt-1">
              <summary className="cursor-pointer text-xs text-[var(--color-accent)] hover:underline">
                Context
              </summary>
              <pre className="mt-1 overflow-auto rounded bg-[var(--color-surface)] p-1.5 text-xs text-[var(--color-text)]">
                {JSON.stringify(context, null, 2)}
              </pre>
            </details>
          )}
        </li>
        );
      })}
    </ul>
  );
}
