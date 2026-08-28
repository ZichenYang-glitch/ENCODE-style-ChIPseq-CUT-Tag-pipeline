import type {
  WorkflowAvailability,
  WorkflowExecutionAvailability,
} from '../../api/types';

export function executionAvailabilityLabel(
  availability: WorkflowExecutionAvailability | null,
): string {
  if (availability === 'available') return 'Runnable';
  if (availability === 'not_configured') return 'Not configured';
  if (availability === 'unavailable') return 'Unavailable';
  return 'Availability unknown';
}

export function executionAvailabilityMessage(
  availability: WorkflowExecutionAvailability | null,
): string {
  if (availability === 'available') {
    return 'Execution runtime is available.';
  }
  if (availability === 'not_configured') {
    return 'Execution runtime is not configured. You can author and validate inputs, but cannot create or start a run.';
  }
  if (availability === 'unavailable') {
    return 'Execution runtime is unavailable. You can author and validate inputs, but cannot create or start a run.';
  }
  return 'Execution availability could not be confirmed. Run creation and start remain disabled.';
}

export function ExecutionAvailabilityBadge({
  availability,
}: {
  availability: WorkflowExecutionAvailability | null;
}) {
  const className =
    availability === 'available'
      ? 'border-[var(--color-success-border)] bg-[var(--color-success-bg)] text-[var(--color-success)]'
      : availability === 'not_configured'
        ? 'border-[var(--color-warning-border)] bg-[var(--color-warning-bg)] text-[var(--color-warning)]'
        : 'border-[var(--color-error-border)] bg-[var(--color-error-bg)] text-[var(--color-error)]';
  return (
    <span
      className={`inline-flex w-fit rounded-[4px] border px-2 py-0.5 text-xs font-medium ${className}`}
      data-testid="execution-availability-badge"
    >
      {executionAvailabilityLabel(availability)}
    </span>
  );
}

export function ExecutionAvailabilityNotice({
  availability,
  showAvailable = false,
}: {
  availability: WorkflowAvailability | WorkflowExecutionAvailability | null;
  showAvailable?: boolean;
}) {
  const execution =
    typeof availability === 'string'
      ? availability
      : availability?.execution ?? null;
  if (execution === 'available' && !showAvailable) return null;
  const ready = execution === 'available';
  return (
    <div
      className={`border-l-[3px] px-3 py-2 text-sm ${
        ready
          ? 'border-[var(--color-success)] bg-[var(--color-success-bg)] text-[var(--color-success)]'
          : 'border-[var(--color-warning)] bg-[var(--color-warning-bg)] text-[var(--color-warning)]'
      }`}
      role="status"
      data-testid="execution-availability-notice"
    >
      {executionAvailabilityMessage(execution)}
    </div>
  );
}
