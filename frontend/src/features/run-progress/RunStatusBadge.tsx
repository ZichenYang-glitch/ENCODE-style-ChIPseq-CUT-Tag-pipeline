interface RunStatusBadgeProps {
  status: string;
}

const statusClasses: Record<string, string> = {
  created: 'bg-gray-100 text-gray-700 border-gray-300',
  validating: 'bg-blue-50 text-blue-700 border-blue-300',
  planned: 'bg-blue-50 text-blue-700 border-blue-300',
  queued: 'bg-blue-50 text-blue-700 border-blue-300',
  running: 'bg-blue-50 text-blue-700 border-blue-300',
  succeeded: 'bg-green-50 text-green-700 border-green-300',
  failed: 'bg-red-50 text-red-700 border-red-300',
  cancelled: 'bg-yellow-50 text-yellow-700 border-yellow-300',
};

export function RunStatusBadge({ status }: RunStatusBadgeProps) {
  const known = status in statusClasses;
  const classes = known
    ? statusClasses[status]
    : 'bg-gray-50 text-gray-600 border-gray-300';
  return (
    <span
      className={`inline-flex items-center rounded border px-2 py-0.5 text-xs font-semibold ${classes}`}
      data-testid="run-status-badge"
    >
      {known ? status : 'unknown'}
    </span>
  );
}
