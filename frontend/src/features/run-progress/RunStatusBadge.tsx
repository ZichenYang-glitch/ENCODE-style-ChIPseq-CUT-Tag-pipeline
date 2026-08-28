import { Badge, type BadgeTone } from '../../components/Badge';

interface RunStatusBadgeProps {
  status: string;
}

const statusTones: Record<string, BadgeTone> = {
  created: 'neutral',
  validating: 'info',
  planned: 'info',
  queued: 'info',
  running: 'info',
  succeeded: 'success',
  failed: 'error',
  cancelled: 'neutral',
};

export function RunStatusBadge({ status }: RunStatusBadgeProps) {
  const known = status in statusTones;
  return (
    <Badge
      tone={known ? statusTones[status] : 'neutral'}
      className="capitalize"
      data-testid="run-status-badge"
    >
      {known ? status : 'unknown'}
    </Badge>
  );
}
