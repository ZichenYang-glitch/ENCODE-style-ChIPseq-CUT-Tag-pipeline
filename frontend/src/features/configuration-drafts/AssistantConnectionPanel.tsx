import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import { connectBlockReason } from './actionGate';
import type { AssistantConnectionView } from './model';
import { connectionDetail, connectionHeading } from './statusCopy';

interface AssistantConnectionPanelProps {
  connection: AssistantConnectionView;
  writesEnabled: boolean;
  onConnect?: () => void;
}

export function AssistantConnectionPanel({
  connection,
  writesEnabled,
  onConnect,
}: AssistantConnectionPanelProps) {
  const connectReason = connectBlockReason(writesEnabled, connection.state);
  const detail = connectionDetail(connection);
  return (
    <Panel title="Assistant connection">
      <div className="space-y-3">
        <p role="status" aria-live="polite">
          <span className="block text-sm font-medium">
            {connectionHeading(connection.state)}
          </span>
          <span className="mt-1 block text-sm text-[var(--color-text-muted)]">
            {detail}
          </span>
        </p>
        <Button
          type="button"
          variant="primary"
          disabled={connectReason !== null || !onConnect}
          aria-describedby={
            connectReason !== null ? 'assistant-connect-reason' : undefined
          }
          onClick={() => {
            if (connectReason !== null) return;
            onConnect?.();
          }}
        >
          Connect assistant
        </Button>
        {connectReason !== null && (
          <p id="assistant-connect-reason" className="text-xs text-[var(--color-text-muted)]">
            {connectReason}
          </p>
        )}
      </div>
    </Panel>
  );
}
