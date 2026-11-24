import { useId, useState, type FormEvent } from 'react';
import { useQuery } from '@tanstack/react-query';
import { Link } from 'react-router-dom';
import { ApiError, fetcher } from '../../api/fetcher';
import { getRun } from '../../api/generated/runs/runs';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';

type Target =
  | { action: 'create'; snapshot_id: string; tags: Record<string, string> }
  | { action: 'cancel'; run_id: string };

interface Confirmation {
  target: Target;
  workflow_id: string;
  confirmation: string;
}

interface RunActionResult {
  run_id: string;
  stop_requested: boolean;
}

const terminal = new Set(['succeeded', 'failed', 'cancelled']);
const controlClass =
  'mt-1 w-full rounded border border-[var(--color-border)] bg-[var(--color-surface)] p-2';

export function PeriRunActions({
  sid,
  connected,
  snapshotId,
  disabled,
}: {
  sid: string;
  connected: boolean;
  snapshotId: string | null;
  disabled: boolean;
}) {
  const actionGroup = useId();
  const [action, setAction] = useState<Target['action']>('create');
  const [runId, setRunId] = useState('');
  const [tagLines, setTagLines] = useState('');
  const [preview, setPreview] = useState<Confirmation | null>(null);
  const [result, setResult] = useState<RunActionResult | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [retry, setRetry] = useState(false);

  if (!connected && result === null) return null;

  const root = `/api/v1/agent/authoring-sessions/${encodeURIComponent(sid)}/run-actions`;
  const stale = preview?.target.action === 'create' &&
    preview.target.snapshot_id !== snapshotId;

  async function review(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (busy || disabled || !connected) return;
    setError(null);
    setRetry(false);
    let target: Target;
    if (action === 'create') {
      if (snapshotId === null) return;
      const tags = parseTags(tagLines);
      if (tags === null) {
        setError('Enter each distinct tag as key=value.');
        return;
      }
      target = { action, snapshot_id: snapshotId, tags };
    } else {
      target = { action, run_id: runId.trim() };
    }
    setBusy(true);
    try {
      setPreview(await post<Confirmation>(`${root}/preview`, { target }));
      setResult(null);
    } catch (failure) {
      setError(message(failure, 'Could not load the confirmation.'));
    } finally {
      setBusy(false);
    }
  }

  async function confirm() {
    if (!preview || busy || disabled || !connected || stale) return;
    setError(null);
    setBusy(true);
    try {
      setResult(await post<RunActionResult>(`${root}/confirm`, {
        target: preview.target,
        confirmation: preview.confirmation,
      }));
      setPreview(null);
      setRetry(false);
    } catch (failure) {
      setRetry(!(failure instanceof ApiError) || failure.status >= 500);
      setError(message(failure, 'The action response was lost. Retry this confirmation.'));
    } finally {
      setBusy(false);
    }
  }

  function changeTarget() {
    setPreview(null);
    setError(null);
    setRetry(false);
  }

  return (
    <Panel title="Peri run actions">
      <div className="space-y-3 text-sm">
        {connected && (preview === null ? (
          <form className="space-y-3" onSubmit={(event) => void review(event)}>
            <fieldset className="space-y-2" disabled={busy || disabled}>
              <legend className="font-medium">Request a run action</legend>
              <label className="flex items-center gap-2">
                <input type="radio" name={actionGroup} checked={action === 'create'}
                  onChange={() => setAction('create')} />
                Create a run from validated inputs
              </label>
              <label className="flex items-center gap-2">
                <input type="radio" name={actionGroup} checked={action === 'cancel'}
                  onChange={() => setAction('cancel')} />
                Cancel a specific run
              </label>
              {action === 'create' ? (
                snapshotId ? (
                  <label className="block">
                    Tags (optional, one key=value per line)
                    <textarea className={controlClass} rows={2} value={tagLines}
                      onChange={(event) => setTagLines(event.target.value)} />
                  </label>
                ) : <p>Validate this draft to obtain a snapshot.</p>
              ) : (
                <label className="block">
                  Run ID
                  <input className={controlClass} value={runId} required maxLength={128}
                    pattern="[A-Za-z0-9][A-Za-z0-9_.\x2d]{0,127}"
                    onChange={(event) => setRunId(event.target.value)} />
                </label>
              )}
              <Button type="submit" variant="secondary" disabled={action === 'create' && !snapshotId}>
                Review run action
              </Button>
            </fieldset>
          </form>
        ) : (
          <div className="space-y-2" aria-label="Review Peri run action" role="group">
            <p className="font-medium">{preview.target.action === 'create' ? 'Create run' : 'Cancel run'}</p>
            <p>Workflow <code>{preview.workflow_id}</code></p>
            {preview.target.action === 'create' ? (
              <>
                <p>Snapshot <code className="break-all">{preview.target.snapshot_id}</code></p>
                <p>Tags <code className="break-all">{JSON.stringify(preview.target.tags)}</code></p>
                <p>The run will be created without preflight or execution.</p>
              </>
            ) : (
              <p>Run <code className="break-all">{preview.target.run_id}</code></p>
            )}
            {stale && <p>Review the current snapshot for a new confirmation.</p>}
            <div className="flex flex-wrap gap-2">
              <Button
                type="button"
                variant={preview.target.action === 'cancel' ? 'danger' : 'primary'}
                disabled={busy || disabled || stale || (error !== null && !retry)}
                onClick={() => void confirm()}
              >
                {retry ? 'Retry confirmation' : preview.target.action === 'create'
                  ? 'Confirm run creation' : 'Confirm cancellation'}
              </Button>
              <Button type="button" variant="secondary" disabled={busy} onClick={changeTarget}>
                Change request
              </Button>
            </div>
          </div>
        ))}
        {error && (connected || result) && <p role="alert">{error}</p>}
        {result && <ConfirmedRun result={result} />}
      </div>
    </Panel>
  );
}

function ConfirmedRun({ result }: { result: RunActionResult }) {
  const readback = useQuery({
    queryKey: ['peri-confirmed-run', result.run_id],
    queryFn: ({ signal }) => getRun(result.run_id, { signal }),
    retry: false,
    refetchOnMount: 'always',
    refetchInterval: ({ state }) =>
      result.stop_requested &&
      !(state.error instanceof ApiError && [401, 403, 404].includes(state.error.status)) &&
      !terminal.has(state.data?.run?.status ?? '') ? 1500 : false,
  });
  const status = readback.data?.run?.status;
  const waitingForStop = result.stop_requested && !terminal.has(status ?? '');
  return (
    <div className="space-y-2">
      <p>
        <Link to={`/runs/${encodeURIComponent(result.run_id)}`} className="underline">
          Open run {result.run_id}
        </Link>
      </p>
      <div role="status">
        {waitingForStop && <p>Stop requested.</p>}
        {readback.isError ? (
          <p>Could not refresh the run status.</p>
        ) : status ? <p>Run status: {status}.</p> : <p>Reading run status…</p>}
      </div>
      {readback.isError && (
        <Button type="button" variant="secondary" onClick={() => void readback.refetch()}>
          Refresh run status
        </Button>
      )}
    </div>
  );
}

function post<T>(path: string, body: unknown): Promise<T> {
  return fetcher<T>(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
}

function parseTags(input: string): Record<string, string> | null {
  const tags = new Map<string, string>();
  for (const line of input.split('\n')) {
    if (!line.trim()) continue;
    const separator = line.indexOf('=');
    if (separator < 1) return null;
    const key = line.slice(0, separator).trim();
    if (!key || tags.has(key)) return null;
    tags.set(key, line.slice(separator + 1).trim());
  }
  return Object.fromEntries(tags);
}

function message(failure: unknown, fallback: string): string {
  return failure instanceof ApiError ? failure.message : fallback;
}
