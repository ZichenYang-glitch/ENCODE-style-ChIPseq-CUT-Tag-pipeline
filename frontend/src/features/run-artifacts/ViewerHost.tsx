import { useEffect, useMemo, useRef, useState } from 'react';
import { Button } from '../../components/Button';
import {
  closeViewerMessage,
  configureViewerMessage,
  newViewerInstanceId,
  postViewerMessage,
  readViewerNotice,
} from '../../viewer/contract';

/** Failure bound for a frame that never confirms its listener. Not a send window. */
export const VIEWER_HANDSHAKE_TIMEOUT_MS = 20_000;

/**
 * Thin same-origin host. The iframe gives the viewer its own document and
 * styles. It does not create a strong security sandbox.
 *
 * Configure is sent only after that document announces that its handler is
 * installed. Receipt is the child's accepted notice, or a later ready or
 * failed notice. Replacing the instance closes the retired frame and stops
 * its timers.
 */
export function ViewerHost({ runId }: { runId: string }) {
  const [attempt, setAttempt] = useState(0);
  const instanceId = useMemo(() => newViewerInstanceId(), [runId, attempt]);
  return (
    <section className="space-y-2" data-testid="viewer-host">
      <div className="flex justify-end">
        <a
          className="text-sm text-[var(--color-link)] underline"
          href={`/viewer.html?run=${encodeURIComponent(runId)}`}
          target="_blank"
          rel="noopener noreferrer"
        >
          Open genome browser full page
        </a>
      </div>
      <ViewerFrame
        key={instanceId}
        runId={runId}
        instanceId={instanceId}
        onRetry={() => setAttempt((value) => value + 1)}
      />
    </section>
  );
}

function ViewerFrame({
  runId,
  instanceId,
  onRetry,
}: {
  runId: string;
  instanceId: string;
  onRetry: () => void;
}) {
  const frameRef = useRef<HTMLIFrameElement>(null);
  const [status, setStatus] = useState<'loading' | 'ready' | 'failed'>('loading');
  const [failure, setFailure] = useState<'timeout' | 'viewer' | null>(null);

  useEffect(() => {
    const owned = frameRef.current?.contentWindow ?? null;
    let alive = true;
    let handshakeOk = false;
    let currentEpoch: string | null = null;
    let phase: 'loading' | 'ready' | 'failed' = 'loading';
    let reason: 'timeout' | 'viewer' | null = null;
    let timer = 0;

    const publish = (
      next: 'loading' | 'ready' | 'failed',
      nextReason: 'timeout' | 'viewer' | null,
    ) => {
      phase = next;
      reason = nextReason;
      setStatus(next);
      setFailure(nextReason);
    };
    const sendConfigure = () => {
      const frame = frameRef.current?.contentWindow;
      if (frame == null) return;
      postViewerMessage(
        frame,
        configureViewerMessage(instanceId, {
          runId,
          region: null,
          artifactIds: null,
          display: { showToolbar: true },
        }),
        window.location.origin,
      );
    };
    const armTimeout = () => {
      window.clearTimeout(timer);
      timer = window.setTimeout(() => {
        if (!alive || handshakeOk) return;
        publish('failed', 'timeout');
      }, VIEWER_HANDSHAKE_TIMEOUT_MS);
    };
    const onMessage = (event: MessageEvent) => {
      if (!alive) return;
      const frame = frameRef.current?.contentWindow ?? null;
      const notice = readViewerNotice(event.data, instanceId, event, {
        origin: window.location.origin,
        source: frame,
      });
      if (notice === null) return;
      if (notice.type === 'helixweave.viewer.listening') {
        const epoch = notice.payload.epoch;
        const reloaded = currentEpoch !== null && epoch !== currentEpoch;
        currentEpoch = epoch;
        if (reloaded) {
          handshakeOk = false;
          publish('loading', null);
          armTimeout();
        }
        if (!handshakeOk) sendConfigure();
        return;
      }
      if (notice.type === 'helixweave.viewer.accepted') {
        if (notice.payload.runId !== runId) return;
        handshakeOk = true;
        window.clearTimeout(timer);
        if (phase === 'failed' && reason === 'timeout') publish('loading', null);
        return;
      }
      if (notice.payload.runId !== runId) return;
      handshakeOk = true;
      window.clearTimeout(timer);
      if (notice.type === 'helixweave.viewer.ready') publish('ready', null);
      else publish('failed', 'viewer');
    };

    window.addEventListener('message', onMessage);
    armTimeout();
    return () => {
      alive = false;
      window.removeEventListener('message', onMessage);
      window.clearTimeout(timer);
      const current = frameRef.current?.contentWindow ?? null;
      if (owned !== null && owned !== current) {
        try {
          postViewerMessage(owned, closeViewerMessage(instanceId), window.location.origin);
        } catch {
          // The retired document is already gone.
        }
      }
    };
  }, [instanceId, runId]);

  return (
    <>
      <iframe
        title="Genome browser"
        src={`/viewer.html#instance=${instanceId}`}
        className="h-[40rem] w-full border border-[var(--color-border)] bg-white"
        ref={frameRef}
        onError={() => {
          setStatus('failed');
          setFailure('viewer');
        }}
      />
      <div className="flex items-center justify-between gap-3">
        <p role="status" className="text-xs text-[var(--color-text-muted)]">
          {status === 'ready'
            ? 'Genome browser is using the authorized track API.'
            : status === 'failed'
              ? failure === 'timeout'
                ? 'Genome browser did not finish opening.'
                : 'Genome browser could not load this run.'
              : 'Opening genome browser.'}
        </p>
        {status === 'failed' ? (
          <Button type="button" variant="quiet" size="sm" onClick={onRetry}>
            Retry genome browser
          </Button>
        ) : null}
      </div>
    </>
  );
}
