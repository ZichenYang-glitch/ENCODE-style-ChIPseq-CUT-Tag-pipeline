import { useEffect, useMemo, useRef, useState } from 'react';
import {
  configureViewerMessage,
  newViewerInstanceId,
  postViewerMessage,
  readViewerNotice,
} from '../../viewer/contract';

/**
 * Thin same-origin host. The iframe gives the viewer its own document and
 * styles. It does not create a strong security sandbox.
 */
export function ViewerHost({ runId }: { runId: string }) {
  const instanceId = useMemo(() => newViewerInstanceId(), [runId]);
  const frameRef = useRef<HTMLIFrameElement>(null);
  const [status, setStatus] = useState<'loading' | 'ready' | 'failed'>('loading');

  useEffect(() => {
    const onMessage = (event: MessageEvent) => {
      const frame = frameRef.current?.contentWindow ?? null;
      const notice = readViewerNotice(event.data, instanceId, event, {
        origin: window.location.origin,
        source: frame,
      });
      if (notice === null || notice.payload.runId !== runId) return;
      setStatus(notice.type === 'helixweave.viewer.ready' ? 'ready' : 'failed');
    };
    window.addEventListener('message', onMessage);
    return () => window.removeEventListener('message', onMessage);
  }, [instanceId, runId]);

  useEffect(() => {
    const send = () => {
      const frame = frameRef.current?.contentWindow;
      if (frame === null || frame === undefined) return;
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
    send();
    const retry = window.setInterval(send, 250);
    const stop = window.setTimeout(() => window.clearInterval(retry), 2_000);
    return () => {
      window.clearInterval(retry);
      window.clearTimeout(stop);
    };
  }, [instanceId, runId]);

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
      <iframe
        key={instanceId}
        title="Genome browser"
        src={`/viewer.html#instance=${instanceId}`}
        className="h-[40rem] w-full border border-[var(--color-border)] bg-white"
        ref={frameRef}
      />
      <p role="status" className="text-xs text-[var(--color-text-muted)]">
        {status === 'ready'
          ? 'Genome browser is using the authorized track API.'
          : status === 'failed'
            ? 'Genome browser could not load this run.'
            : 'Opening genome browser.'}
      </p>
    </section>
  );
}
