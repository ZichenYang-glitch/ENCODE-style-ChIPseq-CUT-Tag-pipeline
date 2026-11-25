import { useCallback, useEffect, useState } from 'react';
import GenomeBrowserPanel from '../features/run-artifacts/GenomeBrowserPanel';
import { bindFramedViewer } from './framedSession';
import {
  VIEWER_PROTOCOL_VERSION,
  framedInstanceId,
  postViewerMessage,
  standaloneRunId,
  type ViewerConfigurePayload,
  type ViewerReferenceStatus,
} from './contract';

function sameConfigure(
  current: ViewerConfigurePayload | null,
  next: ViewerConfigurePayload,
): boolean {
  if (current === null) return false;
  const currentIds = current.artifactIds ?? [];
  const nextIds = next.artifactIds ?? [];
  return (
    current.runId === next.runId &&
    current.region === next.region &&
    current.display.showToolbar === next.display.showToolbar &&
    currentIds.length === nextIds.length &&
    currentIds.every((id, index) => id === nextIds[index])
  );
}

/**
 * Same-origin viewer document. An iframe isolates its document and styles
 * from the platform page. It is not a credential or process sandbox: reads
 * still use the platform session cookie and the authorized track API.
 */
export function ViewerPage() {
  const framed = window.parent !== window;
  const instanceId = framed ? framedInstanceId(window.location.hash) : null;
  const [command, setCommand] = useState<ViewerConfigurePayload | null>(null);
  const [closed, setClosed] = useState(false);

  useEffect(() => {
    if (!framed || instanceId === null) return undefined;
    return bindFramedViewer({
      page: window,
      host: window.parent,
      instanceId,
      onClose: () => {
        setClosed(true);
        setCommand(null);
      },
      onConfigure: (payload) => {
        setClosed(false);
        setCommand((current) => (sameConfigure(current, payload) ? current : payload));
      },
    });
  }, [framed, instanceId]);

  const publish = useCallback(
    (notice: { runId: string; referenceStatus: ViewerReferenceStatus; artifactIds: string[] } | { runId: string; message: string }) => {
      if (!framed || instanceId === null || window.parent === window) return;
      if ('message' in notice) {
        postViewerMessage(
          window.parent,
          {
            version: VIEWER_PROTOCOL_VERSION,
            type: 'helixweave.viewer.failed',
            instanceId,
            payload: notice,
          },
          window.location.origin,
        );
        return;
      }
      postViewerMessage(
        window.parent,
        {
          version: VIEWER_PROTOCOL_VERSION,
          type: 'helixweave.viewer.ready',
          instanceId,
          payload: notice,
        },
        window.location.origin,
      );
    },
    [framed, instanceId],
  );

  const runId = framed ? (command?.runId ?? null) : standaloneRunId(window.location.search);
  const reportFailure = useCallback(
    (message: string) => {
      if (runId !== null) publish({ runId, message });
    },
    [publish, runId],
  );
  // A framed document still waits for the host. A standalone document opens
  // with no run so local tracks and authorized runs can share one window.
  const waitingForHost = closed || (framed && command === null);
  const body = waitingForHost ? (
    <p role="status">Waiting for the platform viewer host.</p>
  ) : (
    <GenomeBrowserPanel
      runId={runId}
      initialRegion={command?.region ?? null}
      artifactIds={command?.artifactIds ?? null}
      showToolbar={command?.display.showToolbar ?? true}
      onReady={publish}
      onFailed={reportFailure}
    />
  );

  return (
    <main className="min-h-screen bg-[var(--color-bg)] p-3 text-[var(--color-text)]">
      {body}
      <p className="mt-4 text-xs text-[var(--color-text-muted)]">
        <a href="/third-party-notices/wuepgg-62.0.2-LICENSE.html">
          WashU Epigenome Browser is included for non-commercial use.
        </a>
      </p>
    </main>
  );
}
