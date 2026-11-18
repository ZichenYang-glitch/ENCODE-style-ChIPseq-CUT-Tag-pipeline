import {
  acceptedViewerMessage,
  listeningViewerMessage,
  postViewerMessage,
  readHostCommand,
  type ViewerConfigurePayload,
} from './contract';

export const VIEWER_LISTENING_INTERVAL_MS = 250;

/**
 * Install the framed page's command handler, then announce that it is
 * listening. The announcement is repeated until configure or close arrives
 * so a host that was not ready for the first post still delivers configure.
 * Calling the returned function stops both the handler and the announcement.
 */
export function bindFramedViewer(options: {
  page: Window;
  host: Window;
  instanceId: string;
  onConfigure: (payload: ViewerConfigurePayload) => void;
  onClose: () => void;
}): () => void {
  const { page, host, instanceId, onConfigure, onClose } = options;
  if (host === page) return () => {};
  let stopped = false;
  let timer: ReturnType<Window['setInterval']> | undefined;
  const origin = page.location.origin;
  const epoch = crypto.randomUUID();
  const stopAnnouncing = () => {
    stopped = true;
    if (timer !== undefined) page.clearInterval(timer);
  };
  const onMessage = (event: MessageEvent) => {
    const command = readHostCommand(event.data, instanceId, event, {
      origin,
      source: host,
    });
    if (command === null) return;
    stopAnnouncing();
    if (command.type === 'helixweave.viewer.close') {
      onClose();
      return;
    }
    onConfigure(command.payload);
    postViewerMessage(host, acceptedViewerMessage(instanceId, command.payload.runId), origin);
  };
  page.addEventListener('message', onMessage);
  const announce = () => {
    if (stopped || host === page) return;
    postViewerMessage(host, listeningViewerMessage(instanceId, epoch), origin);
  };
  announce();
  if (!stopped) timer = page.setInterval(announce, VIEWER_LISTENING_INTERVAL_MS);
  return () => {
    stopAnnouncing();
    page.removeEventListener('message', onMessage);
  };
}
