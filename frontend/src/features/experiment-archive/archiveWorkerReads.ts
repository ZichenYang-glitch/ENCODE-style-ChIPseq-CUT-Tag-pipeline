import workerFetchSource from './archiveWorkerFetch.js?raw';
import { archiveErrorCopy, archiveHttpError } from './archiveView';

export function observeArchiveWorkerReads(
  contentUrl: string,
  indexUrl: string | null,
  onError: (failure: ReturnType<typeof archiveErrorCopy>) => void,
) {
  const createUrl = URL.createObjectURL;
  if (typeof createUrl !== 'function' || typeof Worker === 'undefined') return () => {};
  const channelName = 'helixweave-archive-' + crypto.getRandomValues(new Uint32Array(4)).join('-');
  const channel = new BroadcastChannel(channelName);
  channel.onmessage = ({ data }) => {
    if (!data || !Number.isInteger(data.status)) return;
    onError(archiveHttpError(data.status, typeof data.code === 'string' ? data.code : null));
  };
  const config = { channelName, origin: window.location.origin, members: [contentUrl, indexUrl] };
  const bootstrap = `{\n${workerFetchSource}\ninstallArchiveWorkerFetch(${JSON.stringify(config)});\n}\n`;

  // WashU 62.0.2 creates its classic data workers from JavaScript blobs and
  // immediately revokes their URLs. Supply the member-scoped transport before
  // creation so the worker's range reader and its proxy retry share it.
  const createArchiveUrl: typeof URL.createObjectURL = (object) => {
    const script = object instanceof Blob && object.type.split(';')[0] === 'text/javascript';
    return createUrl.call(URL, script ? new Blob([bootstrap, object], { type: object.type }) : object);
  };
  URL.createObjectURL = createArchiveUrl;
  return () => {
    if (URL.createObjectURL === createArchiveUrl) URL.createObjectURL = createUrl;
    channel.close();
  };
}
