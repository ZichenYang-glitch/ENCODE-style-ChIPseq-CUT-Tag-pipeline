// Runs before the pinned WashU reader in its classic worker. Only the error
// status and public issue code cross back to the archive page.
function installArchiveWorkerFetch({ origin, members, channelName }) {
  const channel = new BroadcastChannel(channelName);
  const fetchBytes = self.fetch.bind(self);
  self.fetch = async (input, init) => {
    const request = input instanceof Request ? input : null;
    const url = new URL(request?.url ?? String(input), origin);
    if (url.origin !== origin || !members.includes(url.href)) {
      throw new TypeError('This archive reader only opens its sealed members.');
    }
    const signal = init?.signal ?? request?.signal;
    let response;
    try {
      response = await fetchBytes(input, {
        ...init, mode: 'same-origin', credentials: 'same-origin', cache: 'no-store', redirect: 'error',
      });
    } catch (error) {
      if (!signal?.aborted && error?.name !== 'AbortError') {
        channel.postMessage({ status: 0, code: 'archive_offline' });
      }
      throw error;
    }
    if (!response.ok && !signal?.aborted) {
      const body = await response.clone().json().catch(() => null);
      const code = body?.issues?.[0]?.code;
      channel.postMessage({ status: response.status, code: typeof code === 'string' ? code : null });
    }
    return response;
  };
}
