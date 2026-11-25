import { randomUUID } from 'node:crypto';

// Installed by the browser proof only. A slice is a view; count bytes when a
// browser read completes. Include repeat reads of the same range.
function measureReads(report, scope) {
  const owners = new WeakMap();
  const ownerOf = (blob) => blob instanceof File ? blob.name : owners.get(blob);
  const slice = Blob.prototype.slice;
  Blob.prototype.slice = function (...args) {
    const view = slice.apply(this, args);
    const name = ownerOf(this);
    if (name) owners.set(view, name);
    return view;
  };
  for (const method of ['arrayBuffer', 'bytes', 'text']) {
    const read = Blob.prototype[method];
    if (!read) continue;
    Blob.prototype[method] = async function (...args) {
      const value = await read.apply(this, args);
      const name = ownerOf(this);
      if (name) report({ name, bytes: value.byteLength ?? this.size, scope, method });
      return value;
    };
  }
  const stream = Blob.prototype.stream;
  Blob.prototype.stream = function (...args) {
    const name = ownerOf(this);
    const input = stream.apply(this, args);
    if (!name) return input;
    return input.pipeThrough(new TransformStream({
      transform(chunk, controller) {
        report({ name, bytes: chunk.byteLength, scope, method: 'stream' });
        controller.enqueue(chunk);
      },
    }));
  };
  for (const method of ['readAsArrayBuffer', 'readAsText']) {
    const read = FileReader.prototype[method];
    FileReader.prototype[method] = function (blob, ...args) {
      const name = ownerOf(blob);
      if (name) {
        this.addEventListener('loadend', () => {
          if (this.result !== null && this.error === null) {
            report({ name, bytes: blob.size, scope, method });
          }
        }, { once: true });
      }
      return read.call(this, blob, ...args);
    };
  }
}

export async function installFileReadMeter(context) {
  const records = [];
  const channel = `washu-test-${randomUUID()}`;
  await context.exposeBinding('__recordWashuFileRead', (_, record) => { records.push(record); });
  // The pinned bundle creates and immediately revokes a Blob URL for each
  // classic worker. Prepend the same observer before any vendor code runs;
  // do not change the worker's messages, file handles, or byte results.
  const workerCode = `(function(){const out = new BroadcastChannel(${JSON.stringify(channel)});(${measureReads})(r => out.postMessage(r), 'worker');})();\n`;
  await context.addInitScript({ content: `(() => {
    const channel = new BroadcastChannel(${JSON.stringify(channel)});
    channel.onmessage = (event) => window.__recordWashuFileRead(event.data);
    (${measureReads})(r => window.__recordWashuFileRead(r), 'document');
    const original = URL.createObjectURL.bind(URL);
    URL.createObjectURL = function (blob) {
      if (blob instanceof Blob && !(blob instanceof File) &&
          /^(text|application)\\/javascript(?:;|$)/i.test(blob.type)) {
        return original(new Blob([${JSON.stringify(workerCode)}, blob], {type:blob.type}));
      }
      return original(blob);
    };
  })();` });
  return {
    mark: () => records.length,
    since: (mark) => records.slice(mark),
  };
}
