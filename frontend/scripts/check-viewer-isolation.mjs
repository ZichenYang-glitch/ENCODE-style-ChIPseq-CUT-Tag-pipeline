import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';

const dist = new URL('../dist/', import.meta.url);
const engine = /wuepgg|epigenome\.wustl\.edu/;

function files() {
  const found = new Map();
  function walk(relative) {
    const absolute = join(dist.pathname, relative);
    let entries;
    try {
      entries = readdirSync(absolute, { withFileTypes: true });
    } catch {
      found.set(relative, readFileSync(absolute));
      return;
    }
    for (const entry of entries) {
      const next = relative === '' ? entry.name : `${relative}/${entry.name}`;
      if (entry.isDirectory()) walk(next);
      else found.set(next, readFileSync(join(dist.pathname, next)));
    }
  }
  walk('');
  return found;
}

function references(bytes) {
  return [...bytes.toString('latin1').matchAll(/assets\/[A-Za-z0-9._~+-]+\.(?:js|css)/g)].map(
    (match) => match[0],
  );
}

const assets = files();
if (!assets.has('index.html') || !assets.has('viewer.html')) {
  throw new Error('platform and viewer documents must both be built');
}

const platform = new Set();
const pending = ['index.html'];
while (pending.length > 0) {
  const path = pending.pop();
  if (platform.has(path) || !assets.has(path)) continue;
  platform.add(path);
  for (const reference of references(assets.get(path))) pending.push(reference);
}

const leaked = [...platform].filter((path) => engine.test(assets.get(path).toString('latin1')));
if (leaked.length > 0) {
  throw new Error(`platform pages load the viewer engine: ${leaked.join(', ')}`);
}

const viewer = [...assets.entries()].filter(
  ([path]) => path === 'viewer.html' || path.startsWith('assets/'),
);
if (!viewer.some(([, bytes]) => engine.test(bytes.toString('latin1')))) {
  throw new Error('viewer build does not contain the WashU engine');
}

console.log(
  JSON.stringify({
    ok: true,
    platformFiles: platform.size,
    viewerDocument: 'viewer.html',
  }),
);
