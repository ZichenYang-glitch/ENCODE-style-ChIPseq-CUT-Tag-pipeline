import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const source = readFileSync(
  join(here, '../../frontend/src/api/generated/catalog/catalog.ts'),
  'utf8',
);

function extractFunction(name) {
  const marker = `export const ${name} = `;
  const start = source.indexOf(marker);
  if (start < 0) {
    throw new Error(`generated client is missing ${name}`);
  }
  const arrow = source.indexOf('=>', start);
  const brace = source.indexOf('{', arrow);
  let depth = 0;
  let end = brace;
  for (; end < source.length; end += 1) {
    const character = source[end];
    if (character === '{') {
      depth += 1;
    } else if (character === '}') {
      depth -= 1;
      if (depth === 0) {
        end += 1;
        break;
      }
    }
  }
  return source.slice(start, end).replace(
    /=\s*\(projectId:\s*string,\s*params\?:\s*[A-Za-z0-9_]+,\)/,
    '= (projectId, params)',
  );
}

const directory = mkdtempSync(join(tmpdir(), 'catalog-client-'));
const modulePath = join(directory, 'catalog.mjs');
writeFileSync(
  modulePath,
  `${extractFunction('getCatalogListExperimentsUrl')}\n${extractFunction('getCatalogListSamplesUrl')}\n`,
);
const client = await import(pathToFileURL(modulePath).href);
const calls = JSON.parse(process.argv[2]);
const urls = calls.map((call) => {
  const build = client[call.fn];
  if (typeof build !== 'function') {
    throw new Error(`generated client is missing ${call.fn}`);
  }
  return build(...call.args);
});
process.stdout.write(JSON.stringify(urls));
