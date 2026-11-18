import { createRequire } from 'node:module';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const require = createRequire(import.meta.url);
const ts = require(join(here, '../../frontend/node_modules/typescript/lib/typescript.js'));
const source = readFileSync(
  join(here, '../../frontend/src/api/generated/catalog/catalog.ts'),
  'utf8',
);
const stripped = source.replace(/(^|\n)import[\s\S]*?;/g, '$1');
const transpiled = ts.transpileModule(stripped, {
  compilerOptions: {
    module: ts.ModuleKind.ES2022,
    target: ts.ScriptTarget.ES2022,
  },
}).outputText;
const directory = mkdtempSync(join(tmpdir(), 'catalog-client-'));
const modulePath = join(directory, 'catalog.mjs');
writeFileSync(modulePath, transpiled);
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
