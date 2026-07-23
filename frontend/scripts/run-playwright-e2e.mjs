import { randomUUID } from 'node:crypto';
import { mkdtempSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const sentinelName = '.encode-platform-playwright-owned';
const runtimeOwner = randomUUID();
const runtimeRoot = mkdtempSync(join(tmpdir(), 'encode-platform-playwright-'));
const sentinelPath = join(runtimeRoot, sentinelName);
writeFileSync(sentinelPath, runtimeOwner, 'utf8');

const cliPath = fileURLToPath(
  new URL('../node_modules/@playwright/test/cli.js', import.meta.url),
);
const runtimeHelperPath = fileURLToPath(
  new URL('../../test/browser/platform_runtime.py', import.meta.url),
);
const child = spawn(process.execPath, [cliPath, 'test', ...process.argv.slice(2)], {
  env: {
    ...process.env,
    ENCODE_PIPELINE_E2E_ROOT: runtimeRoot,
    ENCODE_PIPELINE_E2E_OWNER: runtimeOwner,
  },
  stdio: 'inherit',
});

for (const signal of ['SIGINT', 'SIGTERM']) {
  process.on(signal, () => child.kill(signal));
}

child.once('error', (error) => {
  console.error(error);
  process.exitCode = 1;
});

child.once('exit', (code, signal) => {
  if (code === 0) {
    const cleanup = spawnSync(
      'python3',
      [runtimeHelperPath, '--cleanup-owned-runtime'],
      {
        encoding: 'utf8',
        env: {
          ...process.env,
          ENCODE_PIPELINE_E2E_ROOT: runtimeRoot,
          ENCODE_PIPELINE_E2E_OWNER: runtimeOwner,
        },
        stdio: ['ignore', 'ignore', 'pipe'],
        timeout: 30_000,
      },
    );
    if (cleanup.error || cleanup.status !== 0) {
      console.error('Could not clean the owned Playwright runtime.');
      process.exitCode = 1;
      return;
    }
  }
  if (signal) {
    process.exitCode = signal === 'SIGINT' ? 130 : 143;
    return;
  }
  process.exitCode = code ?? 1;
});
