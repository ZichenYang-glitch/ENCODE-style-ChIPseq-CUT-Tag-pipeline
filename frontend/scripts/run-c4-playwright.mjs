import { randomUUID } from 'node:crypto';
import { cpSync, mkdtempSync, mkdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const evidenceRoot = process.env.HELIXWEAVE_C4_EVIDENCE;
if (!evidenceRoot) {
  console.error('HELIXWEAVE_C4_EVIDENCE is required.');
  process.exit(1);
}
mkdirSync(evidenceRoot, { recursive: true });

const sentinelName = '.encode-platform-playwright-owned';
const runtimeOwner = randomUUID();
const runtimeRoot = mkdtempSync(join(tmpdir(), 'encode-platform-c4-'));
const sentinelPath = join(runtimeRoot, sentinelName);
writeFileSync(sentinelPath, runtimeOwner, 'utf8');
const junitPath = join(evidenceRoot, 'c4-junit.xml');

const cliPath = fileURLToPath(
  new URL('../node_modules/@playwright/test/cli.js', import.meta.url),
);
const runtimeHelperPath = fileURLToPath(
  new URL('../../test/browser/platform_runtime.py', import.meta.url),
);
const configPath = fileURLToPath(
  new URL('../playwright.c4.config.ts', import.meta.url),
);
const child = spawn(
  process.execPath,
  [cliPath, 'test', '--config', configPath, ...process.argv.slice(2)],
  {
    env: {
      ...process.env,
      ENCODE_PIPELINE_E2E_ROOT: runtimeRoot,
      ENCODE_PIPELINE_E2E_OWNER: runtimeOwner,
      HELIXWEAVE_C4_NATIVE: '1',
      HELIXWEAVE_C4_JUNIT: junitPath,
      HELIXWEAVE_C4_EVIDENCE: evidenceRoot,
    },
    stdio: 'inherit',
  },
);

for (const signal of ['SIGINT', 'SIGTERM']) {
  process.on(signal, () => child.kill(signal));
}

function preserveRuntimeEvidence() {
  const logs = join(runtimeRoot, 'logs');
  try {
    cpSync(logs, join(evidenceRoot, 'logs'), { recursive: true, force: true });
  } catch {
    console.error('Could not copy the C4 runtime logs.');
  }
}

child.once('error', (error) => {
  console.error(error);
  process.exitCode = 1;
});

child.once('exit', (code, signal) => {
  preserveRuntimeEvidence();
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
    console.error('Could not clean the owned C4 runtime.');
    process.exitCode = 1;
    return;
  }
  if (signal) {
    process.exitCode = signal === 'SIGINT' ? 130 : 143;
    return;
  }
  process.exitCode = code ?? 1;
});
