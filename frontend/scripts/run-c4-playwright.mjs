import { execFileSync, spawn, spawnSync } from 'node:child_process';
import { createHash, randomUUID } from 'node:crypto';
import { cpSync, mkdtempSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

const evidenceRoot = process.env.HELIXWEAVE_C4_EVIDENCE;
if (!evidenceRoot) {
  console.error('HELIXWEAVE_C4_EVIDENCE is required.');
  process.exit(1);
}
mkdirSync(evidenceRoot, { recursive: true });
const repositoryRoot = fileURLToPath(new URL('../..', import.meta.url));
const git = (args) =>
  execFileSync('git', args, { cwd: repositoryRoot, encoding: 'utf8' }).trim();
const status = execFileSync('git', ['status', '--porcelain'], {
  cwd: repositoryRoot,
  encoding: 'utf8',
});
if (process.env.HELIXWEAVE_C4_REQUIRE_CLEAN === '1' && status.trim()) {
  console.error('C4 requires a clean candidate worktree.');
  process.exit(1);
}
const binaryPath = process.env.HELIXWEAVE_PERI_BINARY;
const binarySha = binaryPath
  ? createHash('sha256').update(readFileSync(binaryPath)).digest('hex')
  : null;
writeFileSync(
  join(evidenceRoot, 'identity.json'),
  `${JSON.stringify(
    {
      git_sha: git(['rev-parse', 'HEAD']),
      git_tree: git(['rev-parse', 'HEAD^{tree}']),
      worktree_status: status,
      served_frontend: 'vite-dev',
      served_frontend_note:
        'The C4 harness starts scripts/run_local_platform.py, which serves this source tree with Vite. It does not serve the packaged production bundle.',
      binary_sha256: binarySha,
    },
    null,
    2,
  )}\n`,
);

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
  writeFileSync(
    join(evidenceRoot, 'cleanup.json'),
    `${JSON.stringify(
      {
        status: cleanup.status,
        signal: cleanup.signal,
        error: cleanup.error ? String(cleanup.error) : null,
        stderr: cleanup.stderr,
      },
      null,
      2,
    )}\n`,
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
