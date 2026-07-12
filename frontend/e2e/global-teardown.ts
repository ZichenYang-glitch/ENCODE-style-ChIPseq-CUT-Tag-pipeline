import { existsSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';

const runtimeRoot = process.env.ENCODE_PIPELINE_E2E_ROOT!;

function processExists(pid: number): boolean {
  try {
    process.kill(pid, 0);
  } catch (error) {
    return (error as NodeJS.ErrnoException).code !== 'ESRCH';
  }
  return true;
}

function signalProcess(pid: number, signal: NodeJS.Signals) {
  try {
    process.kill(pid, signal);
  } catch (error) {
    if ((error as NodeJS.ErrnoException).code !== 'ESRCH') throw error;
  }
}

async function waitForExit(pid: number, timeoutMs: number): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    if (!processExists(pid)) return true;
    await new Promise((resolve) => setTimeout(resolve, 100));
  }
  return !processExists(pid);
}

export default async function globalTeardown() {
  const supervisorPath = join(runtimeRoot, 'supervisor.pid');
  if (!existsSync(supervisorPath)) return;
  const supervisorPid = Number(readFileSync(supervisorPath, 'utf8').trim());
  if (!Number.isSafeInteger(supervisorPid) || supervisorPid <= 1) {
    throw new Error('The local platform supervisor PID is invalid.');
  }
  signalProcess(supervisorPid, 'SIGTERM');
  if (await waitForExit(supervisorPid, 8_000)) return;

  const servicesPath = join(runtimeRoot, 'service-pids.json');
  if (!existsSync(servicesPath)) {
    throw new Error('Service PID evidence is missing during fallback cleanup.');
  }
  const cleanup = spawnSync(
    'python3',
    ['../scripts/run_local_platform.py', '--cleanup-service-pids', servicesPath],
    { cwd: process.cwd(), encoding: 'utf8', timeout: 15_000 },
  );
  if (cleanup.status !== 0) {
    throw new Error(
      `Session cleanup failed: ${cleanup.stderr || cleanup.stdout || 'unknown error'}`,
    );
  }
  signalProcess(supervisorPid, 'SIGKILL');
  if (!(await waitForExit(supervisorPid, 3_000))) {
    throw new Error('The local platform supervisor could not be reaped.');
  }
}
