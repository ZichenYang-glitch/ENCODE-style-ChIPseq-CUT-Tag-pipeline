import { existsSync, readFileSync } from 'node:fs';
import { join } from 'node:path';

const runtimeRoot =
  process.env.ENCODE_PIPELINE_E2E_ROOT ?? '/tmp/encode-platform-playwright';

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

function signalGroup(pid: number, signal: NodeJS.Signals) {
  try {
    process.kill(-pid, signal);
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
  if (existsSync(servicesPath)) {
    const services = JSON.parse(readFileSync(servicesPath, 'utf8')) as Record<
      string,
      number
    >;
    for (const pid of Object.values(services)) signalGroup(pid, 'SIGTERM');
    await new Promise((resolve) => setTimeout(resolve, 1_000));
    for (const pid of Object.values(services)) signalGroup(pid, 'SIGKILL');
  }
  signalProcess(supervisorPid, 'SIGKILL');
  if (!(await waitForExit(supervisorPid, 3_000))) {
    throw new Error('The local platform supervisor could not be reaped.');
  }
}
