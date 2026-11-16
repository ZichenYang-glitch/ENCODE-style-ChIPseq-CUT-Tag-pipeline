import { spawn } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const root = new URL('../..', import.meta.url).pathname;
const head = await git(['rev-parse', 'HEAD']);
const out = process.env.WASHU_ACCEPT_DIR ?? join(tmpdir(), `helixweave-w1-closeout-${head.slice(0, 12)}`);
mkdirSync(out, { recursive: true });

const positive = await runProof(join(out, 'positive'), []);
const negative = await runProof(join(out, 'negative'), ['--negative']);
const index = {
  head,
  out,
  positive: { code: positive.code, tail: positive.tail },
  negative: { code: negative.code, tail: negative.tail },
  note: 'CI browser-e2e does not run this WashU proof.',
};
writeFileSync(join(out, 'ACCEPT.json'), `${JSON.stringify(index, null, 2)}\n`);
console.log(JSON.stringify(index, null, 2));
if (positive.code !== 0 || !positive.tail.includes('WASHU_PROOF ok')) process.exitCode = 1;
if (negative.code === 0 || negative.tail.includes('WASHU_PROOF ok')) process.exitCode = 1;
if (!negative.tail.includes('WASHU_PROOF rejected')) process.exitCode = 1;

function git(args) {
  return new Promise((resolve, reject) => {
    const child = spawn('git', args, { cwd: root });
    let text = '';
    child.stdout.on('data', (chunk) => {
      text += chunk.toString();
    });
    child.on('exit', (code) => {
      if (code === 0) resolve(text.trim());
      else reject(new Error(`git ${args.join(' ')} exited ${code}`));
    });
  });
}

function runProof(dir, args) {
  return new Promise((resolve) => {
    mkdirSync(dir, { recursive: true });
    const child = spawn(process.execPath, ['frontend/scripts/washu-browser-proof.mjs', ...args], {
      cwd: root,
      env: { ...process.env, WASHU_PROOF_DIR: dir },
    });
    let text = '';
    child.stdout.on('data', (chunk) => {
      text += chunk.toString();
      process.stdout.write(chunk);
    });
    child.stderr.on('data', (chunk) => {
      text += chunk.toString();
      process.stderr.write(chunk);
    });
    child.on('exit', (code) => {
      resolve({ code: code ?? 1, tail: text.trim().split('\n').slice(-4).join('\n') });
    });
  });
}
