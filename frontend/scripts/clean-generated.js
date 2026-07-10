import { readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs';
import { resolve, join } from 'node:path';

const generatedDir = resolve(import.meta.dirname, '../src/api/generated');

function cleanFile(filePath) {
  const content = readFileSync(filePath, 'utf-8');
  const cleaned = content
    .split('\n')
    .map((line) => line.replace(/\s+$/, ''))
    .join('\n')
    .replace(/\n+$/, '\n');

  if (content !== cleaned) {
    writeFileSync(filePath, cleaned, 'utf-8');
  }
}

function walk(dir) {
  for (const entry of readdirSync(dir)) {
    const fullPath = join(dir, entry);
    const info = statSync(fullPath);
    if (info.isDirectory()) {
      walk(fullPath);
    } else if (fullPath.endsWith('.ts')) {
      cleanFile(fullPath);
    }
  }
}

walk(generatedDir);
