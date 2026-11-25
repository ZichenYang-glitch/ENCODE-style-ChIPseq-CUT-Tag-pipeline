import { execFileSync, spawnSync } from 'node:child_process'
import { readFileSync } from 'node:fs'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../../../..')
const cache = new Map<string, IndexedBed>()

export function signalBigWig(): Uint8Array {
  return readFileSync(path.join(repoRoot, 'test/fixtures/browser/chr1-signal.bw'))
}

export function bgzfIndex(payload: Uint8Array): Uint8Array {
  const script = [
    'import sys',
    'from browser_bgzf import _bgzf_compress',
    'sys.stdout.buffer.write(_bgzf_compress(sys.stdin.buffer.read()))',
  ].join('\n')
  return execFileSync('python3', ['-c', script], {
    cwd: path.join(repoRoot, 'src/encode_pipeline/services'),
    input: payload,
  })
}

export interface IndexedBed {
  bgzf: Uint8Array
  tabix: Uint8Array
  rawTabix: Uint8Array
}

export function indexedPair(lines: readonly string[]): IndexedBed {
  const key = JSON.stringify(lines)
  const cached = cache.get(key)
  if (cached) return cached
  const script = `
import importlib.util, json, struct, sys
spec = importlib.util.spec_from_file_location("browser_bgzf", sys.argv[1])
mod = importlib.util.module_from_spec(spec)
sys.modules["browser_bgzf"] = mod
spec.loader.exec_module(mod)
indexed = mod.write_indexed_bed(tuple(json.load(sys.stdin)))
raw = mod._decompress_bgzf(indexed.tabix)
sys.stdout.buffer.write(struct.pack("<III", len(indexed.bgzf), len(indexed.tabix), len(raw)))
sys.stdout.buffer.write(indexed.bgzf)
sys.stdout.buffer.write(indexed.tabix)
sys.stdout.buffer.write(raw)
`
  const result = spawnSync(
    'python3',
    ['-c', script, path.join(repoRoot, 'src/encode_pipeline/services/browser_bgzf.py')],
    { input: key, maxBuffer: 32 * 1024 * 1024 },
  )
  if (result.status !== 0) {
    const detail = result.stderr?.toString() || result.stdout?.toString() || 'could not build indexed fixture'
    throw new Error(detail)
  }
  const buffer = result.stdout
  if (!Buffer.isBuffer(buffer) || buffer.length < 12) {
    throw new Error('indexed fixture was empty')
  }
  const bgzfLength = buffer.readUInt32LE(0)
  const tabixLength = buffer.readUInt32LE(4)
  const rawLength = buffer.readUInt32LE(8)
  const bgzfEnd = 12 + bgzfLength
  const tabixEnd = bgzfEnd + tabixLength
  const rawEnd = tabixEnd + rawLength
  if (rawEnd !== buffer.length) throw new Error('indexed fixture lengths did not match')
  const built = {
    bgzf: new Uint8Array(buffer.subarray(12, bgzfEnd)),
    tabix: new Uint8Array(buffer.subarray(bgzfEnd, tabixEnd)),
    rawTabix: new Uint8Array(buffer.subarray(tabixEnd, rawEnd)),
  }
  cache.set(key, built)
  return built
}
