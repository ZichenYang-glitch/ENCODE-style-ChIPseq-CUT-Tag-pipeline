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

export interface IndexedRecordOptions {
  lines: readonly string[]
  header?: string
  trailingNewline?: boolean
  endInDataBlock?: boolean
}

const recordCache = new Map<string, IndexedBed>()

export function indexedRecords(options: IndexedRecordOptions): IndexedBed {
  const key = JSON.stringify({
    lines: options.lines,
    header: options.header ?? '',
    trailingNewline: options.trailingNewline !== false,
    endInDataBlock: options.endInDataBlock === true,
  })
  const cached = recordCache.get(key)
  if (cached) return cached
  const script = String.raw`
import importlib.util, json, struct, sys
spec = importlib.util.spec_from_file_location("browser_bgzf", sys.argv[1])
mod = importlib.util.module_from_spec(spec)
sys.modules["browser_bgzf"] = mod
spec.loader.exec_module(mod)
request = json.load(sys.stdin)
header = request.get("header") or ""
lines = request["lines"]
trailing = bool(request.get("trailingNewline", True))
if header and not header.endswith("\n"):
    header += "\n"
body = bytearray(header.encode("ascii"))
spans = []
for index, line in enumerate(lines):
    if "\n" in line or "\r" in line:
        raise SystemExit("indexed line is invalid")
    parts = line.split("\t")
    if len(parts) < 3:
        raise SystemExit("indexed line is invalid")
    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
    if start < 0 or end <= start:
        raise SystemExit("indexed coordinate is invalid")
    text = line.encode("ascii") + (b"\n" if trailing or index < len(lines) - 1 else b"")
    record_start = len(body)
    body.extend(text)
    spans.append((chrom, start, end, record_start, len(body), text))
if not spans:
    raise SystemExit("indexed fixture was empty")
block_size = mod._MAX_UNCOMPRESSED
blocks = bytearray()
offsets = []
for pos in range(0, len(body), block_size):
    offsets.append(len(blocks))
    blocks.extend(mod._bgzf_block(body[pos:pos + block_size]))

def virtual(pos):
    if pos == len(body):
        if request["endInDataBlock"]:
            return (offsets[-1] << 16) | (pos - (len(offsets) - 1) * block_size)
        return len(blocks) << 16
    block, within = divmod(pos, block_size)
    return (offsets[block] << 16) | within

fixed = tuple(
    mod._LineSpan(chrom, start, end, virtual(begin), virtual(stop), text)
    for chrom, start, end, begin, stop, text in spans
)
tabix = mod._bgzf_compress(mod._write_tabix(fixed))
bgzf = blocks + mod._BGZF_EOF
raw = mod._decompress_bgzf(tabix)
sys.stdout.buffer.write(struct.pack("<III", len(bgzf), len(tabix), len(raw)))
sys.stdout.buffer.write(bgzf)
sys.stdout.buffer.write(tabix)
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
  recordCache.set(key, built)
  return built
}
