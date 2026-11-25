/**
 * Pair a BGZF data file with a tabix index.
 *
 * A TBI has no checksum of the data file. Magic bytes, including a complete
 * index built for a different track, are not enough. The index is parsed in
 * full because it is the small file. The data file is only touched at the
 * header and at most two indexed virtual offsets, each limited to two BGZF
 * blocks.
 */

import { browserFileName } from './contract'
import { isBgzfHeader, isTabixMagic, readFileSlice } from './sliceRead'

const MAX_BGZF_BLOCK = 65536
const DATA_READ_CAP = 512 * 1024
const MAX_INDEX_UNCOMPRESSED = 32 * 1024 * 1024
const MAX_LINE_BYTES = 64 * 1024
const MAX_REFS = 1_000_000
const MAX_BINS = 100_000
const MAX_CHUNKS = 1_000_000
const MAX_INTERVALS = 1_000_000
const PSEUDO_BIN = 37450

export type IndexedLocalFormat = 'longrange' | 'bed' | 'narrowPeak' | 'broadPeak'
export type TabixPairReason = 'wrong-index' | 'not-bgzf' | 'index-mismatch' | 'invalid-longrange'

export type TabixPairVerdict =
  | { ok: true; dataBytesRead: number; maxDataOffset: number }
  | {
      ok: false
      reason: TabixPairReason
      message: string
      dataBytesRead: number
      maxDataOffset: number
    }

interface VirtualOffset {
  block: number
  data: number
}

interface IndexedChunk {
  refName: string
  bin: number
  start: VirtualOffset
  end: VirtualOffset
}

interface ParsedTabix {
  zeroBasedHalfOpen: boolean
  columns: { ref: number; start: number; end: number }
  chunks: IndexedChunk[]
  linearBlocks: number[]
}

type SliceReader = (start: number, length: number) => Promise<Uint8Array>

class DataReadCap extends Error {}
class IndexSizeLimit extends Error {}

const CRC_TABLE = (() => {
  const table = new Uint32Array(256)
  for (let index = 0; index < 256; index += 1) {
    let value = index
    for (let bit = 0; bit < 8; bit += 1) {
      value = (value & 1) === 1 ? 0xedb88320 ^ (value >>> 1) : value >>> 1
    }
    table[index] = value
  }
  return table
})()

export async function tabixPairProblem(
  data: File, index: File, format: IndexedLocalFormat = 'bed',
): Promise<TabixPairVerdict> {
  let dataBytesRead = 0
  let maxDataOffset = 0
  const readData: SliceReader = async (start, length) => {
    if (length > 0 && dataBytesRead + length > DATA_READ_CAP) throw new DataReadCap()
    const bytes = await readFileSlice(data, start, length)
    dataBytesRead += bytes.byteLength
    maxDataOffset = Math.max(maxDataOffset, start + bytes.byteLength)
    return bytes
  }
  try {
    const header = await readData(0, 18)
    if (!isBgzfHeader(header)) {
      return {
        ok: false,
        reason: 'not-bgzf',
        message: `${browserFileName(data)} 不是 BGZF。按窗口读取需要 BGZF 数据文件和 .tbi。`,
        dataBytesRead,
        maxDataOffset,
      }
    }
    if (isTabixMagic(await readFileSlice(index, 0, 4)) && index.size >= 36) {
      return wrongIndex(index, dataBytesRead, maxDataOffset, '是未压缩的 TBI。请选择 tabix 生成的 BGZF 压缩 .tbi。')
    }
    const payload = await readIndexPayload(index)
    const parsed = payload ? parseTabix(payload) : null
    if (!parsed || parsed.chunks.length === 0) {
      return wrongIndex(index, dataBytesRead, maxDataOffset)
    }
    if (!offsetsFit(parsed, data.size)) {
      return indexMismatch(data, index, dataBytesRead, maxDataOffset)
    }
    for (const chunk of sampleChunks(parsed.chunks)) {
      if (!await endsAtRecordBoundary(readData, data.size, chunk.end)) {
        return indexMismatch(data, index, dataBytesRead, maxDataOffset)
      }
      const record = await readIndexedRecord(readData, data.size, chunk.start, parsed.columns)
      if (!record || record.chrom !== chunk.refName) {
        return indexMismatch(data, index, dataBytesRead, maxDataOffset)
      }
      const bin = coordinateBin(record.start, record.end, parsed.zeroBasedHalfOpen)
      if (bin === null || !containsBin(chunk.bin, bin)) {
        return indexMismatch(data, index, dataBytesRead, maxDataOffset)
      }
      if (format === 'longrange' && !validLongrangeMate(record.fourthField)) {
        return {
          ok: false,
          reason: 'invalid-longrange',
          message: `${browserFileName(data)} 的第 4 列缺少可绘制的配对区间和分值。示例为 chr1:100-200,1。请选择预建 longrange，或将峰文件按索引 BED 打开。`,
          dataBytesRead,
          maxDataOffset,
        }
      }
    }
    return { ok: true, dataBytesRead, maxDataOffset }
  } catch (error) {
    if (error instanceof DataReadCap) {
      return indexMismatch(data, index, dataBytesRead, maxDataOffset)
    }
    if (error instanceof IndexSizeLimit) {
      return wrongIndex(index, dataBytesRead, maxDataOffset, '解压后超过 32 MiB。请按染色体拆分数据，并为每份数据重新生成匹配的 .tbi。')
    }
    return wrongIndex(index, dataBytesRead, maxDataOffset)
  }
}

function wrongIndex(
  index: File,
  dataBytesRead: number,
  maxDataOffset: number,
  detail = '不是完整的 tabix 索引。请重新选择数据自带的 .tbi。TBI 不含数据文件的整份摘要。',
): TabixPairVerdict {
  return {
    ok: false,
    reason: 'wrong-index',
    message: `${browserFileName(index)} ${detail}`,
    dataBytesRead,
    maxDataOffset,
  }
}

function indexMismatch(
  data: File,
  index: File,
  dataBytesRead: number,
  maxDataOffset: number,
): TabixPairVerdict {
  return {
    ok: false,
    reason: 'index-mismatch',
    message: `${browserFileName(index)} 的 tabix 结构完整，但和 ${browserFileName(data)} 对不上。染色体、区间或虚拟偏移不一致。TBI 不含整份数据摘要，请换用同一次生成的 .tbi。`,
    dataBytesRead,
    maxDataOffset,
  }
}

async function readIndexPayload(index: File): Promise<Uint8Array | null> {
  const readIndex: SliceReader = (start, length) => readFileSlice(index, start, length)
  const parts: Uint8Array[] = []
  let offset = 0
  let total = 0
  while (offset < index.size) {
    const block = await readBgzfBlock(readIndex, index.size, offset)
    if (!block) return null
    offset += block.size
    if (block.raw.length === 0) {
      if (offset !== index.size) return null
      break
    }
    total += block.raw.length
    if (total > MAX_INDEX_UNCOMPRESSED) throw new IndexSizeLimit()
    parts.push(block.raw)
  }
  if (parts.length === 0 || !isTabixMagic(parts[0] ?? new Uint8Array())) return null
  const payload = new Uint8Array(total)
  let cursor = 0
  for (const part of parts) {
    payload.set(part, cursor)
    cursor += part.byteLength
  }
  return payload
}

function parseTabix(bytes: Uint8Array): ParsedTabix | null {
  if (!isTabixMagic(bytes) || bytes.byteLength < 36) return null
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength)
  const referenceCount = view.getInt32(4, true)
  const formatFlags = view.getInt32(8, true)
  if (referenceCount <= 0 || referenceCount > MAX_REFS || (formatFlags & 0xf) > 2) return null
  const refColumn = view.getInt32(12, true)
  const startColumn = view.getInt32(16, true)
  const endColumn = view.getInt32(20, true)
  if (
    refColumn < 1 ||
    startColumn < 1 ||
    endColumn < 1 ||
    refColumn > 512 ||
    startColumn > 512 ||
    endColumn > 512
  ) {
    return null
  }
  const meta = view.getInt32(24, true)
  const skip = view.getInt32(28, true)
  const nameLength = view.getInt32(32, true)
  if (skip < 0 || meta < 0 || meta > 255 || nameLength < 2 || 36 + nameLength > bytes.byteLength) {
    return null
  }
  const names = splitNames(bytes.subarray(36, 36 + nameLength))
  if (!names || names.length !== referenceCount) return null
  let offset = 36 + nameLength
  const chunks: IndexedChunk[] = []
  const linearBlocks: number[] = []
  for (let reference = 0; reference < referenceCount; reference += 1) {
    if (offset + 4 > bytes.byteLength) return null
    const binCount = view.getUint32(offset, true)
    offset += 4
    if (binCount > MAX_BINS) return null
    const refName = names[reference]
    if (!refName) return null
    for (let binIndex = 0; binIndex < binCount; binIndex += 1) {
      if (offset + 8 > bytes.byteLength) return null
      const bin = view.getUint32(offset, true)
      const chunkCount = view.getUint32(offset + 4, true)
      offset += 8
      if (chunkCount > MAX_CHUNKS || offset + chunkCount * 16 > bytes.byteLength) return null
      if (bin > PSEUDO_BIN) return null
      for (let chunkIndex = 0; chunkIndex < chunkCount; chunkIndex += 1) {
        const start = readVirtual(view, offset)
        const end = readVirtual(view, offset + 8)
        offset += 16
        if (!start || !end) return null
        if (bin < PSEUDO_BIN) chunks.push({ refName, bin, start, end })
      }
    }
    if (offset + 4 > bytes.byteLength) return null
    const intervalCount = view.getUint32(offset, true)
    offset += 4
    if (intervalCount > MAX_INTERVALS || offset + intervalCount * 8 > bytes.byteLength) return null
    for (let interval = 0; interval < intervalCount; interval += 1) {
      const virtual = readVirtual(view, offset)
      offset += 8
      if (!virtual) return null
      if (virtual.block !== 0 || virtual.data !== 0) linearBlocks.push(virtual.block)
    }
  }
  // htslib appends the optional uint64 count of records without a coordinate.
  if (bytes.byteLength !== offset && bytes.byteLength !== offset + 8) return null
  return {
    zeroBasedHalfOpen: (formatFlags & 0x10000) !== 0,
    columns: { ref: refColumn, start: startColumn, end: endColumn },
    chunks,
    linearBlocks,
  }
}

function offsetsFit(parsed: ParsedTabix, fileSize: number): boolean {
  for (const chunk of parsed.chunks) {
    if (chunk.start.block >= fileSize || chunk.end.block > fileSize) return false
    if (chunk.end.block < chunk.start.block) return false
    if (chunk.end.block === chunk.start.block && chunk.end.data < chunk.start.data) return false
  }
  return parsed.linearBlocks.every((block) => block < fileSize)
}

function sampleChunks(chunks: IndexedChunk[]): IndexedChunk[] {
  const sorted = [...chunks].sort((left, right) => compareVirtual(left.start, right.start))
  const first = sorted[0]
  if (!first) return []
  for (let index = sorted.length - 1; index >= 1; index -= 1) {
    const candidate = sorted[index]
    if (!candidate) continue
    if (
      candidate.start.block !== first.start.block ||
      candidate.start.data !== first.start.data ||
      candidate.bin !== first.bin ||
      candidate.refName !== first.refName
    ) {
      return [first, candidate]
    }
  }
  return [first]
}

// Files can share a chromosome, first block and bin. The chunk must also
// end at a real BGZF block or immediately after a record in that block.
async function endsAtRecordBoundary(read: SliceReader, fileSize: number, end: VirtualOffset): Promise<boolean> {
  if (end.block === fileSize) return end.data === 0
  if (end.data === 0) return isBgzfHeader(await read(end.block, 18))
  const block = await readBgzfBlock(read, fileSize, end.block)
  if (!block || end.data > block.raw.length) return false
  if (block.raw[end.data - 1] === 0x0a) return true
  if (end.data !== block.raw.length) return false
  const next = end.block + block.size
  return next === fileSize || (await readBgzfBlock(read, fileSize, next))?.raw.length === 0
}

async function readIndexedRecord(
  read: SliceReader,
  fileSize: number,
  at: VirtualOffset,
  columns: ParsedTabix['columns'],
): Promise<{ chrom: string; start: number; end: number; fourthField?: string } | null> {
  const first = await readBgzfBlock(read, fileSize, at.block)
  if (!first || at.data > first.raw.byteLength) return null
  let raw = first.raw
  if (!hasNewline(raw.subarray(at.data))) {
    const next = await readBgzfBlock(read, fileSize, at.block + first.size)
    if (!next) return null
    const joined = new Uint8Array(raw.byteLength + next.raw.byteLength)
    joined.set(raw)
    joined.set(next.raw, raw.byteLength)
    raw = joined
  }
  const line = raw.subarray(at.data)
  const newline = line.indexOf(0x0a)
  if (newline <= 0 || newline > MAX_LINE_BYTES) return null
  let text: string
  try {
    text = new TextDecoder('utf-8', { fatal: true }).decode(line.subarray(0, newline))
  } catch {
    return null
  }
  if (text.endsWith('\r')) text = text.slice(0, -1)
  const fields = text.split('\t')
  const chrom = fields[columns.ref - 1]
  const start = parseCoordinate(fields[columns.start - 1])
  const end = parseCoordinate(fields[columns.end - 1])
  if (!chrom || start === null || end === null) return null
  return { chrom, start, end, fourthField: fields[3] }
}

function validLongrangeMate(field: string | undefined): boolean {
  if (!field) return false
  const match = /^[\w.]+:(\d+)-(\d+),([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)$/i.exec(field)
  if (!match) return false
  const start = Number(match[1])
  const end = Number(match[2])
  return Number.isSafeInteger(start) && Number.isSafeInteger(end) && start < end && Number.isFinite(Number(match[3]))
}

async function readBgzfBlock(
  read: SliceReader,
  fileSize: number,
  offset: number,
): Promise<{ raw: Uint8Array; size: number } | null> {
  if (offset < 0 || offset >= fileSize) return null
  const prefix = await read(offset, Math.min(64, fileSize - offset))
  if (!isBgzfHeader(prefix)) return null
  const extraLength = prefix[10]! | (prefix[11]! << 8)
  const extraEnd = 12 + extraLength
  if (extraEnd > MAX_BGZF_BLOCK) return null
  const extra = extraEnd <= prefix.byteLength ? prefix : await read(offset, extraEnd)
  if (extra.byteLength < extraEnd) return null
  let cursor = 12
  let blockSize = -1
  while (cursor + 4 <= extraEnd) {
    const subLength = extra[cursor + 2]! | (extra[cursor + 3]! << 8)
    if (extra[cursor] === 0x42 && extra[cursor + 1] === 0x43 && subLength >= 2 && cursor + 6 <= extra.byteLength) {
      blockSize = (extra[cursor + 4]! | (extra[cursor + 5]! << 8)) + 1
      break
    }
    cursor += 4 + subLength
  }
  if (blockSize < 28 || blockSize > MAX_BGZF_BLOCK || offset + blockSize > fileSize) return null
  const block = blockSize <= extra.byteLength ? extra.subarray(0, blockSize) : await read(offset, blockSize)
  if (block.byteLength < blockSize) return null
  let raw: Uint8Array
  try {
    raw = await inflateRaw(block.subarray(12 + extraLength, blockSize - 8))
  } catch {
    return null
  }
  const view = new DataView(block.buffer, block.byteOffset, block.byteLength)
  const crc = view.getUint32(blockSize - 8, true)
  const size = view.getUint32(blockSize - 4, true)
  if (crc !== crc32(raw) || size !== (raw.byteLength & 0xffffffff)) return null
  return { raw, size: blockSize }
}

async function inflateRaw(compressed: Uint8Array): Promise<Uint8Array> {
  if (typeof DecompressionStream !== 'function') {
    throw new Error('deflate-raw is unavailable')
  }
  const source = new Response(new Uint8Array(compressed)).body!
  return new Uint8Array(await new Response(
    source.pipeThrough(new DecompressionStream('deflate-raw')),
  ).arrayBuffer())
}

function splitNames(bytes: Uint8Array): string[] | null {
  const names: string[] = []
  let start = 0
  for (let index = 0; index < bytes.byteLength; index += 1) {
    if (bytes[index] !== 0) continue
    if (index === start) return null
    try {
      names.push(new TextDecoder('utf-8', { fatal: true }).decode(bytes.subarray(start, index)))
    } catch {
      return null
    }
    start = index + 1
  }
  if (start !== bytes.byteLength) return null
  return names
}

function readVirtual(view: DataView, offset: number): VirtualOffset | null {
  const low = view.getUint32(offset, true)
  const high = view.getUint32(offset + 4, true)
  const block = high * 0x10000 + (low >>> 16)
  if (!Number.isSafeInteger(block)) return null
  return { block, data: low & 0xffff }
}

function compareVirtual(left: VirtualOffset, right: VirtualOffset): number {
  if (left.block !== right.block) return left.block - right.block
  return left.data - right.data
}

// htslib can merge sparse child bins into their parent. The indexed
// chunk may therefore cover a broader interval than its first feature.
function containsBin(indexed: number, feature: number): boolean {
  while (feature > indexed) feature = (feature - 1) >>> 3
  return feature === indexed
}

function coordinateBin(start: number, end: number, zeroBasedHalfOpen: boolean): number | null {
  let begin = start
  let finish = end
  if (!zeroBasedHalfOpen) {
    if (start < 1) return null
    begin = start - 1
  }
  if (!Number.isSafeInteger(begin) || !Number.isSafeInteger(finish) || finish <= begin) return null
  return reg2bin(begin, finish)
}

function reg2bin(start: number, end: number): number {
  const last = end - 1
  if (start >> 14 === last >> 14) return Math.floor(((1 << 15) - 1) / 7) + (start >> 14)
  if (start >> 17 === last >> 17) return Math.floor(((1 << 12) - 1) / 7) + (start >> 17)
  if (start >> 20 === last >> 20) return Math.floor(((1 << 9) - 1) / 7) + (start >> 20)
  if (start >> 23 === last >> 23) return Math.floor(((1 << 6) - 1) / 7) + (start >> 23)
  if (start >> 26 === last >> 26) return Math.floor(((1 << 3) - 1) / 7) + (start >> 26)
  return 0
}

function parseCoordinate(value: string | undefined): number | null {
  if (!value || !/^(?:0|[1-9][0-9]{0,14})$/.test(value)) return null
  const parsed = Number(value)
  return Number.isSafeInteger(parsed) ? parsed : null
}

function hasNewline(bytes: Uint8Array): boolean {
  return bytes.indexOf(0x0a) !== -1
}

function crc32(bytes: Uint8Array): number {
  let crc = 0xffffffff
  for (const byte of bytes) crc = CRC_TABLE[(crc ^ byte) & 0xff]! ^ (crc >>> 8)
  return (crc ^ 0xffffffff) >>> 0
}
