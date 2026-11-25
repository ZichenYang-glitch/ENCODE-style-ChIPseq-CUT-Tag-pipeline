import { PLAIN_PEAK_MAX_BYTES, assemblyKey } from './contract'

export const HEADER_BYTES = 18
export const SCAN_WINDOW_BYTES = 64 * 1024
const MAX_LINE_BYTES = 64 * 1024
const CHROM = /^[A-Za-z0-9_.:-]{1,128}$/
const INTEGER = /^(?:0|[1-9][0-9]{0,14})$/

export async function readFileSlice(
  file: Blob,
  start: number,
  length: number,
): Promise<Uint8Array> {
  if (length <= 0 || start < 0 || start >= file.size) return new Uint8Array()
  const end = Math.min(file.size, start + length)
  return bytesOf(file.slice(start, end))
}

async function bytesOf(blob: Blob): Promise<Uint8Array> {
  if (typeof blob.arrayBuffer === 'function') {
    return new Uint8Array(await blob.arrayBuffer())
  }
  // jsdom's Blob has slice, but the sliced object has no arrayBuffer. The
  // track package reads that slice with FileReader, which is available there.
  const reader = new FileReader()
  const buffer = await new Promise<ArrayBuffer>((resolve, reject) => {
    reader.onerror = () => {
      reader.abort()
      reject(reader.error ?? new Error('slice read failed'))
    }
    reader.onload = () => {
      const result = reader.result
      if (result instanceof ArrayBuffer) resolve(result)
      else reject(new Error('slice read failed'))
    }
    reader.readAsArrayBuffer(blob)
  })
  return new Uint8Array(buffer)
}

export function isBgzfHeader(bytes: Uint8Array): boolean {
  if (bytes.length < 16) return false
  if (bytes[0] !== 0x1f || bytes[1] !== 0x8b || bytes[2] !== 8) return false
  if ((bytes[3] & 0x04) === 0) return false
  return bytes[12] === 0x42 && bytes[13] === 0x43
}

export function isTabixMagic(bytes: Uint8Array): boolean {
  return (
    bytes.length >= 4 &&
    bytes[0] === 0x54 &&
    bytes[1] === 0x42 &&
    bytes[2] === 0x49 &&
    bytes[3] === 0x01
  )
}

/** A .tbi is either raw TBI bytes or a BGZF block whose payload starts with TBI. */
export async function looksLikeTabix(file: Blob): Promise<boolean> {
  const head = await readFileSlice(file, 0, HEADER_BYTES)
  if (isTabixMagic(head)) return true
  if (!isBgzfHeader(head)) return false
  if (typeof DecompressionStream !== 'function') return true
  try {
    const inflated = await inflateFirstBgzfBlock(file)
    return inflated !== null && isTabixMagic(inflated)
  } catch {
    return false
  }
}

async function inflateFirstBgzfBlock(file: Blob): Promise<Uint8Array | null> {
  const prefix = await readFileSlice(file, 0, 64)
  if (!isBgzfHeader(prefix) || prefix.length < 18) return null
  const xlen = prefix[10]! | (prefix[11]! << 8)
  const extraEnd = 12 + xlen
  const extra = extraEnd <= prefix.length ? prefix : await readFileSlice(file, 0, extraEnd)
  let cursor = 12
  let blockSize = -1
  while (cursor + 4 <= extraEnd && cursor + 4 <= extra.length) {
    const si1 = extra[cursor]
    const si2 = extra[cursor + 1]
    const slen = extra[cursor + 2]! | (extra[cursor + 3]! << 8)
    if (si1 === 0x42 && si2 === 0x43 && slen >= 2 && cursor + 6 <= extra.length) {
      blockSize = (extra[cursor + 4]! | (extra[cursor + 5]! << 8)) + 1
      break
    }
    cursor += 4 + slen
  }
  if (blockSize < 28 || blockSize > file.size) return null
  const block = await readFileSlice(file, 0, blockSize)
  if (block.length < blockSize) return null
  const compressed = block.slice(12 + xlen, block.length - 8)
  const stream = new DecompressionStream('deflate-raw')
  const writer = stream.writable.getWriter()
  await writer.write(compressed)
  await writer.close()
  return new Uint8Array(await new Response(stream.readable).arrayBuffer())
}

export function isBigWigMagic(bytes: Uint8Array): boolean {
  if (bytes.length < 4) return false
  const little =
    bytes[0] === 0x26 && bytes[1] === 0xfc && bytes[2] === 0x8f && bytes[3] === 0x88
  const big =
    bytes[0] === 0x88 && bytes[1] === 0x8f && bytes[2] === 0xfc && bytes[3] === 0x26
  return little || big
}

export function isBigBedMagic(bytes: Uint8Array): boolean {
  if (bytes.length < 4) return false
  const little =
    bytes[0] === 0xeb && bytes[1] === 0xf2 && bytes[2] === 0x89 && bytes[3] === 0x87
  const big =
    bytes[0] === 0x87 && bytes[1] === 0x89 && bytes[2] === 0xf2 && bytes[3] === 0xeb
  return little || big
}

export interface PlainPeakScan {
  format: 'BED3' | 'BED4' | 'narrowPeak' | 'broadPeak'
  dataRows: number
  declaredAssembly?: string
}

export async function scanPlainPeak(
  file: File,
  extensionFormat: PlainPeakScan['format'],
): Promise<{ ok: PlainPeakScan } | { error: string }> {
  if (file.size > PLAIN_PEAK_MAX_BYTES) {
    return {
      error: `${file.name} 有 ${file.size} 字节，超过纯文本峰文件上限 ${PLAIN_PEAK_MAX_BYTES} 字节（8 MiB）。请改用带 .tbi 的 BGZF 峰文件，按窗口读取。这个文件没有被读入。`,
    }
  }
  if (file.size === 0) {
    return { error: `${file.name} 是空文件。` }
  }
  const decoder = new TextDecoder('utf-8', { fatal: true })
  let offset = 0
  let pending = ''
  let lineNumber = 1
  let dataRows = 0
  let maxColumns = 0
  let declaredAssembly: string | undefined
  const finishLine = (raw: string): { error: string } | null => {
    const verdict = inspectPeakLine(raw, lineNumber, file.name)
    lineNumber += 1
    if (typeof verdict === 'object') return verdict
    if (verdict === 'skip') {
      // UCSC track and browser directives can explicitly name an assembly.
      const declaration = /^track\s/i.test(raw)
        ? /(?:^|\s)(?:db|genome)=["']?([A-Za-z0-9_.-]+)/i.exec(raw)
        : /^browser\s+genome\s+([A-Za-z0-9_.-]+)/i.exec(raw)
      const value = declaration?.[1]
      if (value && declaredAssembly && assemblyKey(value) !== assemblyKey(declaredAssembly)) {
        return { error: `${file.name} 的文件头分别声明 ${declaredAssembly} 和 ${value}。请核对参考组装。` }
      }
      if (value) declaredAssembly = value
    }
    if (verdict === 'data') {
      dataRows += 1
      maxColumns = Math.max(maxColumns, raw.replace(/\r$/, '').split('\t').length)
    }
    return null
  }
  while (offset < file.size) {
    const bytes = await readFileSlice(file, offset, SCAN_WINDOW_BYTES)
    if (bytes.byteLength === 0) break
    offset += bytes.byteLength
    if (bytes.includes(0)) {
      return { error: `${file.name} 不是纯文本 BED。` }
    }
    let text: string
    try {
      text = decoder.decode(bytes, { stream: offset < file.size })
    } catch {
      return { error: `${file.name} 不是 UTF-8 文本。` }
    }
    pending += text
    if (pending.length > MAX_LINE_BYTES && !pending.includes('\n')) {
      return { error: `${file.name} 第 ${lineNumber} 行超过 ${MAX_LINE_BYTES} 字节。` }
    }
    const parts = pending.split('\n')
    pending = parts.pop() ?? ''
    if (pending.length > MAX_LINE_BYTES) {
      return { error: `${file.name} 第 ${lineNumber} 行超过 ${MAX_LINE_BYTES} 字节。` }
    }
    for (const part of parts) {
      const problem = finishLine(part)
      if (problem) return problem
    }
  }
  if (pending) {
    const problem = finishLine(pending)
    if (problem) return problem
  }
  if (dataRows === 0) {
    return { error: `${file.name} 没有可绘制的区间。` }
  }
  const format =
    extensionFormat === 'BED3' || extensionFormat === 'BED4'
      ? maxColumns === 3
        ? 'BED3'
        : 'BED4'
      : extensionFormat
  return { ok: { format, dataRows, declaredAssembly } }
}

function inspectPeakLine(
  raw: string,
  lineNumber: number,
  fileName: string,
): 'skip' | 'data' | { error: string } {
  const line = raw.endsWith('\r') ? raw.slice(0, -1) : raw
  if (line === '') return 'skip'
  const head = line.split('\t', 1)[0] ?? ''
  if (
    line.startsWith('#') ||
    head === 'track' ||
    head.startsWith('track ') ||
    head === 'browser' ||
    head.startsWith('browser ')
  ) {
    return 'skip'
  }
  const columns = line.split('\t')
  const invalid = `${fileName} 第 ${lineNumber} 行不是支持的 BED3/BED4 峰行。需要制表符分隔、至少三列，以及半开整数坐标。`
  if (columns.length < 3) return { error: invalid }
  if (!CHROM.test(columns[0] ?? '')) return { error: invalid }
  if (!INTEGER.test(columns[1] ?? '') || !INTEGER.test(columns[2] ?? '')) {
    return { error: invalid }
  }
  const start = Number(columns[1])
  const end = Number(columns[2])
  if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || end <= start) {
    return { error: invalid }
  }
  return 'data'
}
