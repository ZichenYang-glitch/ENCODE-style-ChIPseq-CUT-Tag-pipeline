import { PLAIN_PEAK_MAX_BYTES, assemblyKey, browserFileName } from './contract'

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

export async function hasBigWigIndexes(file: File, header: Uint8Array): Promise<boolean> {
  if (header.length < 64) return false
  const little = header[0] === 0x26
  const view = new DataView(header.buffer, header.byteOffset, header.byteLength)
  const headerEnd = 64 + 24 * view.getUint16(6, little)
  const chromTree = Number(view.getBigUint64(8, little))
  const data = Number(view.getBigUint64(16, little))
  const index = Number(view.getBigUint64(24, little))
  if (
    chromTree < headerEnd || chromTree + 32 > file.size ||
    data < headerEnd || data + 8 > file.size ||
    index < headerEnd || index + 48 > file.size
  ) return false
  const [chromMagic, indexMagic] = await Promise.all([
    readFileSlice(file, chromTree, 4),
    readFileSlice(file, index, 4),
  ])
  return new DataView(chromMagic.buffer, chromMagic.byteOffset, 4).getUint32(0, little) === 0x78ca8c91 &&
    new DataView(indexMagic.buffer, indexMagic.byteOffset, 4).getUint32(0, little) === 0x2468ace0
}

export function isPeakDirective(raw: string): boolean {
  const line = raw.endsWith('\r') ? raw.slice(0, -1) : raw
  if (line === '') return true
  const head = line.split('\t', 1)[0] ?? ''
  return (
    line.startsWith('#') ||
    head === 'track' ||
    head.startsWith('track ') ||
    head === 'browser' ||
    head.startsWith('browser ')
  )
}

// UCSC track and browser directives can explicitly name an assembly.
export function directiveAssembly(raw: string): string | undefined {
  const line = raw.endsWith('\r') ? raw.slice(0, -1) : raw
  if (!/^track\s/i.test(line)) {
    return /^browser\s+genome\s+([A-Za-z0-9_.-]+)/i.exec(line)?.[1]
  }
  // Consume each complete value so a quoted label can contain text such as "with db=hg38".
  const attributes = /\s([^\s=]+)=("(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\S*)/g
  for (const [, key, value] of line.matchAll(attributes)) {
    if (/^(?:db|genome)$/i.test(key)) {
      const assembly = /^["']?([A-Za-z0-9_.-]+)/.exec(value)
      if (assembly) return assembly[1]
    }
  }
  return undefined
}

export interface PlainPeakScan {
  format: 'BED3' | 'BED4' | 'narrowPeak' | 'broadPeak'
  dataRows: number
  declaredAssembly?: string
  dataFile: File
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
  const dataParts: string[] = []
  let rows: string[] = []
  let skippedLines = false
  const finishLine = (raw: string): { error: string } | null => {
    const verdict = inspectPeakLine(raw, lineNumber, file.name)
    lineNumber += 1
    if (typeof verdict === 'object') return verdict
    if (verdict === 'skip') {
      skippedLines = true
      const value = directiveAssembly(raw)
      if (value && declaredAssembly && assemblyKey(value) !== assemblyKey(declaredAssembly)) {
        return { error: `${file.name} 的文件头分别声明 ${declaredAssembly} 和 ${value}。请核对参考组装。` }
      }
      if (value) declaredAssembly = value
    }
    if (verdict === 'data') {
      dataRows += 1
      rows.push(raw)
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
    if (rows.length) dataParts.push(rows.join('\n') + '\n')
    rows = []
  }
  if (pending) {
    const problem = finishLine(pending)
    if (problem) return problem
    dataParts.push(...rows)
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
  // WashU infers the delimiter from the first rows. UCSC directives stay in the
  // assembly check above and the native reader receives just the validated intervals.
  const dataFile = skippedLines
    ? new File(dataParts, browserFileName(file), { type: file.type, lastModified: file.lastModified })
    : file
  return { ok: { format, dataRows, declaredAssembly, dataFile } }
}

function inspectPeakLine(
  raw: string,
  lineNumber: number,
  fileName: string,
): 'skip' | 'data' | { error: string } {
  const line = raw.endsWith('\r') ? raw.slice(0, -1) : raw
  if (isPeakDirective(line)) return 'skip'
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
