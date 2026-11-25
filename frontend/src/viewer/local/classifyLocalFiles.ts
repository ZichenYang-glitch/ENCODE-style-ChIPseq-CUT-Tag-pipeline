import {
  PLAIN_PEAK_MAX_BYTES,
  READ_BIGWIG_SLICE,
  READ_TABIX_SLICE,
  READ_TEXT_FILE,
  assemblyKey,
  browserFileName,
  nativeTabixFiles,
  type PreparedLocalFile,
} from './contract'
import {
  hasBigWigIndexes,
  isBigBedMagic,
  isBigWigMagic,
  readFileSlice,
  scanPlainPeak,
} from './sliceRead'
import { tabixPairProblem, type IndexedLocalFormat } from './tabixPair'

export interface LocalSelectionIssue {
  id: string
  reason: string
  message: string
}

export interface LocalChoiceOption {
  id: string
  label: string
  prepared: PreparedLocalFile
  checkPair?: () => Promise<LocalSelectionIssue | null>
}

export interface LocalSelectionChoice {
  id: string
  reason: 'ambiguous-format' | 'ambiguous-index'
  message: string
  options: LocalChoiceOption[]
}

export interface ClassifiedLocalFiles {
  ready: PreparedLocalFile[]
  issues: LocalSelectionIssue[]
  choices: LocalSelectionChoice[]
}

type FileKind =
  | 'longrange'
  | 'bigwig'
  | 'plain-bed'
  | 'plain-narrow'
  | 'plain-broad'
  | 'indexed-bed'
  | 'indexed-narrow'
  | 'indexed-broad'
  | 'ambiguous-gz'
  | 'bedpe'
  | 'bare-longrange'
  | 'unknown'

interface NamedFile {
  file: File
  name: string
}

export async function classifyLocalFiles(
  input: Iterable<File>,
  assembly?: string,
): Promise<ClassifiedLocalFiles> {
  const named = Array.from(input).map((file) => ({
    file,
    name: browserFileName(file),
  }))
  const indexes: NamedFile[] = []
  const others: NamedFile[] = []
  for (const item of named) {
    if (item.name.toLowerCase().endsWith('.tbi')) indexes.push(item)
    else others.push(item)
  }
  const indexesByName = new Map<string, File[]>()
  for (const item of indexes) {
    const key = item.name.toLowerCase()
    const list = indexesByName.get(key) ?? []
    list.push(item.file)
    indexesByName.set(key, list)
  }
  const consumed = new Set<File>()
  const ready: PreparedLocalFile[] = []
  const issues: LocalSelectionIssue[] = []
  const choices: LocalSelectionChoice[] = []
  const indexedNames = new Map<string, number>()
  for (const item of others.filter(item => needsIndex(fileKind(item.name)))) {
    const key = `${item.name.toLowerCase()}.tbi`
    indexedNames.set(key, (indexedNames.get(key) ?? 0) + 1)
  }
  const unmatched = indexes
    .filter(item => !indexedNames.has(item.name.toLowerCase()))
    .map(item => item.file)

  for (const item of others) {
    const kind = fileKind(item.name)
    if (!needsIndex(kind)) {
      try {
        const prepared = await prepareUnindexed(item, kind, assembly)
        if ('issue' in prepared) issues.push(prepared.issue)
        else ready.push(prepared)
      } catch {
        issues.push(readError(item.file))
      }
      continue
    }
    const indexName = `${item.name.toLowerCase()}.tbi`
    const exact = indexesByName.get(indexName) ?? []
    const candidates = exact.length ? exact : unmatched
    for (const file of candidates) consumed.add(file)
    if (!candidates.length) {
      issues.push(issue('missing-index',
        `${item.name} 缺少 .tbi。请同时选择数据和对应索引，文件名不同时可以手动配对。`))
      continue
    }
    const formats = kind === 'ambiguous-gz'
      ? (['longrange', 'bed'] as const)
      : [fixedIndexedFormat(kind)]
    if (exact.length === 1 && formats.length === 1 && indexedNames.get(indexName) === 1) {
      const problem = await checkIndex(item.file, exact[0], formats[0])
      if (problem) issues.push(problem)
      else ready.push(optionFor(item, exact[0], formats[0]).prepared)
      continue
    }
    const options = candidates.flatMap(index => formats.map(format => ({
      ...optionFor(item, index, format),
      checkPair: () => checkIndex(item.file, index, format),
    })))
    choices.push({
      id: trackId(),
      reason: formats.length > 1 ? 'ambiguous-format' : 'ambiguous-index',
      message: `请为 ${item.name} 指定对应索引和格式。选定后会在浏览器中核对。`,
      options,
    })
  }

  for (const item of indexes) {
    if (consumed.has(item.file)) continue
    issues.push(
      issue(
        'unmatched-index',
        `${item.name} 没有同名的数据文件。请把它和对应的 .longrange.gz 或 .bed.gz 一起选择。`,
      ),
    )
  }
  return { ready, issues, choices }
}

function fileKind(name: string): FileKind {
  const lower = name.toLowerCase()
  if (lower.endsWith('.bedpe') || lower.endsWith('.bedpe.gz') || lower.endsWith('.bedpe.bgz')) {
    return 'bedpe'
  }
  if (lower.endsWith('.longrange.gz') || lower.endsWith('.longrange.bgz')) return 'longrange'
  if (lower.endsWith('.bw') || lower.endsWith('.bigwig')) return 'bigwig'
  if (lower.endsWith('.narrowpeak.gz') || lower.endsWith('.narrowpeak.bgz')) return 'indexed-narrow'
  if (lower.endsWith('.broadpeak.gz') || lower.endsWith('.broadpeak.bgz')) return 'indexed-broad'
  if (lower.endsWith('.bed.gz') || lower.endsWith('.bed.bgz')) return 'indexed-bed'
  if (lower.endsWith('.narrowpeak')) return 'plain-narrow'
  if (lower.endsWith('.broadpeak')) return 'plain-broad'
  if (lower.endsWith('.bed')) return 'plain-bed'
  if (lower.endsWith('.gz') || lower.endsWith('.bgz')) return 'ambiguous-gz'
  if (lower.includes('longrange')) return 'bare-longrange'
  return 'unknown'
}

function needsIndex(kind: FileKind): boolean {
  return (
    kind === 'longrange' ||
    kind === 'indexed-bed' ||
    kind === 'indexed-narrow' ||
    kind === 'indexed-broad' ||
    kind === 'ambiguous-gz'
  )
}

async function prepareUnindexed(
  item: NamedFile,
  kind: FileKind,
  assembly?: string,
): Promise<PreparedLocalFile | { issue: LocalSelectionIssue }> {
  if (kind === 'bigwig') return prepareBigWig(item)
  if (kind === 'plain-bed') return preparePlain(item, 'BED3', assembly)
  if (kind === 'plain-narrow') return preparePlain(item, 'narrowPeak', assembly)
  if (kind === 'plain-broad') return preparePlain(item, 'broadPeak', assembly)
  if (kind === 'bedpe') {
    return {
      issue: issue(
        'bedpe-not-converted',
        `${item.name} 是原始 BEDPE。本地轨道不把整份转成 longrange。请改用已经建好索引的 .longrange.gz 与 .tbi。`,
      ),
    }
  }
  if (kind === 'bare-longrange') {
    return {
      issue: issue(
        'missing-index',
        `${item.name} 不是带 .tbi 的 BGZF longrange。预建 longrange 不会整份读入或转换。`,
      ),
    }
  }
  return {
    issue: issue(
      'unsupported',
      `${item.name} 不是支持的本地轨道。支持预建 longrange（.longrange.gz + .tbi）、bigWig，以及 8 MiB 以内的 BED3/BED4、narrowPeak、broadPeak；索引峰文件用 .bed.gz + .tbi。`,
    ),
  }
}

async function prepareBigWig(
  item: NamedFile,
): Promise<PreparedLocalFile | { issue: LocalSelectionIssue }> {
  const header = await readFileSlice(item.file, 0, 64)
  if (isBigBedMagic(header)) {
    return {
      issue: issue(
        'bigbed-not-bigwig',
        `${item.name} 是 bigBed，不是 bigWig。请换一个 bigWig 文件。`,
      ),
    }
  }
  if (!isBigWigMagic(header)) {
    return {
      issue: issue(
        'not-bigwig',
        `${item.name} 的文件头不是 bigWig。请换一个 .bw 或 .bigWig。`,
      ),
    }
  }
  if (!(await hasBigWigIndexes(item.file, header))) {
    return {
      issue: issue(
        'invalid-bigwig',
        `${item.name} 缺少可读的 bigWig 文件头或信号索引。请重新选择完整的 .bw 或 .bigWig 文件。`,
      ),
    }
  }
  return {
    format: 'bigwig',
    washUType: 'bigwig',
    readMethod: READ_BIGWIG_SLICE,
    isText: false,
    label: item.name,
    dataFile: item.file,
  }
}

async function preparePlain(
  item: NamedFile,
  extensionFormat: 'BED3' | 'narrowPeak' | 'broadPeak',
  assembly?: string,
): Promise<PreparedLocalFile | { issue: LocalSelectionIssue }> {
  const scanned = await scanPlainPeak(item.file, extensionFormat)
  if ('error' in scanned) {
    const reason =
      item.file.size === 0
        ? 'empty-file'
        : item.file.size > PLAIN_PEAK_MAX_BYTES
          ? 'plain-too-large'
          : 'plain-invalid'
    return { issue: issue(reason, scanned.error) }
  }
  const declared = scanned.ok.declaredAssembly
  if (assembly && declared && assemblyKey(declared) !== assemblyKey(assembly)) {
    return { issue: issue('assembly-conflict',
      `${item.name} 的文件头声明 ${declared}，当前选择为 ${assembly}。请选择相符的组装或文件。`) }
  }
  return {
    format: scanned.ok.format,
    washUType: 'bed',
    readMethod: READ_TEXT_FILE,
    isText: true,
    label: item.name,
    dataFile: scanned.ok.dataFile,
  }
}

function readError(file: File): LocalSelectionIssue {
  return issue('unreadable-file', `浏览器无法读取 ${browserFileName(file)}。请确认文件可访问，再重新选择。`)
}

async function checkIndex(data: File, index: File, format: IndexedLocalFormat): Promise<LocalSelectionIssue | null> {
  try {
    const verdict = await tabixPairProblem(data, index, format)
    return verdict.ok ? null : issue(verdict.reason, verdict.message)
  } catch {
    return issue('unreadable-pair',
      `无法读取 ${browserFileName(data)} 与 ${browserFileName(index)}。请重新选择可访问的数据文件和索引。`)
  }
}

function fixedIndexedFormat(kind: FileKind): IndexedLocalFormat {
  if (kind === 'longrange') return 'longrange'
  if (kind === 'indexed-narrow') return 'narrowPeak'
  if (kind === 'indexed-broad') return 'broadPeak'
  return 'bed'
}

function optionFor(
  item: NamedFile,
  index: File,
  format: IndexedLocalFormat,
): LocalChoiceOption {
  const prepared: PreparedLocalFile = {
    format: format === 'longrange' ? 'longrange' : `indexed-${format}`,
    washUType: format === 'longrange' ? 'longrange' : 'bed',
    readMethod: READ_TABIX_SLICE,
    isText: false,
    label: item.name,
    dataFile: item.file,
    files: nativeTabixFiles(item.file, index),
  }
  const indexLabel = `${browserFileName(index)}（${index.size} 字节）`
  const action =
    format === 'longrange'
      ? '作为 longrange'
      : format === 'bed'
        ? '作为索引 BED'
        : `作为索引 ${format}`
  return {
    id: trackId(),
    label: `${item.name} + ${indexLabel}，${action}`,
    prepared,
  }
}

function issue(reason: string, message: string): LocalSelectionIssue {
  return { id: trackId(), reason, message }
}

function trackId(): string {
  if (typeof crypto !== 'undefined' && typeof crypto.randomUUID === 'function') {
    return crypto.randomUUID()
  }
  return `local-${Date.now().toString(16)}-${Math.random().toString(16).slice(2)}`
}
