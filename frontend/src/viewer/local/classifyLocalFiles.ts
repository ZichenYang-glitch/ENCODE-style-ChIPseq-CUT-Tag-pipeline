import {
  PLAIN_PEAK_MAX_BYTES,
  READ_BIGWIG_SLICE,
  READ_TABIX_SLICE,
  READ_TEXT_FILE,
  browserFileName,
  nativeTabixFiles,
  type PreparedLocalFile,
} from './contract'
import {
  HEADER_BYTES,
  isBgzfHeader,
  isBigBedMagic,
  isBigWigMagic,
  looksLikeTabix,
  readFileSlice,
  scanPlainPeak,
} from './sliceRead'

export interface LocalSelectionIssue {
  id: string
  reason: string
  message: string
}

export interface LocalChoiceOption {
  id: string
  label: string
  prepared: PreparedLocalFile
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
  const seenIndexKeys = new Set<string>()

  for (const item of others) {
    const kind = fileKind(item.name)
    if (!needsIndex(kind)) {
      const prepared = await prepareUnindexed(item, kind)
      if ('issue' in prepared) issues.push(prepared.issue)
      else ready.push(prepared)
      continue
    }
    const key = `${item.name.toLowerCase()}.tbi`
    if (seenIndexKeys.has(key)) continue
    seenIndexKeys.add(key)
    const group = others.filter(
      (candidate) =>
        needsIndex(fileKind(candidate.name)) &&
        `${candidate.name.toLowerCase()}.tbi` === key,
    )
    const indexFiles = indexesByName.get(key) ?? []
    for (const indexFile of indexFiles) consumed.add(indexFile)
    await resolveIndexedGroup(group, indexFiles, ready, issues, choices)
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
): Promise<PreparedLocalFile | { issue: LocalSelectionIssue }> {
  if (kind === 'bigwig') return prepareBigWig(item)
  if (kind === 'plain-bed') return preparePlain(item, 'BED3')
  if (kind === 'plain-narrow') return preparePlain(item, 'narrowPeak')
  if (kind === 'plain-broad') return preparePlain(item, 'broadPeak')
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
  const header = await readFileSlice(item.file, 0, 4)
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
  return {
    format: scanned.ok.format,
    washUType: 'bed',
    readMethod: READ_TEXT_FILE,
    isText: true,
    label: item.name,
    dataFile: item.file,
  }
}

async function resolveIndexedGroup(
  group: NamedFile[],
  indexFiles: File[],
  ready: PreparedLocalFile[],
  issues: LocalSelectionIssue[],
  choices: LocalSelectionChoice[],
): Promise<void> {
  if (indexFiles.length === 0) {
    for (const item of group) {
      issues.push(
        issue(
          'missing-index',
          `${item.name} 缺少同名 .tbi。预建索引轨道不会整份解压或转换。请同时选择索引，或移除该文件。`,
        ),
      )
    }
    return
  }
  const viable: LocalChoiceOption[] = []
  const failures: LocalSelectionIssue[] = []
  for (const item of group) {
    for (const indexFile of indexFiles) {
      const problem = await indexedProblem(item.file, indexFile)
      if (problem) {
        failures.push(problem)
        continue
      }
      const kind = fileKind(item.name)
      if (kind === 'ambiguous-gz') {
        viable.push(optionFor(item, indexFile, 'longrange'))
        viable.push(optionFor(item, indexFile, 'bed'))
      } else {
        viable.push(optionFor(item, indexFile, fixedIndexedFormat(kind)))
      }
    }
  }
  if (viable.length === 0) {
    issues.push(...failures)
    return
  }
  if (failures.length > 0) issues.push(...failures)
  if (viable.length === 1) {
    const only = viable[0]
    if (only) ready.push(only.prepared)
    return
  }
  const formats = new Set(viable.map((option) => option.prepared.washUType))
  choices.push({
    id: trackId(),
    reason: formats.size > 1 ? 'ambiguous-format' : 'ambiguous-index',
    message:
      formats.size > 1
        ? '这些文件有不止一种读法。请选择要添加的那一种；其余不会读取。'
        : '有多个索引候选。请选择一个；未选择则不读取。',
    options: viable,
  })
}

async function indexedProblem(data: File, index: File): Promise<LocalSelectionIssue | null> {
  const [dataHeader, indexOk] = await Promise.all([
    readFileSlice(data, 0, HEADER_BYTES),
    looksLikeTabix(index),
  ])
  if (!indexOk) {
    return issue(
      'wrong-index',
      `${browserFileName(index)} 不是 tabix 索引。需要未压缩的 TBI 头，或 BGZF 里以 TBI 开头的 .tbi。`,
    )
  }
  if (!isBgzfHeader(dataHeader)) {
    return issue(
      'not-bgzf',
      `${browserFileName(data)} 不是 BGZF。按窗口读取需要 BGZF 数据文件和 .tbi。`,
    )
  }
  return null
}

function fixedIndexedFormat(kind: FileKind): 'longrange' | 'bed' | 'narrowPeak' | 'broadPeak' {
  if (kind === 'longrange') return 'longrange'
  if (kind === 'indexed-narrow') return 'narrowPeak'
  if (kind === 'indexed-broad') return 'broadPeak'
  return 'bed'
}

function optionFor(
  item: NamedFile,
  index: File,
  format: 'longrange' | 'bed' | 'narrowPeak' | 'broadPeak',
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
