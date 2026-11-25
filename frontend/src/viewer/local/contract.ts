/**
 * Models for GenomeHub in wuepgg 62.0.2.
 *
 * wuepgg3-track 9.0.3 takes the local reader only when `fileObj` is a File
 * and `url` is `''`. Any other URL is the remote reader, and that reader can
 * fall through to a public proxy. Indexed longrange and BED pass `files` as
 * [data, index]. The track package treats the longer name as the .tbi.
 * bigWig passes the File on `fileObj` alone. Plain peaks set `isText` after
 * a byte cap, because that reader parses the whole file.
 */

export const PLAIN_PEAK_MAX_BYTES = 8 * 1024 * 1024

export const READ_TABIX_SLICE = 'wuepgg3-track-tabix-file-slice'
export const READ_BIGWIG_SLICE = 'wuepgg3-track-bigwig-file-slice'
export const READ_TEXT_FILE = 'wuepgg3-track-text-file'

export interface LocalWashuTrackModel {
  type: string
  name: string
  url: ''
  id: string
  fileObj: File
  files?: File[]
  isText: boolean
  metadata: {
    genome: string
    assemblySignature: string
    readMethod: string
    format: string
    localOnly: true
  }
  [key: string]: unknown
}

export interface LocalWashuTrack {
  id: string
  label: string
  assembly: string
  model: LocalWashuTrackModel
}

export interface PreparedLocalFile {
  format: string
  washUType: 'longrange' | 'bigwig' | 'bed'
  readMethod: string
  isText: boolean
  label: string
  dataFile: File
  files?: File[]
}

export function browserFileName(file: File): string {
  const leaf = file.name.split(/[/\\]/).pop() ?? ''
  const trimmed = leaf.trim()
  return trimmed === '' ? 'local-track' : trimmed
}

// Equivalent public names of the same browser coordinate scaffold.
export function assemblyKey(value: string): string {
  const key = value.trim().toLowerCase()
  if (/^(?:hg38|grch38)(?:\.p[0-9]+)?$/.test(key)) return 'hg38'
  if (/^(?:hg19|grch37)(?:\.p[0-9]+)?$/.test(key)) return 'hg19'
  if (key === 'grcm38' || key === 'grcm38.p6') return 'mm10'
  if (key === 'grcm39') return 'mm39'
  return key
}

export function assemblySignature(assembly: string): string {
  const key = assemblyKey(assembly)
  if (key === 'hg38') return '用户指定 GRCh38'
  if (key === 'mm10') return '用户指定 GRCm38'
  return `用户指定 ${assembly.trim()}`
}

export function nativeTabixFiles(data: File, index: File): File[] {
  if (index.name.length >= data.name.length) return [data, index]
  const renamed = new File([index], `${data.name}.tbi`, {
    type: index.type,
    lastModified: index.lastModified,
  })
  return [data, renamed]
}

export function buildLocalTrack(
  prepared: PreparedLocalFile,
  assembly: string,
): LocalWashuTrack {
  const id = newTrackId()
  const signature = assemblySignature(assembly)
  const model: LocalWashuTrackModel = {
    type: prepared.washUType,
    name: prepared.label,
    url: '',
    id,
    fileObj: prepared.dataFile,
    isText: prepared.isText,
    ...(prepared.isText && { textConfig: {} }),
    metadata: {
      genome: assembly,
      assemblySignature: signature,
      readMethod: prepared.readMethod,
      format: prepared.format,
      localOnly: true,
    },
  }
  if (prepared.files) model.files = prepared.files
  return {
    id,
    label: prepared.label,
    assembly,
    model,
  }
}

function newTrackId(): string {
  if (typeof crypto !== 'undefined' && typeof crypto.randomUUID === 'function') {
    return crypto.randomUUID()
  }
  return `local-${Date.now().toString(16)}-${Math.random().toString(16).slice(2)}`
}
