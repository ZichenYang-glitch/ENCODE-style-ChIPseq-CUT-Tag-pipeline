import { afterEach, describe, expect, it, vi } from 'vitest'
import {
  PLAIN_PEAK_MAX_BYTES,
  READ_BIGWIG_SLICE,
  READ_TABIX_SLICE,
  READ_TEXT_FILE,
  assemblySignature,
  nativeTabixFiles,
} from './contract'
import { classifyLocalFiles } from './classifyLocalFiles'

function bgzfHeader(): Uint8Array {
  return new Uint8Array([
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00,
  ])
}

function tbiBytes(): Uint8Array {
  return new Uint8Array([0x54, 0x42, 0x49, 0x01, 0, 0, 0, 0])
}

function fileBytes(bytes: Uint8Array): ArrayBuffer {
  const copy = new ArrayBuffer(bytes.byteLength)
  new Uint8Array(copy).set(bytes)
  return copy
}

function fileNamed(name: string, bytes: Uint8Array): File {
  return new File([fileBytes(bytes)], name)
}

afterEach(() => {
  vi.restoreAllMocks()
})

describe('assemblySignature', () => {
  it('signs hg38 and GRCh38 as the user-declared reference', () => {
    expect(assemblySignature('hg38')).toBe('用户指定 GRCh38')
    expect(assemblySignature('GRCh38')).toBe('用户指定 GRCh38')
    expect(assemblySignature('GRCh38.p14')).toBe('用户指定 GRCh38')
    expect(assemblySignature('mm10')).toBe('用户指定 mm10')
  })
})

describe('nativeTabixFiles', () => {
  it('keeps the original data file and makes the index name the longer one', () => {
    const data = fileNamed('abcdefghij.longrange.gz', bgzfHeader())
    const index = fileNamed('x.tbi', tbiBytes())
    const [kept, renamed] = nativeTabixFiles(data, index)
    expect(kept).toBe(data)
    expect(renamed.name).toBe('abcdefghij.longrange.gz.tbi')
    expect(renamed.name.length).toBeGreaterThan(kept.name.length)
    expect(renamed.size).toBe(index.size)
  })
})

describe('classifyLocalFiles', () => {
  it('pairs prebuilt longrange and bigWig without reading them whole or calling fetch', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
    const longrange = fileNamed('pet.longrange.gz', bgzfHeader())
    const index = fileNamed('pet.longrange.gz.tbi', tbiBytes())
    const signal = fileNamed(
      'signal.bigWig',
      new Uint8Array([0x26, 0xfc, 0x8f, 0x88, 0, 0, 0, 0]),
    )
    const result = await classifyLocalFiles([index, signal, longrange])
    expect(result.issues).toEqual([])
    expect(result.choices).toEqual([])
    expect(result.ready.map((track) => track.washUType)).toEqual(['bigwig', 'longrange'])
    const paired = result.ready[1]
    expect(paired?.readMethod).toBe(READ_TABIX_SLICE)
    expect(paired?.isText).toBe(false)
    expect(paired?.dataFile).toBe(longrange)
    expect(paired?.files?.[0]).toBe(longrange)
    expect(paired?.files?.[1]).toBe(index)
    expect(paired?.files?.[0]?.name.length).toBeLessThan(paired?.files?.[1]?.name.length ?? 0)
    expect(result.ready[0]?.readMethod).toBe(READ_BIGWIG_SLICE)
    expect(result.ready[0]?.files).toBeUndefined()
    expect(fetchSpy).not.toHaveBeenCalled()
  })

  it('reports a missing or wrong longrange index and does not convert the file', async () => {
    const missing = await classifyLocalFiles([
      fileNamed('pet.longrange.gz', bgzfHeader()),
    ])
    expect(missing.ready).toEqual([])
    expect(missing.issues.map((item) => item.reason)).toEqual(['missing-index'])

    const wrong = await classifyLocalFiles([
      fileNamed('pet.longrange.gz', bgzfHeader()),
      fileNamed('pet.longrange.gz.tbi', new Uint8Array([1, 2, 3, 4])),
    ])
    expect(wrong.ready).toEqual([])
    expect(wrong.issues.map((item) => item.reason)).toEqual(['wrong-index'])

    const raw = await classifyLocalFiles([
      fileNamed('pairs.bedpe', new TextEncoder().encode('chr1\t1\t2\tchr1\t3\t4\n')),
    ])
    expect(raw.ready).toEqual([])
    expect(raw.issues[0]?.reason).toBe('bedpe-not-converted')
  })

  it('asks when a gzip pair could be longrange or indexed BED', async () => {
    const result = await classifyLocalFiles([
      fileNamed('sample.gz', bgzfHeader()),
      fileNamed('sample.gz.tbi', tbiBytes()),
    ])
    expect(result.ready).toEqual([])
    expect(result.choices).toHaveLength(1)
    expect(result.choices[0]?.reason).toBe('ambiguous-format')
    expect(result.choices[0]?.options.map((option) => option.prepared.washUType)).toEqual([
      'longrange',
      'bed',
    ])
  })

  it('asks when several indexes match one longrange file', async () => {
    const result = await classifyLocalFiles([
      fileNamed('pet.longrange.gz', bgzfHeader()),
      fileNamed('pet.longrange.gz.tbi', tbiBytes()),
      fileNamed('pet.longrange.gz.tbi', tbiBytes()),
    ])
    expect(result.ready).toEqual([])
    expect(result.choices[0]?.reason).toBe('ambiguous-index')
    expect(result.choices[0]?.options).toHaveLength(2)
  })

  it('accepts BED3 and BED4 and rejects a short row without uploading', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
    const bed = fileNamed(
      'egfr.bed',
      new TextEncoder().encode(
        'track name=egfr\n# comment\nchr7\t55019017\t55019117\tEGFR\nchr7\t55020000\t55020100\n',
      ),
    )
    const result = await classifyLocalFiles([bed])
    expect(result.issues).toEqual([])
    expect(result.ready[0]).toMatchObject({
      format: 'BED4',
      washUType: 'bed',
      readMethod: READ_TEXT_FILE,
      isText: true,
      label: 'egfr.bed',
    })
    const broken = await classifyLocalFiles([
      fileNamed('bad.bed', new TextEncoder().encode('chr1\t10\n')),
    ])
    expect(broken.ready).toEqual([])
    expect(broken.issues[0]?.reason).toBe('plain-invalid')
    expect(fetchSpy).not.toHaveBeenCalled()
  })

  it('scans a plain peak in windows and refuses a file past the cap before reading it', async () => {
    const slice = vi.spyOn(File.prototype, 'slice')
    const line = 'chr1\t1\t2\tpeak\n'
    const repeats = Math.ceil((70 * 1024) / line.length)
    const text = line.repeat(repeats)
    const scanned = await classifyLocalFiles([
      fileNamed('wide.narrowPeak', new TextEncoder().encode(text)),
    ])
    expect(scanned.ready[0]?.format).toBe('narrowPeak')
    const windows = slice.mock.calls.map((call) => Number(call[1]) - Number(call[0]))
    expect(windows.length).toBeGreaterThan(1)
    expect(Math.max(...windows)).toBeLessThanOrEqual(64 * 1024)
    slice.mockRestore()

    const tooBig = new File([new ArrayBuffer(PLAIN_PEAK_MAX_BYTES + 1)], 'huge.bed')
    const sliceAgain = vi.spyOn(File.prototype, 'slice')
    const rejected = await classifyLocalFiles([tooBig])
    expect(rejected.ready).toEqual([])
    expect(rejected.issues[0]?.reason).toBe('plain-too-large')
    expect(sliceAgain).not.toHaveBeenCalled()
  })

  it('accepts an indexed BED that is larger than the plain-text cap', async () => {
    const padded = new Uint8Array(PLAIN_PEAK_MAX_BYTES + 50)
    padded.set(bgzfHeader())
    const result = await classifyLocalFiles([
      fileNamed('peaks.bed.gz', padded),
      fileNamed('peaks.bed.gz.tbi', tbiBytes()),
    ])
    expect(result.issues).toEqual([])
    expect(result.ready[0]).toMatchObject({
      format: 'indexed-bed',
      washUType: 'bed',
      readMethod: READ_TABIX_SLICE,
      isText: false,
    })
  })

  it('rejects a bigBed renamed as bigWig and an unmatched index', async () => {
    const result = await classifyLocalFiles([
      fileNamed('peaks.bw', new Uint8Array([0xeb, 0xf2, 0x89, 0x87])),
      fileNamed('orphan.tbi', tbiBytes()),
    ])
    expect(result.ready).toEqual([])
    expect(result.issues.map((item) => item.reason)).toEqual([
      'bigbed-not-bigwig',
      'unmatched-index',
    ])
  })
})
