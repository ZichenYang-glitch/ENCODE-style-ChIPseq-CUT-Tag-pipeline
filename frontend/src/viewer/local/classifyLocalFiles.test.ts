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
import { indexedPair, signalBigWig } from './indexFixture'
import { readFileSlice } from './sliceRead'

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
    expect(assemblySignature('mm10')).toBe('用户指定 GRCm38')
    expect(assemblySignature('GRCm38.p6')).toBe('用户指定 GRCm38')
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
  it('pairs a real prebuilt longrange and a bigWig without calling fetch', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
    const built = indexedPair([
      'chr7\t55019017\t55019117\tchr7:55210000-55211000,1',
    ])
    const longrange = fileNamed('pet.longrange.gz', built.bgzf)
    const index = fileNamed('pet.longrange.gz.tbi', built.tabix)
    const signal = fileNamed('signal.bigWig', signalBigWig())
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

    const magicOnly = await classifyLocalFiles([
      fileNamed('pet.longrange.gz', bgzfHeader()),
      fileNamed('pet.longrange.gz.tbi', tbiBytes()),
    ])
    expect(magicOnly.ready).toEqual([])
    expect(magicOnly.issues.map((item) => item.reason)).toEqual(['wrong-index'])

    const raw = await classifyLocalFiles([
      fileNamed('pairs.bedpe', new TextEncoder().encode('chr1\t1\t2\tchr1\t3\t4\n')),
    ])
    expect(raw.ready).toEqual([])
    expect(raw.issues[0]?.reason).toBe('bedpe-not-converted')
  })

  it('asks when a gzip pair could be longrange or indexed BED', async () => {
    const built = indexedPair(['chr1\t10\t20\tpeak'])
    const result = await classifyLocalFiles([
      fileNamed('sample.gz', built.bgzf),
      fileNamed('sample.gz.tbi', built.tabix),
    ])
    expect(result.ready).toEqual([])
    expect(result.choices).toHaveLength(1)
    expect(result.choices[0]?.reason).toBe('ambiguous-format')
    expect(result.choices[0]?.options.map((option) => option.washUType)).toEqual([
      'longrange',
      'bed',
    ])
  })

  it('asks when several indexes match one longrange file', async () => {
    const built = indexedPair([
      'chr7\t55019017\t55019117\tchr7:55210000-55211000,1',
    ])
    const result = await classifyLocalFiles([
      fileNamed('pet.longrange.gz', built.bgzf),
      fileNamed('pet.longrange.gz.tbi', built.tabix),
      fileNamed('pet.longrange.gz.tbi', built.tabix),
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
    const built = indexedPair(['chr1\t10\t20\tpeak'])
    const padded = new Uint8Array(PLAIN_PEAK_MAX_BYTES + 50)
    padded.set(built.bgzf)
    const result = await classifyLocalFiles([
      fileNamed('peaks.bed.gz', padded),
      fileNamed('peaks.bed.gz.tbi', built.tabix),
    ])
    expect(result.issues).toEqual([])
    expect(result.ready[0]).toMatchObject({
      format: 'indexed-bed',
      washUType: 'bed',
      readMethod: READ_TABIX_SLICE,
      isText: false,
    })
  })

  it('rejects a real index that belongs to another track', async () => {
    const here = indexedPair(['chr7\t55019017\t55019117\tchr7:55210000-55211000,1'])
    const there = indexedPair(['chr1\t80000000\t80000100\tchr1:80010000-80010100,1'])
    const result = await classifyLocalFiles([
      fileNamed('pet.longrange.gz', here.bgzf),
      fileNamed('pet.longrange.gz.tbi', there.tabix),
    ])
    expect(result.ready).toEqual([])
    expect(result.choices).toEqual([])
    expect(result.issues.map((item) => item.reason)).toEqual(['index-mismatch'])
    expect(result.issues[0]?.message).toContain('不含整份数据摘要')
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

  it('rejects truncated bigWigs but preserves a valid neighbouring signal', async () => {
    const signal = signalBigWig()
    const result = await classifyLocalFiles([
      fileNamed('only-magic.bw', signal.subarray(0, 4)),
      fileNamed('missing-index.bw', signal.subarray(0, 170)),
      fileNamed('good.bw', signal),
    ])
    expect(result.ready.map(file => file.label)).toEqual(['good.bw'])
    expect(result.issues.map(issue => issue.reason)).toEqual(['invalid-bigwig', 'invalid-bigwig'])
  })

  it('gives the text reader the intervals that passed validation, retaining their names and coordinates', async () => {
    const row = 'chr1\t1020\t1080\tLOCAL_PEAK'
    const original = fileNamed('headers.bed', new TextEncoder().encode(
      'track name=peaks db=hg38\r\n# comment\r\n' + row + '\r\n\r\nchr1\t2000\t2100\t第二个峰',
    ))
    const [track] = (await classifyLocalFiles([original], 'hg38')).ready
    expect(track.dataFile.name).toBe(original.name)
    expect(new TextDecoder().decode(await readFileSlice(track.dataFile, 0, track.dataFile.size)))
      .toBe(row + '\r\nchr1\t2000\t2100\t第二个峰')
  })
})
