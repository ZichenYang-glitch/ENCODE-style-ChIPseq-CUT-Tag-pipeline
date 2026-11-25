/**
 * @vitest-environment node
 */
import { describe, expect, it, vi } from 'vitest'
import { indexedPair } from './indexFixture'
import { tabixPairProblem } from './tabixPair'

function fileNamed(name: string, bytes: Uint8Array): File {
  const copy = new ArrayBuffer(bytes.byteLength)
  new Uint8Array(copy).set(bytes)
  return new File([copy], name)
}

function chr1Line(start: number): string {
  return `chr1\t${start}\t${start + 40}\tchr1:${start + 1000}-${start + 1100},1`
}

function manyChr1(count: number): string[] {
  return Array.from({ length: count }, (_, index) => chr1Line(1_000_000 + index * 50))
}

describe('tabixPairProblem', () => {
  it('accepts a matching index from a bounded window and a raw TBI', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
    const built = indexedPair(manyChr1(12_000))
    const data = fileNamed('wide.longrange.gz', built.bgzf)
    const verdict = await tabixPairProblem(data, fileNamed('wide.longrange.gz.tbi', built.tabix))
    expect(verdict.ok).toBe(true)
    if (verdict.ok) {
      expect(verdict.dataBytesRead).toBeGreaterThan(0)
      expect(verdict.dataBytesRead).toBeLessThan(data.size / 2)
      expect(verdict.maxDataOffset).toBeLessThanOrEqual(data.size)
    }
    const raw = await tabixPairProblem(data, fileNamed('wide.longrange.gz.tbi', built.rawTabix))
    expect(raw.ok).toBe(true)
    expect(fetchSpy).not.toHaveBeenCalled()
    fetchSpy.mockRestore()
  })

  it('accepts the optional uint64 unplaced-record count emitted by htslib', async () => {
    const built = indexedPair([chr1Line(1_000_000)])
    const raw = new Uint8Array(built.rawTabix.length + 8)
    raw.set(built.rawTabix)
    new DataView(raw.buffer).setBigUint64(built.rawTabix.length, 42n, true)
    const data = fileNamed('htslib.longrange.gz', built.bgzf)
    expect((await tabixPairProblem(data, fileNamed('htslib.tbi', raw))).ok).toBe(true)
    expect((await tabixPairProblem(data, fileNamed('cut.tbi', raw.subarray(0, -1)))).ok).toBe(false)
  })

  it('accepts a sparse tabix chunk promoted into a parent bin by htslib', async () => {
    const built = indexedPair([chr1Line(1_000_000)])
    const raw = Uint8Array.from(built.rawTabix)
    const view = new DataView(raw.buffer)
    const firstBin = 36 + view.getUint32(32, true) + 4
    const leafBin = view.getUint32(firstBin, true)
    expect(leafBin).toBeGreaterThan(4680)
    view.setUint32(firstBin, (leafBin - 1) >>> 3, true)
    const result = await tabixPairProblem(
      fileNamed('sparse.longrange.gz', built.bgzf),
      fileNamed('sparse.longrange.gz.tbi', raw),
    )
    expect(result.ok).toBe(true)
  })

  it('does not read padding pasted after a matching BGZF payload', async () => {
    const built = indexedPair([chr1Line(1_000_000)])
    const padded = new Uint8Array(built.bgzf.byteLength + 1024 * 1024)
    padded.set(built.bgzf)
    padded.fill(0x5a, built.bgzf.byteLength)
    const verdict = await tabixPairProblem(
      fileNamed('pad.longrange.gz', padded),
      fileNamed('pad.longrange.gz.tbi', built.tabix),
    )
    expect(verdict.ok).toBe(true)
    if (verdict.ok) {
      expect(verdict.maxDataOffset).toBeLessThanOrEqual(built.bgzf.byteLength)
      expect(verdict.dataBytesRead).toBeLessThan(padded.byteLength / 2)
    }
  })

  it('rejects the wrong BGZF block end even when start, chromosome and bin coincide', async () => {
    const built = indexedPair([chr1Line(1_000_000)])
    const raw = Uint8Array.from(built.rawTabix)
    const view = new DataView(raw.buffer)
    const firstChunkEnd = 36 + view.getUint32(32, true) + 4 + 8 + 8
    view.setBigUint64(firstChunkEnd, 5n << 16n, true)
    const result = await tabixPairProblem(fileNamed('near.gz', built.bgzf), fileNamed('other.tbi', raw))
    expect(result.ok).toBe(false)
    if (!result.ok) expect(result.reason).toBe('index-mismatch')
  })

  it('rejects a real index from another chromosome without reading the whole data file', async () => {
    const dataBuilt = indexedPair(manyChr1(12_000))
    const other = indexedPair(['chr7\t55019017\t55019117\tchr7:55210000-55211000,1'])
    const data = fileNamed('chr1.longrange.gz', dataBuilt.bgzf)
    const verdict = await tabixPairProblem(
      data,
      fileNamed('chr1.longrange.gz.tbi', other.tabix),
    )
    expect(verdict.ok).toBe(false)
    if (!verdict.ok) {
      expect(verdict.reason).toBe('index-mismatch')
      expect(verdict.message).toContain('不含整份数据摘要')
      expect(verdict.dataBytesRead).toBeGreaterThan(0)
      expect(verdict.dataBytesRead).toBeLessThan(data.size / 2)
    }
  })

  it('rejects a real index for the same chromosome at a different interval', async () => {
    const near = indexedPair([chr1Line(1_000_000)])
    const far = indexedPair([chr1Line(80_000_000)])
    const verdict = await tabixPairProblem(
      fileNamed('near.longrange.gz', near.bgzf),
      fileNamed('near.longrange.gz.tbi', far.tabix),
    )
    expect(verdict.ok).toBe(false)
    if (!verdict.ok) expect(verdict.reason).toBe('index-mismatch')
  })

  it('rejects a truncated or magic-only TBI before reading past the data header', async () => {
    const built = indexedPair(manyChr1(4_000))
    const data = fileNamed('wide.longrange.gz', built.bgzf)
    const truncated = built.tabix.subarray(0, Math.max(24, Math.floor(built.tabix.byteLength / 2)))
    const cut = await tabixPairProblem(data, fileNamed('wide.longrange.gz.tbi', truncated))
    expect(cut.ok).toBe(false)
    if (!cut.ok) {
      expect(cut.reason).toBe('wrong-index')
      expect(cut.message).toContain('不含数据文件的整份摘要')
      expect(cut.dataBytesRead).toBeLessThanOrEqual(18)
    }
    const magic = await tabixPairProblem(
      data,
      fileNamed('wide.longrange.gz.tbi', new Uint8Array([0x54, 0x42, 0x49, 0x01, 0, 0, 0, 0])),
    )
    expect(magic.ok).toBe(false)
    if (!magic.ok) {
      expect(magic.reason).toBe('wrong-index')
      expect(magic.dataBytesRead).toBeLessThanOrEqual(18)
    }
    const damaged = new Uint8Array(built.tabix)
    damaged[Math.min(40, damaged.byteLength - 1)] ^= 0xff
    const corrupt = await tabixPairProblem(data, fileNamed('wide.longrange.gz.tbi', damaged))
    expect(corrupt.ok).toBe(false)
    if (!corrupt.ok) expect(['wrong-index', 'index-mismatch']).toContain(corrupt.reason)
  })

  it('rejects a long index whose virtual offsets do not fit the data file', async () => {
    const small = indexedPair([chr1Line(1_000_000)])
    const wide = indexedPair(manyChr1(12_000))
    const verdict = await tabixPairProblem(
      fileNamed('small.longrange.gz', small.bgzf),
      fileNamed('small.longrange.gz.tbi', wide.tabix),
    )
    expect(verdict.ok).toBe(false)
    if (!verdict.ok) {
      expect(verdict.reason).toBe('index-mismatch')
      expect(verdict.dataBytesRead).toBeLessThanOrEqual(18)
    }
  })
})
