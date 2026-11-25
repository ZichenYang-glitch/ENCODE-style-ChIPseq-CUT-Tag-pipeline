/**
 * @vitest-environment node
 */
import { readFileSync } from 'node:fs'
import { afterEach, describe, expect, it, vi } from 'vitest'
import { PLAIN_PEAK_MAX_BYTES } from './contract'
import { classifyLocalFiles } from './classifyLocalFiles'
import { indexedRecords } from './indexFixture'
import { tabixPairProblem } from './tabixPair'

const reviewBytes = JSON.parse(
  readFileSync(new URL('./indexed-review-bytes.json', import.meta.url), 'utf8'),
) as Record<string, string>

function reviewFile(name: string): File {
  const encoded = reviewBytes[name]
  if (!encoded) throw new Error(`missing review fixture ${name}`)
  return new File([Uint8Array.from(Buffer.from(encoded, 'base64'))], name)
}

function fileNamed(name: string, bytes: Uint8Array): File {
  return new File([Uint8Array.from(bytes)], name)
}

afterEach(() => {
  vi.restoreAllMocks()
})

describe('htslib indexed files from the WashU review', () => {
  it('reads a complete final record whether or not it ends in LF', async () => {
    for (const name of [
      'final-line.bed.gz',
      'newline.bed.gz',
      'matching-eof.bed.gz',
      'matching-lf.bed.gz',
    ]) {
      const data = reviewFile(name)
      const verdict = await tabixPairProblem(data, reviewFile(`${name}.tbi`))
      expect(verdict.ok, name).toBe(true)
      if (verdict.ok) {
        expect(verdict.dataBytesRead, name).toBeGreaterThan(0)
        expect(verdict.dataBytesRead, name).toBeLessThan(64 * 1024)
        expect(verdict.maxDataOffset, name).toBeLessThanOrEqual(data.size)
      }
    }
    const swapped = await tabixPairProblem(
      reviewFile('final-line.bed.gz'),
      reviewFile('newline.bed.gz.tbi'),
    )
    expect(swapped.ok).toBe(false)
    if (!swapped.ok) expect(swapped.reason).toBe('index-mismatch')
  })

  it.each(['hg38', 'mm10'])(
    'keeps an undeclared final record on %s and blocks a conflicting db before a track exists',
    async (assembly) => {
      const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
      const names = [
        'human.bed.gz',
        'human.bed.gz.tbi',
        'mouse.bed.gz',
        'mouse.bed.gz.tbi',
        'pet.longrange.gz',
        'pet.longrange.gz.tbi',
        'mouse.longrange.gz',
        'mouse.longrange.gz.tbi',
        'mouse.bed',
        'final-line.bed.gz',
        'final-line.bed.gz.tbi',
        'mm39-peak.bed',
        'mm39-peak.bed.gz',
        'mm39-peak.bed.gz.tbi',
        'mm39-pet.longrange.gz',
        'mm39-pet.longrange.gz.tbi',
        'matching-eof.bed.gz',
        'matching-eof.bed.gz.tbi',
      ]
      const result = await classifyLocalFiles(names.map(reviewFile), assembly)
      const ready = result.ready.map((item) => item.label)
      const blocked = result.issues.map((item) => item.reason)
      if (assembly === 'hg38') {
        expect(ready).toEqual(['human.bed.gz', 'pet.longrange.gz', 'final-line.bed.gz', 'matching-eof.bed.gz'])
        expect(blocked).toEqual(Array(6).fill('assembly-conflict'))
        expect(result.issues.some((item) => item.message.includes('mouse.bed.gz') && item.message.includes('mm10'))).toBe(true)
      } else {
        expect(ready).toEqual([
          'mouse.bed.gz',
          'pet.longrange.gz',
          'mouse.longrange.gz',
          'mouse.bed',
          'final-line.bed.gz',
          'matching-eof.bed.gz',
        ])
        expect(blocked).toEqual(Array(4).fill('assembly-conflict'))
        expect(result.issues.some((item) => item.message.includes('human.bed.gz') && item.message.includes('hg38'))).toBe(true)
      }
      expect(result.choices).toEqual([])
      expect(result.issues.every((item) => item.reason === 'assembly-conflict')).toBe(true)
      expect(result.issues.some((item) => item.message.includes('mm39'))).toBe(true)
      expect(fetchSpy).not.toHaveBeenCalled()
    },
  )
})

describe('declared assemblies on plain BED, BGZF BED, and BGZF longrange', () => {
  const peak = 'chr1\t1000\t1200\tLocalPeak'
  const pair = 'chr1\t1000\t1200\tchr1:3000-3200,1'

  async function open(
    format: 'plain-bed' | 'bgzf-bed' | 'bgzf-longrange',
    header: string | undefined,
    assembly: string,
    trailingNewline = true,
  ) {
    if (format === 'plain-bed') {
      const text = `${header ? `${header}\n` : ''}${peak}${trailingNewline ? '\n' : ''}`
      return classifyLocalFiles([new File([text], 'peaks.bed')], assembly)
    }
    const built = indexedRecords({
      header,
      lines: [format === 'bgzf-longrange' ? pair : peak],
      trailingNewline,
    })
    const name = format === 'bgzf-longrange' ? 'pairs.longrange.gz' : 'peaks.bed.gz'
    return classifyLocalFiles(
      [fileNamed(name, built.bgzf), fileNamed(`${name}.tbi`, built.tabix)],
      assembly,
    )
  }

  it.each([
    ['plain-bed', undefined, 'hg38', true],
    ['plain-bed', undefined, 'mm10', true],
    ['plain-bed', 'track db=hg38', 'hg38', true],
    ['plain-bed', 'track db=GRCh38', 'hg38', true],
    ['plain-bed', 'track db=mm10', 'mm10', true],
    ['plain-bed', 'track db=GRCm38', 'mm10', true],
    ['plain-bed', 'browser genome GRCm38.p6', 'mm10', true],
    ['plain-bed', 'track db=mm10', 'hg38', false],
    ['plain-bed', 'track db=mm39', 'mm10', false],
    ['plain-bed', 'track db=hg38', 'mm10', false],
    ['plain-bed', 'browser genome mm39', 'mm10', false],
    ['bgzf-bed', undefined, 'hg38', true],
    ['bgzf-bed', undefined, 'mm10', true],
    ['bgzf-bed', 'track db=hg38', 'hg38', true],
    ['bgzf-bed', 'track db=GRCh38', 'hg38', true],
    ['bgzf-bed', 'track db=mm10', 'mm10', true],
    ['bgzf-bed', 'track db=GRCm38', 'mm10', true],
    ['bgzf-bed', 'browser genome GRCm38.p6', 'mm10', true],
    ['bgzf-bed', 'track db=mm10', 'hg38', false],
    ['bgzf-bed', 'track db=mm39', 'mm10', false],
    ['bgzf-bed', 'track db=mm39', 'hg38', false],
    ['bgzf-bed', 'track db=hg38', 'mm10', false],
    ['bgzf-longrange', undefined, 'hg38', true],
    ['bgzf-longrange', undefined, 'mm10', true],
    ['bgzf-longrange', 'track db=hg38', 'hg38', true],
    ['bgzf-longrange', 'track db=mm10', 'mm10', true],
    ['bgzf-longrange', 'track db=GRCm38', 'mm10', true],
    ['bgzf-longrange', 'track db=mm10', 'hg38', false],
    ['bgzf-longrange', 'track db=mm39', 'mm10', false],
    ['bgzf-longrange', 'browser genome hg38', 'mm10', false],
  ] as const)(
    '%s header %s on %s is drawable=%s',
    async (format, header, assembly, drawable) => {
      const result = await open(format, header, assembly)
      expect(result.choices).toEqual([])
      if (drawable) {
        expect(result.issues).toEqual([])
        expect(result.ready).toHaveLength(1)
      } else {
        expect(result.ready).toEqual([])
        expect(result.issues.map((item) => item.reason)).toEqual(['assembly-conflict'])
        expect(result.issues[0]?.message).toContain(assembly)
      }
    },
  )

  it('accepts a final line without LF and still blocks a conflicting declaration', async () => {
    for (const format of ['plain-bed', 'bgzf-bed', 'bgzf-longrange'] as const) {
      for (const assembly of ['hg38', 'mm10']) {
        const kept = await open(format, undefined, assembly, false)
        expect(kept.issues, `${format} ${assembly}`).toEqual([])
        expect(kept.ready, `${format} ${assembly}`).toHaveLength(1)
        const blocked = await open(format, 'track db=mm39', assembly, false)
        expect(blocked.ready, `${format} ${assembly}`).toEqual([])
        expect(blocked.issues.map((item) => item.reason), `${format} ${assembly}`).toEqual([
          'assembly-conflict',
        ])
      }
    }
  })

  it('reads the last of several indexed records when only that record lacks LF', async () => {
    for (const trailingNewline of [true, false]) {
      const built = indexedRecords({
        lines: [
          'chr1\t1000\t1100\tchr1:3000-3200,1',
          'chr1\t5000\t5100\tchr1:7000-7200,1',
        ],
        trailingNewline,
      })
      const verdict = await tabixPairProblem(
        fileNamed('pairs.longrange.gz', built.bgzf),
        fileNamed('pairs.longrange.gz.tbi', built.tabix),
        'longrange',
      )
      expect(verdict.ok, `trailingNewline=${trailingNewline}`).toBe(true)
      if (verdict.ok) expect(verdict.maxDataOffset).toBeLessThanOrEqual(built.bgzf.byteLength)
    }
  })

  it('still reports a missing index, a mismatched index, and an ambiguous name', async () => {
    const declared = indexedRecords({
      header: 'track db=mm39',
      lines: [peak],
    })
    const missing = await classifyLocalFiles([fileNamed('peaks.bed.gz', declared.bgzf)], 'mm10')
    expect(missing.ready).toEqual([])
    expect(missing.issues.map((item) => item.reason)).toEqual(['missing-index'])

    const other = indexedRecords({
      lines: ['chr1\t80000000\t80000100\tchr1:80010000-80010100,1'],
    })
    const mismatch = await classifyLocalFiles(
      [
        fileNamed('peaks.bed.gz', declared.bgzf),
        fileNamed('peaks.bed.gz.tbi', other.tabix),
      ],
      'mm10',
    )
    expect(mismatch.ready).toEqual([])
    expect(mismatch.issues.map((item) => item.reason)).toEqual(['index-mismatch'])

    const ambiguous = indexedRecords({
      header: 'track db=hg38',
      lines: [pair],
    })
    const choice = await classifyLocalFiles(
      [fileNamed('sample.gz', ambiguous.bgzf), fileNamed('sample.gz.tbi', ambiguous.tabix)],
      'mm10',
    )
    expect(choice.ready).toEqual([])
    expect(choice.issues).toEqual([])
    expect(choice.choices.map((item) => item.reason)).toEqual(['ambiguous-format'])
    const reasons = await Promise.all(
      choice.choices[0]?.options.map(async (option) => (await option.checkPair?.())?.reason) ?? [],
    )
    expect(reasons).toEqual(['assembly-conflict', 'assembly-conflict'])
  })

  it('treats equivalent header names as one assembly and blocks contradictory headers', async () => {
    const same = await open('bgzf-bed', 'track db=mm10\ntrack genome=GRCm38', 'mm10')
    expect(same.issues).toEqual([])
    expect(same.ready).toHaveLength(1)
    const split = await open('bgzf-longrange', 'track db=hg38\nbrowser genome mm10', 'hg38')
    expect(split.ready).toEqual([])
    expect(split.issues[0]?.reason).toBe('assembly-conflict')
    expect(split.issues[0]?.message).toContain('hg38')
    expect(split.issues[0]?.message).toContain('mm10')
  })

  it('checks assembly directives beyond one BGZF block and signals an unfinished header', async () => {
    // Real htslib pairs with 3,300 or 25,000 leading comment lines and a track directive.
    const verify = (name: string) => classifyLocalFiles(
      [reviewFile(name), reviewFile(name + '.tbi')], 'hg38',
    )
    const matching = await verify('long-header-to-hg38.bed.gz')
    expect(matching.ready).toHaveLength(1)
    expect(matching.issues).toEqual([])
    const foreign = await verify('long-header-to-mm10.bed.gz')
    expect(foreign.ready).toEqual([])
    expect(foreign.issues.map(item => item.reason)).toEqual(['assembly-conflict'])
    const incomplete = await verify('over-limit-header.bed.gz')
    expect(incomplete.ready).toEqual([])
    expect(incomplete.issues.map(item => item.reason)).toEqual(['unreadable-header'])
  })

  it('does not read the padded tail or call the network while checking a declaration', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
    const slice = vi.spyOn(File.prototype, 'slice')
    const built = indexedRecords({
      header: 'track db=hg38',
      lines: [pair],
      trailingNewline: false,
    })
    const padded = new Uint8Array(PLAIN_PEAK_MAX_BYTES + 4096)
    padded.set(built.bgzf)
    const blocked = await classifyLocalFiles(
      [
        fileNamed('pairs.longrange.gz', padded),
        fileNamed('pairs.longrange.gz.tbi', built.tabix),
      ],
      'mm10',
    )
    expect(blocked.ready).toEqual([])
    expect(blocked.issues.map((item) => item.reason)).toEqual(['assembly-conflict'])
    const allowed = await classifyLocalFiles(
      [
        fileNamed('pairs.longrange.gz', padded),
        fileNamed('pairs.longrange.gz.tbi', built.tabix),
      ],
      'GRCh38',
    )
    expect(allowed.issues).toEqual([])
    expect(allowed.ready.map((item) => item.label)).toEqual(['pairs.longrange.gz'])
    const ends = slice.mock.calls.map((call) => Number(call[1]))
    expect(Math.max(...ends)).toBeLessThan(padded.byteLength / 2)
    expect(fetchSpy).not.toHaveBeenCalled()
  })
})
