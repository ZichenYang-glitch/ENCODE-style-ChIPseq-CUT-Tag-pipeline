/**
 * @vitest-environment node
 *
 * Checks the local readers shipped in wuepgg 62.0.2 / wuepgg3-track 9.0.3,
 * then reads one window of a real BGZF longrange file and the fixture bigWig
 * through the same @gmod libraries, counting bytes from Blob.slice.
 */
import { spawnSync } from 'node:child_process'
import { createRequire } from 'node:module'
import { mkdtempSync, readFileSync, readdirSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'
import { afterAll, describe, expect, it, vi } from 'vitest'
import { looksLikeTabix } from './sliceRead'

const require = createRequire(import.meta.url)
const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../../../..')
const frontendRoot = path.join(repoRoot, 'frontend')
const fixtureDir = mkdtempSync(path.join(tmpdir(), 'local-track-'))

function arrayBuffer(bytes: Uint8Array): ArrayBuffer {
  const copy = new ArrayBuffer(bytes.byteLength)
  new Uint8Array(copy).set(bytes)
  return copy
}

class SliceCounter {
  readBytes = 0
  readFileBytes = 0

  constructor(private readonly bytes: Uint8Array) {}

  async read(
    buffer: Buffer,
    offset = 0,
    length = 0,
    position = 0,
  ): Promise<{ bytesRead: number; buffer: Buffer }> {
    if (!length) return { bytesRead: 0, buffer }
    const start = position
    const end = Math.min(this.bytes.length, start + length)
    const sliced = new Uint8Array(
      await new Blob([arrayBuffer(this.bytes)]).slice(start, end).arrayBuffer(),
    )
    this.readBytes += sliced.byteLength
    return { bytesRead: Buffer.from(sliced).copy(buffer, offset), buffer }
  }

  async readFile(): Promise<Buffer> {
    const sliced = new Uint8Array(
      await new Blob([arrayBuffer(this.bytes)]).slice(0, this.bytes.byteLength).arrayBuffer(),
    )
    this.readFileBytes += sliced.byteLength
    return Buffer.from(sliced)
  }
}

class BigWigSliceCounter {
  readBytes = 0
  maxEnd = 0

  constructor(private readonly bytes: Uint8Array) {}

  async read(length: number, position = 0): Promise<Uint8Array> {
    const start = position
    const end = Math.min(this.bytes.byteLength, start + length)
    this.maxEnd = Math.max(this.maxEnd, end)
    const sliced = new Uint8Array(
      await new Blob([arrayBuffer(this.bytes)]).slice(start, end).arrayBuffer(),
    )
    this.readBytes += sliced.byteLength
    return sliced
  }
}

function writeLongrangeFixture(): { data: Uint8Array; index: Uint8Array } {
  const dataPath = path.join(fixtureDir, 'pet.longrange.gz')
  const indexPath = path.join(fixtureDir, 'pet.longrange.gz.tbi')
  const script = `
import importlib.util, sys
spec = importlib.util.spec_from_file_location("browser_bgzf", sys.argv[1])
mod = importlib.util.module_from_spec(spec)
sys.modules["browser_bgzf"] = mod
spec.loader.exec_module(mod)
lines = []
for i in range(2500):
    start = 1000 + i * 100
    lines.append(f"chr1\\t{start}\\t{start + 40}\\tchr1:{start + 5000}-{start + 5100},1")
lines.append("chr7\\t55019017\\t55019117\\tchr7:55210000-55211000,1")
indexed = mod.write_indexed_bed(tuple(lines))
open(sys.argv[2], "wb").write(indexed.bgzf)
open(sys.argv[3], "wb").write(indexed.tabix)
`
  const result = spawnSync(
    'python3',
    [
      '-c',
      script,
      path.join(repoRoot, 'src/encode_pipeline/services/browser_bgzf.py'),
      dataPath,
      indexPath,
    ],
    { encoding: 'utf8' },
  )
  if (result.status !== 0) {
    throw new Error(result.stderr || result.stdout || 'could not build the longrange fixture')
  }
  return {
    data: new Uint8Array(readFileSync(dataPath)),
    index: new Uint8Array(readFileSync(indexPath)),
  }
}

afterAll(() => {
  rmSync(fixtureDir, { recursive: true, force: true })
})

describe('wuepgg local slice readers', () => {
  it('uses File.slice for tabix and bigWig and does not read a whole prebuilt longrange', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch')
    const wuepgg = JSON.parse(
      readFileSync(path.join(frontendRoot, 'node_modules/wuepgg/package.json'), 'utf8'),
    ) as { version: string }
    const trackPackage = JSON.parse(
      readFileSync(path.join(frontendRoot, 'node_modules/wuepgg3-track/package.json'), 'utf8'),
    ) as { version: string }
    const tabixPackage = JSON.parse(
      readFileSync(path.join(frontendRoot, 'node_modules/@gmod/tabix/package.json'), 'utf8'),
    ) as { version: string }
    const bbiPackage = JSON.parse(
      readFileSync(path.join(frontendRoot, 'node_modules/@gmod/bbi/package.json'), 'utf8'),
    ) as { version: string }
    expect(wuepgg.version).toBe('62.0.2')
    expect(trackPackage.version).toBe('9.0.3')
    expect(tabixPackage.version).toBe('1.6.1')
    expect(bbiPackage.version).toBe('7.1.0')

    const bundleName = readdirSync(path.join(frontendRoot, 'node_modules/wuepgg3-track/dist')).find(
      (name) => name.startsWith('index-') && name.endsWith('.mjs'),
    )
    expect(bundleName).toBeTruthy()
    const bundle = readFileSync(
      path.join(frontendRoot, 'node_modules/wuepgg3-track/dist', bundleName ?? ''),
      'utf8',
    )
    expect(bundle).toContain('this.blob.slice(t, s)')
    expect(bundle).toContain('this.blob.slice(n, i)')
    expect(bundle).toContain('D.fileObj instanceof File')
    expect(bundle).toContain('D.url === ""')

    const { data, index } = writeLongrangeFixture()
    expect(await looksLikeTabix(new Blob([arrayBuffer(index)]))).toBe(true)
    const { TabixIndexedFile } = require('@gmod/tabix') as {
      TabixIndexedFile: new (args: {
        filehandle: SliceCounter
        tbiFilehandle: SliceCounter
      }) => {
        getLines: (
          chrom: string,
          start: number,
          end: number,
          callback: (line: string) => void,
        ) => Promise<void>
      }
    }
    const dataHandle = new SliceCounter(data)
    const indexHandle = new SliceCounter(index)
    const tabix = new TabixIndexedFile({
      filehandle: dataHandle,
      tbiFilehandle: indexHandle,
    })
    const lines: string[] = []
    await tabix.getLines('chr7', 55_019_000, 55_020_000, (line) => {
      lines.push(line)
    })
    expect(lines).toEqual(['chr7\t55019017\t55019117\tchr7:55210000-55211000,1'])
    expect(dataHandle.readFileBytes).toBe(0)
    expect(dataHandle.readBytes).toBeGreaterThan(0)
    expect(dataHandle.readBytes).toBeLessThan(data.byteLength)
    expect(indexHandle.readFileBytes).toBe(index.byteLength)

    const signal = readFileSync(path.join(repoRoot, 'test/fixtures/browser/chr1-signal.bw'))
    const padded = new Uint8Array(signal.byteLength + 1024 * 1024)
    padded.set(signal)
    padded.fill(0xab, signal.byteLength)
    const signalHandle = new BigWigSliceCounter(padded)
    const { BigWig } = require('@gmod/bbi') as {
      BigWig: new (args: { filehandle: BigWigSliceCounter }) => {
        getFeatures: (
          chrom: string,
          start: number,
          end: number,
        ) => Promise<Array<{ start: number; end: number; score: number }>>
      }
    }
    const features = await new BigWig({ filehandle: signalHandle }).getFeatures('chr1', 0, 200_000)
    expect(features.length).toBeGreaterThan(0)
    expect(features.every((feature) => Number.isFinite(feature.score))).toBe(true)
    expect(signalHandle.readBytes).toBeGreaterThan(0)
    expect(signalHandle.readBytes).toBeLessThan(padded.byteLength / 2)
    expect(signalHandle.maxEnd).toBeLessThan(padded.byteLength)
    expect(fetchSpy).not.toHaveBeenCalled()

    const cloneDeep = require('lodash/cloneDeep.js') as (value: unknown) => {
      trackModelArr: Array<{ fileObj: File }>
    }
    const localFile = new File([arrayBuffer(data)], 'pet.longrange.gz')
    const localIndex = new File([arrayBuffer(index)], 'pet.longrange.gz.tbi')
    const model = {
      type: 'longrange',
      name: 'pet.longrange.gz',
      url: '',
      fileObj: localFile,
      files: [localFile, localIndex],
      isText: false,
      metadata: { genome: 'hg38' },
    }
    expect(() => cloneDeep({ trackModelArr: [model] })).not.toThrow()
    // Node 22 structuredClone turns File into Blob. The packaged Chromium
    // worker is where a local track still arrives as a File.
    const posted = structuredClone(model)
    expect(posted.url).toBe('')
    expect(posted.fileObj).toBeInstanceOf(Blob)
    expect(posted.files[0]).toBeInstanceOf(Blob)
    expect(posted.files[1]).toBeInstanceOf(Blob)
    fetchSpy.mockRestore()
  })
})
