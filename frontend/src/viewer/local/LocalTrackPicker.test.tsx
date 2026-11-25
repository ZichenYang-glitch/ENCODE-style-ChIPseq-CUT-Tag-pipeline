import { render, screen, waitFor, within } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { describe, expect, it, vi } from 'vitest'
import { indexedPair, indexedRecords, signalBigWig } from './indexFixture'
import LocalTrackPicker, { type LocalWashuTrack } from './LocalTrackPicker'

function fileBytes(bytes: Uint8Array): ArrayBuffer {
  const copy = new ArrayBuffer(bytes.byteLength)
  new Uint8Array(copy).set(bytes)
  return copy
}

function localFile(name: string, bytes: Uint8Array): File {
  return new File([fileBytes(bytes)], name)
}

function bgzfHeader(): Uint8Array {
  return new Uint8Array([
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00,
  ])
}

describe('LocalTrackPicker', () => {
  it('adds longrange, bigWig, and BED, then removes one track', async () => {
    const user = userEvent.setup()
    const seen: LocalWashuTrack[][] = []
    render(
      <LocalTrackPicker
        assembly="hg38"
        onTracksChange={(tracks) => {
          seen.push(tracks)
        }}
      />,
    )
    expect(screen.getByTestId('local-track-assembly')).toHaveTextContent('用户指定 GRCh38')
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockRejectedValue(new Error('network'))
    const built = indexedPair([
      'chr7\t55019017\t55019117\tchr7:55210000-55211000,1',
    ])
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('pet.longrange.gz', built.bgzf),
      localFile('pet.longrange.gz.tbi', built.tabix),
      localFile('signal.bw', signalBigWig()),
      localFile('egfr.bed', new TextEncoder().encode('chr7\t55019017\t55019117\tEGFR\n')),
    ])
    await waitFor(() => expect(screen.getAllByTestId('local-track')).toHaveLength(3))
    const rows = screen.getAllByTestId('local-track')
    expect(rows.map((row) => row.getAttribute('data-track-type'))).toEqual([
      'longrange',
      'bigwig',
      'bed',
    ])
    expect(rows.map((row) => row.getAttribute('data-read-method'))).toEqual([
      'wuepgg3-track-tabix-file-slice',
      'wuepgg3-track-bigwig-file-slice',
      'wuepgg3-track-text-file',
    ])
    const latest = seen.at(-1) ?? []
    expect(latest).toHaveLength(3)
    for (const track of latest) {
      expect(track.assembly).toBe('hg38')
      expect(track.model.url).toBe('')
      expect(track.model.fileObj).toBeInstanceOf(File)
      expect(track.model.metadata.localOnly).toBe(true)
      expect(track.model.metadata.assemblySignature).toBe('用户指定 GRCh38')
      expect(track.model.name).toBe(track.label)
      expect(JSON.stringify(track.model.metadata)).not.toMatch(/[/\\]home[/\\]/)
    }
    const longrange = latest.find((track) => track.model.type === 'longrange')
    expect(longrange?.model.files).toHaveLength(2)
    expect(longrange?.model.isText).toBe(false)
    expect(longrange?.model.options).toEqual({ clampHeight: true })
    expect(fetchSpy).not.toHaveBeenCalled()

    await user.click(screen.getByRole('button', { name: 'Remove local track signal.bw' }))
    await waitFor(() => expect(screen.getAllByTestId('local-track')).toHaveLength(2))
    expect(screen.queryByText('signal.bw')).not.toBeInTheDocument()
    expect(seen.at(-1)?.map((track) => track.label)).toEqual(['pet.longrange.gz', 'egfr.bed'])
    fetchSpy.mockRestore()
  })

  it('keeps an existing file on its original assembly when the reference changes', async () => {
    const user = userEvent.setup()
    const onTracksChange = vi.fn()
    const view = render(
      <LocalTrackPicker assembly="hg38" onTracksChange={onTracksChange} />,
    )
    await user.upload(
      screen.getByLabelText('Local track files'),
      localFile('egfr.bed', new TextEncoder().encode('chr7\t55019017\t55019117\n')),
    )
    await waitFor(() => expect(screen.getByTestId('local-track')).toBeInTheDocument())
    const added = onTracksChange.mock.calls.at(-1)?.[0] as LocalWashuTrack[]
    view.rerender(<LocalTrackPicker assembly="mm10" onTracksChange={onTracksChange} />)
    const row = screen.getByTestId('local-track')
    expect(row).toHaveAttribute('data-assembly', 'hg38')
    expect(within(row).getByText('用户指定 GRCh38')).toBeInTheDocument()
    expect(screen.getByTestId('local-track-assembly')).toHaveTextContent('用户指定 GRCm38')
    expect(screen.getByRole('status')).toHaveAttribute('data-reason', 'assembly-not-relabeled')
    expect(screen.getByRole('status')).toHaveTextContent('未改标')
    const after = onTracksChange.mock.calls.at(-1)?.[0] as LocalWashuTrack[]
    expect(after).toHaveLength(1)
    expect(after[0]?.assembly).toBe('hg38')
    expect(after[0]?.label).toBe('egfr.bed')
    expect(after[0]?.model.name).toBe('egfr.bed')
    expect(after[0]?.model.metadata.genome).toBe('hg38')
    expect(after[0]?.id).toBe(added[0]?.id)
  })

  it('lets the user choose an ambiguous pair and recover by skipping it', async () => {
    const user = userEvent.setup()
    const onTracksChange = vi.fn()
    render(<LocalTrackPicker assembly="hg38" onTracksChange={onTracksChange} />)
    const ambiguous = indexedPair(['chr1\t10\t20\tchr1:30-40,1'])
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('sample.gz', ambiguous.bgzf),
      localFile('sample.gz.tbi', ambiguous.tabix),
    ])
    const choice = await screen.findByTestId('local-track-choice')
    expect(choice).toHaveAttribute('data-reason', 'ambiguous-format')
    expect(screen.queryByTestId('local-track')).not.toBeInTheDocument()
    await user.click(within(choice).getByRole('button', { name: /作为 longrange/ }))
    const row = await screen.findByTestId('local-track')
    expect(row).toHaveAttribute('data-track-type', 'longrange')
    expect(row).toHaveAttribute('data-read-method', 'wuepgg3-track-tabix-file-slice')

    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('other.gz', ambiguous.bgzf),
      localFile('other.gz.tbi', ambiguous.tabix),
    ])
    const second = await screen.findByTestId('local-track-choice')
    await user.click(within(second).getByRole('button', { name: '不添加' }))
    await waitFor(() => expect(screen.queryByTestId('local-track-choice')).not.toBeInTheDocument())
    expect(screen.getAllByTestId('local-track')).toHaveLength(1)
  })

  it('lets a peak misidentified as longrange be retried as indexed BED', async () => {
    const user = userEvent.setup()
    render(<LocalTrackPicker assembly="hg38" onTracksChange={() => undefined} />)
    const peaks = indexedPair(['chr1\t10\t20\tpeak'])
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('sample.gz', peaks.bgzf),
      localFile('sample.gz.tbi', peaks.tabix),
      localFile('signal.bw', signalBigWig()),
    ])
    const choice = await screen.findByTestId('local-track-choice')
    await user.click(within(choice).getByRole('button', { name: /作为 longrange/ }))
    expect(await within(choice).findByRole('alert')).toHaveAttribute('data-reason', 'invalid-longrange')
    expect(screen.getByTestId('local-track')).toHaveTextContent('signal.bw')
    await user.click(within(choice).getByRole('button', { name: /作为索引 BED/ }))
    await waitFor(() => expect(screen.getAllByTestId('local-track')).toHaveLength(2))
    expect(screen.queryByRole('alert')).not.toBeInTheDocument()
    expect(screen.queryByTestId('local-track-choice')).not.toBeInTheDocument()
  })

  it('withholds conflicting indexed files and keeps the assembly the user selected', async () => {
    const user = userEvent.setup()
    const onTracksChange = vi.fn()
    render(<LocalTrackPicker assembly="mm10" onTracksChange={onTracksChange} />)
    const goodPeak = indexedRecords({
      lines: ['chr1\t1000\t1200\tLocalPeak'],
      trailingNewline: false,
    })
    const goodPairs = indexedRecords({
      lines: ['chr1\t1000\t1200\tchr1:3000-3200,1'],
    })
    const foreignPeak = indexedRecords({
      header: 'track db=mm39',
      lines: ['chr1\t1000\t1200\tForeignPeak'],
    })
    const foreignPairs = indexedRecords({
      header: 'browser genome hg38',
      lines: ['chr1\t1000\t1200\tchr1:3000-3200,1'],
      trailingNewline: false,
    })
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('foreign.bed', new TextEncoder().encode('track db=mm39\nchr1\t10\t20\tForeign')),
      localFile('local.bed', new TextEncoder().encode('chr1\t10\t20\tLocal')),
      localFile('foreign.bed.gz', foreignPeak.bgzf),
      localFile('foreign.bed.gz.tbi', foreignPeak.tabix),
      localFile('local.bed.gz', goodPeak.bgzf),
      localFile('local.bed.gz.tbi', goodPeak.tabix),
      localFile('foreign.longrange.gz', foreignPairs.bgzf),
      localFile('foreign.longrange.gz.tbi', foreignPairs.tabix),
      localFile('local.longrange.gz', goodPairs.bgzf),
      localFile('local.longrange.gz.tbi', goodPairs.tabix),
    ])
    await waitFor(() => expect(screen.getAllByTestId('local-track')).toHaveLength(3))
    expect(screen.getAllByRole('alert')).toHaveLength(3)
    for (const alert of screen.getAllByRole('alert')) {
      expect(alert).toHaveAttribute('data-reason', 'assembly-conflict')
    }
    for (const row of screen.getAllByTestId('local-track')) {
      expect(row).toHaveAttribute('data-assembly', 'mm10')
      expect(row).toHaveTextContent('用户指定 GRCm38')
      expect(row.textContent).not.toContain('foreign')
    }
    const handed = onTracksChange.mock.calls.at(-1)?.[0] as LocalWashuTrack[]
    expect(handed.map((track) => track.label)).toEqual([
      'local.bed',
      'local.bed.gz',
      'local.longrange.gz',
    ])
    expect(handed.every((track) => track.assembly === 'mm10')).toBe(true)
    expect(handed.every((track) => track.model.metadata.assemblySignature === '用户指定 GRCm38')).toBe(true)
  })

  it('blocks an indexed mm10 declaration opened on hg38 before it becomes a track', async () => {
    const user = userEvent.setup()
    const onTracksChange = vi.fn()
    render(<LocalTrackPicker assembly="hg38" onTracksChange={onTracksChange} />)
    const foreign = indexedRecords({
      header: 'track db=mm10',
      lines: ['chr7\t55103200\t55104100\tMousePeak'],
    })
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('mouse.bed.gz', foreign.bgzf),
      localFile('mouse.bed.gz.tbi', foreign.tabix),
      localFile('local.bed', new TextEncoder().encode('chr7\t10\t20\tLocalPeak')),
    ])
    const alert = await screen.findByRole('alert')
    expect(alert).toHaveAttribute('data-reason', 'assembly-conflict')
    const row = await screen.findByTestId('local-track')
    expect(row).toHaveAttribute('data-assembly', 'hg38')
    expect(row).toHaveTextContent('用户指定 GRCh38')
    expect(row).toHaveTextContent('local.bed')
    const handed = onTracksChange.mock.calls.at(-1)?.[0] as LocalWashuTrack[]
    expect(handed.map((track) => track.label)).toEqual(['local.bed'])
  })

  it('reports an assembly conflict when an ambiguous gzip is confirmed', async () => {
    const user = userEvent.setup()
    const onTracksChange = vi.fn()
    render(<LocalTrackPicker assembly="mm10" onTracksChange={onTracksChange} />)
    const built = indexedRecords({
      header: 'track db=hg38',
      lines: ['chr1\t1000\t1200\tchr1:3000-3200,1'],
    })
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('sample.gz', built.bgzf),
      localFile('sample.gz.tbi', built.tabix),
    ])
    const choice = await screen.findByTestId('local-track-choice')
    expect(choice).toHaveAttribute('data-reason', 'ambiguous-format')
    await user.click(within(choice).getByRole('button', { name: /作为索引 BED/ }))
    expect(await within(choice).findByRole('alert')).toHaveAttribute('data-reason', 'assembly-conflict')
    expect(screen.queryByTestId('local-track')).not.toBeInTheDocument()
    for (const call of onTracksChange.mock.calls) expect(call[0]).toEqual([])
  })

  it('shows a recoverable error for a missing index and keeps a valid peak', async () => {
    const user = userEvent.setup()
    render(
      <LocalTrackPicker assembly="hg38" onTracksChange={() => undefined} />,
    )
    await user.upload(screen.getByLabelText('Local track files'), [
      localFile('pet.longrange.gz', bgzfHeader()),
      localFile('peak.bed', new TextEncoder().encode('chr1\t10\t20\tpeak\n')),
    ])
    const alert = await screen.findByRole('alert')
    expect(alert).toHaveAttribute('data-reason', 'missing-index')
    expect(await screen.findByTestId('local-track')).toHaveAttribute('data-format', 'BED4')
    await user.click(screen.getByRole('button', { name: 'Dismiss local track message' }))
    expect(screen.queryByRole('alert')).not.toBeInTheDocument()
    expect(screen.getByTestId('local-track')).toBeInTheDocument()
  })
})
