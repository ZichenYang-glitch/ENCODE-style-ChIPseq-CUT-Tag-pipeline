import { describe, expect, it } from 'vitest'
import { classifyLocalFiles } from './classifyLocalFiles'

describe('recovering one local selection while retaining its neighbours', () => {
  it('asks which renamed TBI to use even when there is only one candidate', async () => {
    const indexed = new File(['reads'], 'pet.longrange.gz')
    const index = new File(['invalid'], 'original-index.tbi')
    const result = await classifyLocalFiles([indexed, index])
    expect(result.ready).toEqual([])
    expect(result.issues).toEqual([])
    const [choice] = result.choices
    expect(choice.reason).toBe('ambiguous-index')
    expect(choice.options[0].label).toContain('original-index.tbi')
    expect(await choice.options[0].checkPair?.()).toBeTruthy()
  })

  it('rejects an explicitly conflicting BED reference and keeps the compatible file', async () => {
    const wrong = new File(['track name="peaks" db=hg19\nchr7\t5\t10\n'], 'other.bed')
    const good = new File(['browser genome GRCh38.p14\nchr7\t5\t10\n'], 'peaks.bed')
    const result = await classifyLocalFiles([wrong, good], 'hg38')
    expect(result.ready.map(item => item.label)).toEqual(['peaks.bed'])
    expect(result.issues[0]).toMatchObject({ reason: 'assembly-conflict' })
    expect(result.issues[0].message).toContain('hg19')
  })

  it('reports a file that disappeared from disk without discarding the readable selection', async () => {
    const stale = new File(['stale'], 'old.bw')
    Object.defineProperty(stale, 'slice', { value() { throw new Error('unavailable') } })
    const result = await classifyLocalFiles([stale, new File(['chr7\t5\t10\n'], 'peak.bed')])
    expect(result.ready.map(item => item.label)).toEqual(['peak.bed'])
    expect(result.issues[0].reason).toBe('unreadable-file')
  })
})
