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
    expect(await choice.options[0].prepare()).toHaveProperty('issue')
  })

  it.each([
    ['hg38', 'GRCh38.p14', 'hg19'],
    ['mm10', 'GRCm38.p6', 'mm39'],
  ])('rejects a conflicting BED and retains the compatible %s file', async (window, same, other) => {
    const wrong = new File([`track name="peaks" db=${other}\nchr7\t5\t10\n`], 'other.bed')
    const good = new File([`browser genome ${same}\nchr7\t5\t10\n`], 'peaks.bed')
    const result = await classifyLocalFiles([wrong, good], window)
    expect(result.ready.map(item => item.label)).toEqual(['peaks.bed'])
    expect(result.issues[0]).toMatchObject({ reason: 'assembly-conflict' })
    expect(result.issues[0].message).toContain(other)
  })

  it('reports a file that disappeared from disk without discarding the readable selection', async () => {
    const stale = new File(['stale'], 'old.bw')
    Object.defineProperty(stale, 'slice', { value() { throw new Error('unavailable') } })
    const result = await classifyLocalFiles([stale, new File(['chr7\t5\t10\n'], 'peak.bed')])
    expect(result.ready.map(item => item.label)).toEqual(['peak.bed'])
    expect(result.issues[0].reason).toBe('unreadable-file')
  })
})
