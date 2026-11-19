import { describe, expect, it } from 'vitest';
import {
  readCarriedHandoff,
  readDraftCarry,
  readExperimentBinding,
  sampleRevisionSetsEqual,
} from './handoff';

const PROJECT = `prj_${'1'.repeat(32)}`;
const LEGACY = `prj_${'0'.repeat(32)}`;
const REV_A = `smpr_${'a'.repeat(32)}`;
const REV_B = `smpr_${'b'.repeat(32)}`;
const REFERENCE = `refpr_${'1'.repeat(32)}`;
const WORKFLOW = 'encode-style-chipseq-cuttag-atac-mnase';

const saved = {
  workflowId: WORKFLOW,
  projectId: PROJECT,
  config: { threads: 4, outdir: 'results' },
  options: { strict_inputs: false },
  samples: [{ sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' }],
  sampleRevisionIds: [REV_A, REV_B],
  inputSelections: [],
  referenceProfileRevisionId: REFERENCE,
  workflowSamplesPresent: true,
};

describe('draft carry', () => {
  it('copies config, project, and sample revisions without an experiment', () => {
    const result = readDraftCarry(saved);
    expect(result.handoff).toMatchObject({
      workflowId: WORKFLOW,
      projectId: PROJECT,
      config: { threads: 4, outdir: 'results' },
      options: { strict_inputs: false },
      sampleRevisionIds: [REV_A, REV_B],
      referenceProfileRevisionId: REFERENCE,
      samples: [{ sample: 'S1', fastq_1: '/data/S1.fastq.gz', layout: 'SE' }],
    });
    expect(result.handoff).not.toHaveProperty('experimentRevisionId');
    expect(JSON.stringify(result.handoff)).not.toContain('workflow_samples');
    expect(result.notes.join(' ')).toMatch(/does not bind an experiment/);
    expect(result.notes.join(' ')).toMatch(/not an experiment binding/);
    expect(result.notes.join(' ')).toMatch(/workflow sample plan stays on the saved draft/);
    expect(result.notes.join(' ')).toMatch(/sample sheet was copied/);
  });

  it('keeps sample revision order and treats a reversed set as equal', () => {
    const result = readDraftCarry(saved);
    expect(result.handoff?.sampleRevisionIds).toEqual([REV_A, REV_B]);
    expect(
      sampleRevisionSetsEqual(result.handoff?.sampleRevisionIds ?? [], [REV_B, REV_A]),
    ).toBe(true);
  });

  it('refuses an empty, duplicate, or legacy project selection', () => {
    expect(readDraftCarry({ ...saved, sampleRevisionIds: [] }).handoff).toBeNull();
    expect(readDraftCarry({ ...saved, sampleRevisionIds: [REV_A, REV_A] }).handoff).toBeNull();
    expect(readDraftCarry({ ...saved, projectId: LEGACY }).handoff).toBeNull();
    expect(readDraftCarry({ ...saved, config: ['threads'] }).handoff).toBeNull();
  });

  it('does not copy a sample path or a workflow plan into the sheet', () => {
    const result = readDraftCarry({ ...saved, samples: '/data/samples.tsv' });
    expect(result.handoff?.samples).toBeNull();
    expect(result.handoff?.sampleRevisionIds).toEqual([REV_A, REV_B]);
    expect(result.notes.join(' ')).toMatch(/sample path was not copied/i);
  });

  it('rejects a carried handoff for another workflow', () => {
    const result = readDraftCarry(saved);
    expect(readCarriedHandoff({
      experimentValidationHandoff: result.handoff,
    }, 'bulk-rnaseq')).toBeNull();
    expect(readCarriedHandoff({
      experimentValidationHandoff: result.handoff,
    }, WORKFLOW)?.config).toEqual({ threads: 4, outdir: 'results' });
  });

  it('reads only an unresolved or frozen experiment binding', () => {
    expect(readExperimentBinding({ status: 'unresolved' })).toEqual({
      status: 'unresolved',
    });
    const revision = `expr_${'e'.repeat(32)}`;
    expect(readExperimentBinding({
      status: 'experiment_snapshot',
      experiment_revision_id: revision,
    })).toEqual({
      status: 'experiment_snapshot',
      experimentRevisionId: revision,
    });
    expect(readExperimentBinding({ status: 'experiment_snapshot' })).toBeNull();
    expect(sampleRevisionSetsEqual([REV_A, REV_B], [REV_A, `smpr_${'c'.repeat(32)}`])).toBe(false);
  });
});
