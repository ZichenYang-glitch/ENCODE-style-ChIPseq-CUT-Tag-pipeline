import { describe, expect, it } from 'vitest';
import { createStubWorkflowClient } from './client';

describe('stub workflow client', () => {
  it('fails closed instead of inventing resolved Project/Sample provenance', async () => {
    const response = await createStubWorkflowClient().validateWorkflow(
      'encode-style-chipseq-cuttag-atac-mnase',
      {
        config: {},
        samples: 'samples.tsv',
        options: {},
      },
      {
        project_id: 'prj_11111111111111111111111111111111',
        sample_revision_ids: [
          'smpr_22222222222222222222222222222222',
        ],
      },
    );

    expect(response.ok).toBe(false);
    expect(response.snapshot).toBeNull();
    expect(response.issues).toEqual([
      expect.objectContaining({
        code: 'DATA_BINDING_SELECTION_UNAVAILABLE',
        path: 'project_id',
        source: 'registry',
      }),
    ]);
  });
});
