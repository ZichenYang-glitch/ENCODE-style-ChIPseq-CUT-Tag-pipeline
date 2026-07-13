import type {
  ValidateWorkflowResponse,
  WorkflowInputs,
  WorkflowSchema,
  WorkflowSummary,
} from '../api/types';

export const stubWorkflows: WorkflowSummary[] = [
  {
    metadata: {
      workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
      name: 'ENCODE-style ChIP-seq/CUT\u0026Tag/ATAC/MNase',
      version: '0.3.0',
      description: 'Adapter wrapper for the current ENCODE-style workflow.',
      engines: ['snakemake'],
      tags: ['chipseq', 'cuttag', 'atac', 'mnase', 'encode-style'],
    },
    capabilities: {
      supports: ['validation'],
    },
  },
];

export const stubWorkflowSchemas: Record<string, WorkflowSchema> = {
  'encode-style-chipseq-cuttag-atac-mnase': {
    config_schema: {
      schema_kind: 'hints',
      required: ['samples'],
      properties: {
        samples: {
          type: 'string',
          description: 'Path to the sample sheet TSV.',
        },
      },
    },
    sample_schema: {
      schema_kind: 'hints',
      format: 'tsv',
      required_columns: ['sample', 'fastq_1'],
    },
    option_schema: {
      schema_kind: 'hints',
      properties: {
        strict_inputs: {
          type: 'boolean',
          default: false,
        },
      },
    },
  },
};

export const stubValidationResponses: Record<
  string,
  (inputs: WorkflowInputs) => ValidateWorkflowResponse
> = {
  'encode-style-chipseq-cuttag-atac-mnase': (inputs) => {
    const issues = [];
    if (!inputs.samples) {
      issues.push({
        code: 'ENCODE_SAMPLES_INVALID',
        message: 'Sample sheet is invalid.',
        severity: 'error' as const,
        path: 'samples',
        source: 'samples',
        technical_message: 'Missing samples field.',
        hint: 'Provide a path to a TSV sample sheet.',
        context: {},
      });
    }
    if (issues.length > 0) {
      return {
        ok: false,
        workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
        value: null,
        snapshot: null,
        issues,
      };
    }
    return {
      ok: true,
      workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
      value: {
        config: inputs.config,
        samples: [],
      },
      snapshot: {
        snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
        workflow_id: 'encode-style-chipseq-cuttag-atac-mnase',
        schema_version: '1.0.0',
        adapter_version: '0.3.0',
        payload_digest: 'a'.repeat(64),
        validated_at: '2026-07-14T00:00:00Z',
        expires_at: '2026-07-14T00:30:00Z',
      },
      issues: [],
    };
  },
};
