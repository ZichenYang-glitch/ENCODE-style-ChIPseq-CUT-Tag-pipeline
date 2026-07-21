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
    schema_version: '1.1.0',
    upstream_identity: null,
    availability: {
      authoring: 'available',
      execution: 'available',
      reason_code: 'WORKFLOW_EXECUTION_READY',
    },
  },
  {
    metadata: {
      workflow_id: 'bulk-rnaseq',
      name: 'Bulk RNA-seq',
      version: '0.3.0',
      description:
        'Offline-composable HelixWeave adapter for pinned nf-core/rnaseq 3.26.0.',
      engines: ['nextflow'],
      tags: ['bulk-rnaseq', 'illumina', 'nf-core'],
    },
    schema_version: '1.0.0',
    capabilities: {
      supports: ['validation', 'input_authoring'],
    },
    upstream_identity: {
      name: 'nf-core/rnaseq',
      version: '3.26.0',
      revision: 'e7ca46272c8f9d5ceee3f71759f4ba551d3217a4',
    },
    availability: {
      authoring: 'available',
      execution: 'not_configured',
      reason_code: 'WORKFLOW_EXECUTION_NOT_CONFIGURED',
    },
  },
];

export const stubWorkflowSchemas: Record<string, WorkflowSchema> = {
  'encode-style-chipseq-cuttag-atac-mnase': {
    schema_version: '1.1.0',
    input_modes: {
      config: ['object'],
      samples: ['inline_rows', 'server_path'],
      options: ['object'],
    },
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
  'bulk-rnaseq': {
    schema_version: '1.0.0',
    input_modes: {
      config: ['object'],
      samples: ['inline_rows'],
      options: ['object'],
    },
    config_schema: {
      type: 'object',
      title: 'Bulk RNA-seq adapter config',
      description:
        'The complete adapter-owned contract is loaded by the authoring workspace.',
    },
    sample_schema: {
      type: 'array',
      title: 'Bulk RNA-seq sample rows',
      description:
        'Inline rows support SE, PE, repeated lanes, strandedness, and fixed ILLUMINA platform identity.',
    },
    option_schema: {
      type: 'object',
      title: 'Bulk RNA-seq adapter options',
      additionalProperties: false,
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
