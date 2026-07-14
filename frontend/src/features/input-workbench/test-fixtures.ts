import type { WorkflowSchemaResponse } from '../../api/generated/models';

export const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

export function createAuthoringSchemaFixture(): WorkflowSchemaResponse {
  return {
    schema_version: '1.1.0',
    schema_dialect: 'https://json-schema.org/draft/2020-12/schema',
    coverage: {
      config: 'partial',
      samples: 'complete',
      options: 'complete',
    },
    authoring_modes: {
      config: ['schema_form', 'yaml'],
      samples: ['tsv_upload', 'inline_table'],
      options: ['schema_form'],
    },
    input_modes: {
      config: ['object'],
      samples: ['inline_rows', 'server_path'],
      options: ['object'],
    },
    limits: {
      max_request_bytes: 2_097_152,
      max_sample_rows: 1_000,
      max_sample_columns: 64,
      max_sample_column_name_length: 128,
      max_sample_cell_length: 4_096,
    },
    config_schema: {
      $schema: 'https://json-schema.org/draft/2020-12/schema',
      $id: 'https://encode-pipeline.org/schemas/test/config/1.1.0',
      type: 'object',
      title: 'Test config',
      properties: {
        outdir: { type: 'string', const: 'results', default: 'results' },
        threads: { type: 'integer', minimum: 1, default: 8 },
        replicate_analysis: {
          type: 'object',
          title: 'Replicate analysis',
          default: { enabled: true },
          properties: {
            enabled: {
              type: 'boolean',
              title: 'Replicate analysis enabled',
              default: true,
            },
          },
          required: ['enabled'],
          additionalProperties: false,
        },
        chipseq_idr: {
          type: 'object',
          title: 'ChIP-seq IDR',
          default: { enabled: false },
          properties: {
            enabled: {
              type: 'boolean',
              title: 'ChIP-seq IDR enabled',
              default: false,
            },
          },
          required: ['enabled'],
          additionalProperties: false,
        },
        qc: {
          type: 'object',
          default: {},
          properties: {
            summary: { type: 'boolean', default: true },
          },
          additionalProperties: true,
        },
      },
      additionalProperties: true,
    },
    sample_schema: {
      $schema: 'https://json-schema.org/draft/2020-12/schema',
      $id: 'https://encode-pipeline.org/schemas/test/samples/1.1.0',
      type: 'array',
      minItems: 1,
      maxItems: 1_000,
      items: {
        type: 'object',
        required: ['sample', 'fastq_1', 'layout'],
        additionalProperties: false,
        properties: {
          sample: {
            type: 'string',
            minLength: 1,
            maxLength: 4_096,
            description: 'Sample identifier.',
          },
          fastq_1: {
            type: 'string',
            minLength: 1,
            maxLength: 4_096,
            description: 'R1 FASTQ path.',
          },
          layout: { type: 'string', enum: ['PE', 'SE'] },
          fastq_2: { type: 'string', maxLength: 4_096 },
        },
      },
    },
    option_schema: {
      $schema: 'https://json-schema.org/draft/2020-12/schema',
      $id: 'https://encode-pipeline.org/schemas/test/options/1.1.0',
      type: 'object',
      properties: {
        strict_inputs: { type: 'boolean', default: false },
      },
      additionalProperties: false,
    },
  };
}
