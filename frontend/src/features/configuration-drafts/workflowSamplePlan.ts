import type { SampleColumn } from '../input-workbench/schemaContract';

export interface WorkflowSampleRead {
  read: number;
  input_file_revision_id: string;
}

export interface WorkflowSamplePlan {
  sample_revision_id: string;
  assay: string;
  target: string;
  peak_mode: string;
  layout: string;
  role: string;
  control_sample_revision_id: string | null;
  reads: WorkflowSampleRead[];
}

export const WORKFLOW_SAMPLE_SCHEMA_KEYS = [
  'assay',
  'target',
  'peak_mode',
  'layout',
  'role',
] as const;

export function planFieldsFromSchema(columns: readonly SampleColumn[]): SampleColumn[] {
  const byKey = new Map(columns.map((column) => [column.key, column]));
  return WORKFLOW_SAMPLE_SCHEMA_KEYS.flatMap((key) => {
    const column = byKey.get(key);
    return column ? [column] : [];
  });
}

export function readNumbersForLayout(layout: string): number[] {
  return layout === 'PE' ? [1, 2] : [1];
}

export function emptyWorkflowSamplePlan(sampleRevisionId: string): WorkflowSamplePlan {
  return {
    sample_revision_id: sampleRevisionId,
    assay: '',
    target: '',
    peak_mode: '',
    layout: '',
    role: '',
    control_sample_revision_id: null,
    reads: [{ read: 1, input_file_revision_id: '' }],
  };
}

export function parseWorkflowSamples(value: unknown): WorkflowSamplePlan[] {
  if (!Array.isArray(value)) return [];
  return value.flatMap((item) => {
    if (!item || typeof item !== 'object') return [];
    const record = item as Record<string, unknown>;
    if (typeof record.sample_revision_id !== 'string') return [];
    const reads = Array.isArray(record.reads)
      ? record.reads.flatMap((read) => {
          if (!read || typeof read !== 'object') return [];
          const row = read as Record<string, unknown>;
          if (typeof row.read !== 'number' || typeof row.input_file_revision_id !== 'string') {
            return [];
          }
          return [{ read: row.read, input_file_revision_id: row.input_file_revision_id }];
        })
      : [];
    return [
      {
        sample_revision_id: record.sample_revision_id,
        assay: typeof record.assay === 'string' ? record.assay : '',
        target: typeof record.target === 'string' ? record.target : '',
        peak_mode: typeof record.peak_mode === 'string' ? record.peak_mode : '',
        layout: typeof record.layout === 'string' ? record.layout : '',
        role: typeof record.role === 'string' ? record.role : '',
        control_sample_revision_id:
          typeof record.control_sample_revision_id === 'string'
            ? record.control_sample_revision_id
            : null,
        reads: reads.length > 0 ? reads : [{ read: 1, input_file_revision_id: '' }],
      },
    ];
  });
}

export function workflowSamplesForSave(plans: readonly WorkflowSamplePlan[]): WorkflowSamplePlan[] {
  return plans.map((plan) => ({
    ...plan,
    control_sample_revision_id:
      plan.role === 'control' ? null : plan.control_sample_revision_id,
    reads: readNumbersForLayout(plan.layout).map((read) => ({
      read,
      input_file_revision_id:
        plan.reads.find((item) => item.read === read)?.input_file_revision_id ?? '',
    })),
  }));
}

export function removeWorkflowSamplePlan(
  plans: readonly WorkflowSamplePlan[],
  sampleRevisionId: string,
): WorkflowSamplePlan[] {
  return plans
    .filter((plan) => plan.sample_revision_id !== sampleRevisionId)
    .map((plan) =>
      plan.control_sample_revision_id === sampleRevisionId
        ? { ...plan, control_sample_revision_id: null }
        : plan,
    );
}
