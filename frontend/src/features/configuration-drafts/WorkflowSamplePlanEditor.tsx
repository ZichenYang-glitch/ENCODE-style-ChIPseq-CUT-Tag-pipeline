import type { ProjectInputRevisionResponse, ProjectSampleSummaryResponse } from '../../api/generated/models';
import type { SampleColumn } from '../input-workbench/schemaContract';
import {
  planFieldsFromSchema,
  readNumbersForLayout,
  type WorkflowSamplePlan,
} from './workflowSamplePlan';

export type InputCatalogStatus = 'loading' | 'empty' | 'error' | 'ready';

interface WorkflowSamplePlanEditorProps {
  columns: readonly SampleColumn[];
  plans: readonly WorkflowSamplePlan[];
  samples: readonly ProjectSampleSummaryResponse[];
  revisions: readonly ProjectInputRevisionResponse[];
  pageRevisions?: readonly ProjectInputRevisionResponse[];
  catalogStatus: InputCatalogStatus;
  catalogError: string | null;
  hasPrevious: boolean;
  hasNext: boolean;
  onAdd: (sampleRevisionId: string) => void;
  onRemove: (sampleRevisionId: string) => void;
  onChange: (sampleRevisionId: string, plan: WorkflowSamplePlan) => void;
  onPrevious?: () => void;
  onNext?: () => void;
  onRetry?: () => void;
}

export function WorkflowSamplePlanEditor({
  columns,
  plans,
  samples,
  revisions,
  pageRevisions,
  catalogStatus,
  catalogError,
  hasPrevious,
  hasNext,
  onAdd,
  onRemove,
  onChange,
  onPrevious,
  onNext,
  onRetry,
}: WorkflowSamplePlanEditorProps) {
  const fields = planFieldsFromSchema(columns);
  const listed = new Set(revisions.map((item) => item.input_file_revision_id));
  const selectedReads = plans.flatMap((plan) =>
    plan.reads.map((read) => read.input_file_revision_id).filter(Boolean),
  );
  const unlisted = selectedReads.filter((id) => catalogStatus === 'ready' && !listed.has(id));
  const available = samples.filter(
    (sample) => !plans.some((plan) => plan.sample_revision_id === sample.sample_revision_id),
  );
  const controls = plans.filter((plan) => plan.role === 'control');
  return (
    <section className="space-y-3" aria-labelledby="workflow-sample-plan-title">
      <h2 id="workflow-sample-plan-title" className="text-sm font-medium">
        Workflow sample plan
      </h2>
      <p className="text-xs text-[var(--color-text-muted)]">
        Choose sample and input revisions from this project. Paths are not entered here.
      </p>
      <div className="space-y-2" aria-label="Admitted input revisions">
        <h3 className="text-sm font-medium">Admitted input revisions</h3>
        {catalogStatus === 'loading' && <p role="status">Loading admitted input revisions.</p>}
        {catalogStatus === 'empty' && (
          <p role="status">This project has no admitted input revisions.</p>
        )}
        {catalogStatus === 'error' && (
          <div className="space-y-2">
            <p role="alert">{catalogError ?? 'The admitted input catalog could not be loaded.'}</p>
            <button type="button" className="min-h-11 text-sm underline" onClick={() => onRetry?.()}>
              Try the input catalog again
            </button>
          </div>
        )}
        {catalogStatus === 'ready' && (
          <ul className="space-y-1">
            {(pageRevisions ?? revisions).map((item) => (
              <li key={item.input_file_revision_id} className="text-sm">
                {item.stable_key} · revision {item.revision_number}
              </li>
            ))}
          </ul>
        )}
        {unlisted.length > 0 && (
          <p role="status">
            {unlisted.length === 1
              ? 'A selected input revision is no longer listed.'
              : `${unlisted.length} selected input revisions are no longer listed.`}
          </p>
        )}
        {(hasPrevious || hasNext) && (
          <div className="flex gap-2">
            <button type="button" disabled={!hasPrevious} onClick={() => onPrevious?.()}>
              Previous input revisions
            </button>
            <button type="button" disabled={!hasNext} onClick={() => onNext?.()}>
              Next input revisions
            </button>
          </div>
        )}
      </div>
      <label className="block text-sm">
        Add a sample revision to the plan
        <select
          className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
          value=""
          onChange={(event) => {
            if (event.target.value) onAdd(event.target.value);
          }}
        >
          <option value="">Select a project sample</option>
          {available.map((sample) => (
            <option key={sample.sample_revision_id} value={sample.sample_revision_id}>
              {sample.display_name}
            </option>
          ))}
        </select>
      </label>
      {plans.map((plan, index) => {
        const sample = samples.find((item) => item.sample_revision_id === plan.sample_revision_id);
        return (
          <article key={plan.sample_revision_id} className="space-y-2" aria-label={`Sample plan ${index + 1}`}>
            <p className="text-sm font-medium">{sample?.display_name ?? 'Sample revision no longer listed'}</p>
            {fields.map((field) =>
              field.enumValues ? (
                <label key={field.key} className="block text-sm">
                  {field.label}
                  <select
                    aria-label={`${sample?.display_name ?? 'Sample'} ${field.label}`}
                    className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
                    value={planFieldValue(plan, field.key)}
                    onChange={(event) =>
                      onChange(
                        plan.sample_revision_id,
                        planWithField(plan, field.key, event.target.value),
                      )
                    }
                  >
                    <option value="">Select…</option>
                    {field.enumValues.filter(Boolean).map((option) => (
                      <option key={option} value={option}>
                        {option}
                      </option>
                    ))}
                  </select>
                </label>
              ) : (
                <label key={field.key} className="block text-sm">
                  {field.label}
                  <input
                    aria-label={`${sample?.display_name ?? 'Sample'} ${field.label}`}
                    className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
                    value={planFieldValue(plan, field.key)}
                    onChange={(event) =>
                      onChange(
                        plan.sample_revision_id,
                        planWithField(plan, field.key, event.target.value),
                      )
                    }
                  />
                </label>
              ),
            )}
            {readNumbersForLayout(plan.layout).map((read) => (
              <label key={read} className="block text-sm">
                Read {read}
                <select
                  aria-label={`${sample?.display_name ?? 'Sample'} read ${read}`}
                  className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
                  value={plan.reads.find((item) => item.read === read)?.input_file_revision_id ?? ''}
                  onChange={(event) => {
                    const next = readNumbersForLayout(plan.layout).map((number) => ({
                      read: number,
                      input_file_revision_id:
                        number === read
                          ? event.target.value
                          : plan.reads.find((item) => item.read === number)?.input_file_revision_id ?? '',
                    }));
                    onChange(plan.sample_revision_id, { ...plan, reads: next });
                  }}
                >
                  <option value="">Select an admitted input revision</option>
                  {revisions.map((item) => (
                    <option key={item.input_file_revision_id} value={item.input_file_revision_id}>
                      {item.stable_key}
                    </option>
                  ))}
                </select>
              </label>
            ))}
            {plan.role !== 'control' && (
              <label className="block text-sm">
                Control sample revision
                <select
                  aria-label={`${sample?.display_name ?? 'Sample'} control sample revision`}
                  className="mt-1 block h-11 w-full rounded-[4px] border border-[var(--color-border-strong)] bg-[var(--color-surface)] px-2 text-sm sm:h-9"
                  value={plan.control_sample_revision_id ?? ''}
                  onChange={(event) =>
                    onChange(plan.sample_revision_id, {
                      ...plan,
                      control_sample_revision_id: event.target.value || null,
                    })
                  }
                >
                  <option value="">No control sample</option>
                  {controls
                    .filter((item) => item.sample_revision_id !== plan.sample_revision_id)
                    .map((item) => (
                      <option key={item.sample_revision_id} value={item.sample_revision_id}>
                        {samples.find((sample) => sample.sample_revision_id === item.sample_revision_id)
                          ?.display_name ?? item.sample_revision_id}
                      </option>
                    ))}
                </select>
              </label>
            )}
            <button type="button" className="min-h-11 text-sm underline" onClick={() => onRemove(plan.sample_revision_id)}>
              Remove sample plan {index + 1}
            </button>
          </article>
        );
      })}
    </section>
  );
}

function planFieldValue(plan: WorkflowSamplePlan, key: string): string {
  if (key === 'assay' || key === 'target' || key === 'peak_mode' || key === 'layout' || key === 'role') {
    return plan[key];
  }
  return '';
}

function planWithField(plan: WorkflowSamplePlan, key: string, value: string): WorkflowSamplePlan {
  if (key === 'role' && value === 'control') {
    return { ...plan, role: value, control_sample_revision_id: null };
  }
  if (key === 'assay' || key === 'target' || key === 'peak_mode' || key === 'layout' || key === 'role') {
    return { ...plan, [key]: value };
  }
  return plan;
}
