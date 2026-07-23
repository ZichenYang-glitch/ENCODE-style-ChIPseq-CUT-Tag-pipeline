import { useEffect } from 'react';
import * as Tabs from '@radix-ui/react-tabs';
import { FileCode2, ListChecks, Settings2, TableProperties } from 'lucide-react';
import { useSearchParams } from 'react-router-dom';
import { Panel } from '../../components/Panel';
import type { WorkflowAvailability } from '../../api/types';
import { ConfigEditor, type ConfigMode } from './ConfigEditor';
import { DraftReview } from './DraftReview';
import { OptionsEditor } from './OptionsEditor';
import { SampleEditor } from './SampleEditor';
import type { WorkbenchSchema } from './schemaContract';
import { useInputDraft } from './useInputDraft';
import { ValidatedSubmission } from './ValidatedSubmission';

type WorkbenchStep = 'config' | 'samples' | 'options' | 'review';

const STEPS: Array<{
  value: WorkbenchStep;
  label: string;
  icon: typeof FileCode2;
}> = [
  { value: 'config', label: 'Config', icon: FileCode2 },
  { value: 'samples', label: 'Samples', icon: TableProperties },
  { value: 'options', label: 'Options', icon: Settings2 },
  { value: 'review', label: 'Review', icon: ListChecks },
];

function readStep(value: string | null): WorkbenchStep {
  return STEPS.some((step) => step.value === value)
    ? (value as WorkbenchStep)
    : 'config';
}

function readMode(value: string | null): ConfigMode {
  return value === 'yaml' ? 'yaml' : 'form';
}

interface InputWorkbenchProps {
  workflowId: string;
  schema: WorkbenchSchema;
  availability: WorkflowAvailability | null;
}

export function InputWorkbench({
  workflowId,
  schema,
  availability,
}: InputWorkbenchProps) {
  const [searchParams, setSearchParams] = useSearchParams();
  const draft = useInputDraft(schema);
  const step = readStep(searchParams.get('step'));
  const requestedMode = readMode(searchParams.get('mode'));
  const mode =
    draft.state.yamlIssue !== null && requestedMode === 'form'
      ? 'yaml'
      : requestedMode;

  useEffect(() => {
    if (mode === 'yaml') draft.syncYaml();
  }, [mode, draft]);

  function setStep(nextStep: string) {
    if (nextStep === 'review' && draft.state.yamlIssue !== null) return;
    const next = new URLSearchParams(searchParams);
    next.set('step', readStep(nextStep));
    if (nextStep !== 'config') next.delete('mode');
    setSearchParams(next);
  }

  function setMode(nextMode: ConfigMode) {
    if (nextMode === 'form' && draft.state.yamlIssue !== null) return;
    if (nextMode === 'yaml') draft.syncYaml();
    const next = new URLSearchParams(searchParams);
    next.set('step', 'config');
    next.set('mode', nextMode);
    setSearchParams(next);
  }

  return (
    <section className="input-workbench min-w-0 flex-1">
      <Panel title="Input authoring" className="min-w-0 overflow-hidden">
        <header className="flex min-w-0 flex-col gap-2 border-b border-[var(--color-border)] pb-3 sm:flex-row sm:items-start sm:justify-between">
          <div className="min-w-0">
            <h2 className="text-base font-semibold">Input workbench</h2>
            <p className="mt-1 break-all font-mono text-xs text-[var(--color-text-muted)]">
              {workflowId}
            </p>
          </div>
          <div className="flex flex-wrap gap-2 text-xs">
            <span className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2 py-1">
              Schema {schema.contract.schema_version}
            </span>
            <span className="rounded border border-amber-200 bg-[var(--color-warning-bg)] px-2 py-1 text-[var(--color-warning)]">
              Draft only · not scientifically validated
            </span>
          </div>
        </header>
        <p className="border-b border-[var(--color-border)] py-2 text-xs text-[var(--color-text-muted)]">
          This draft stays only in this page session; a browser refresh clears this draft.
        </p>

        <Tabs.Root value={step} onValueChange={setStep} className="min-w-0 pt-3">
          <Tabs.List
            className="grid grid-cols-2 gap-1 border-b border-[var(--color-border)] sm:flex"
            aria-label="Input workbench sections"
          >
            {STEPS.map(({ value, label, icon: Icon }) => (
              <Tabs.Trigger
                key={value}
                value={value}
                disabled={value === 'review' && draft.state.yamlIssue !== null}
                className="inline-flex min-w-0 items-center justify-center gap-1.5 border-b-2 border-transparent px-3 py-2 text-sm font-medium text-[var(--color-text-muted)] outline-none transition-colors hover:text-[var(--color-text)] focus-visible:ring-2 focus-visible:ring-[var(--color-accent)] data-[state=active]:border-[var(--color-accent)] data-[state=active]:text-[var(--color-accent)] disabled:cursor-not-allowed disabled:opacity-50"
              >
                <Icon aria-hidden="true" size={15} />
                {label}
              </Tabs.Trigger>
            ))}
          </Tabs.List>

          <Tabs.Content value="config" className="min-w-0 pt-4 outline-none">
            <ConfigEditor
              schema={schema}
              draft={draft}
              mode={mode}
              onModeChange={setMode}
            />
          </Tabs.Content>
          <Tabs.Content value="samples" className="min-w-0 pt-4 outline-none">
            <SampleEditor schema={schema} draft={draft} />
          </Tabs.Content>
          <Tabs.Content value="options" className="min-w-0 pt-4 outline-none">
            <OptionsEditor schema={schema} draft={draft} />
          </Tabs.Content>
          <Tabs.Content
            value="review"
            forceMount
            className="min-w-0 pt-4 outline-none data-[state=inactive]:hidden"
          >
            <DraftReview schema={schema} draft={draft} />
            <ValidatedSubmission
              workflowId={workflowId}
              draft={draft}
              availability={availability}
            />
          </Tabs.Content>
        </Tabs.Root>
      </Panel>
    </section>
  );
}
