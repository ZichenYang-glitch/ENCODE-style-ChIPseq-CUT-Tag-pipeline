import { useEffect, useMemo, useState } from 'react';
import * as Tabs from '@radix-ui/react-tabs';
import { ArrowLeft, FileCode2, ListChecks, Settings2, TableProperties } from 'lucide-react';
import { Link, useSearchParams } from 'react-router-dom';
import { Panel } from '../../components/Panel';
import type { WorkflowAvailability } from '../../api/types';
import type { ReferenceProfileSummary } from '../../api/runTypes';
import { ConfigEditor, type ConfigMode } from './ConfigEditor';
import { DraftReview } from './DraftReview';
import { OptionsEditor } from './OptionsEditor';
import { ReferenceProfileSelector } from './ReferenceProfileSelector';
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
  workflowName?: string;
  schema: WorkbenchSchema;
  availability: WorkflowAvailability | null;
  referenceProfiles: readonly ReferenceProfileSummary[] | null;
  referenceProfilesLoading: boolean;
  referenceProfilesRefreshing: boolean;
  referenceProfilesError: string | null;
  onRefreshReferenceProfiles: () => void;
}

export function InputWorkbench({
  workflowId,
  workflowName,
  schema,
  availability,
  referenceProfiles,
  referenceProfilesLoading,
  referenceProfilesRefreshing,
  referenceProfilesError,
  onRefreshReferenceProfiles,
}: InputWorkbenchProps) {
  const [searchParams, setSearchParams] = useSearchParams();
  const draft = useInputDraft(schema);
  const [referenceSelectionMessage, setReferenceSelectionMessage] = useState<
    string | null
  >(null);
  const step = readStep(searchParams.get('step'));
  const requestedMode = readMode(searchParams.get('mode'));
  const mode =
    draft.state.yamlIssue !== null && requestedMode === 'form'
      ? 'yaml'
      : requestedMode;

  useEffect(() => {
    if (mode === 'yaml') draft.syncYaml();
  }, [mode, draft]);

  const selectedReference = useMemo(
    () =>
      referenceProfiles?.find(
        (profile) =>
          profile.revision_id === draft.state.referenceProfileRevisionId,
      ) ?? null,
    [referenceProfiles, draft.state.referenceProfileRevisionId],
  );

  useEffect(() => {
    if (referenceProfiles === null) return;
    const selectedRevisionId = draft.state.referenceProfileRevisionId;
    if (
      selectedRevisionId !== null &&
      !referenceProfiles.some(
        (profile) => profile.revision_id === selectedRevisionId,
      )
    ) {
      draft.setReferenceProfileRevisionId(null);
      setReferenceSelectionMessage(
        'The selected reference is no longer enabled or current. Select and validate an available revision again.',
      );
      return;
    }
    if (
      selectedRevisionId === null &&
      referenceProfiles.length === 1 &&
      referenceSelectionMessage === null
    ) {
      draft.setReferenceProfileRevisionId(referenceProfiles[0].revision_id);
      setReferenceSelectionMessage(null);
    }
  }, [
    draft.setReferenceProfileRevisionId,
    draft.state.referenceProfileRevisionId,
    referenceSelectionMessage,
    referenceProfiles,
  ]);

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

  const stepReady: Record<WorkbenchStep, boolean> = {
    config:
      draft.configValid &&
      draft.state.yamlIssue === null &&
      draft.state.configFormIssue === null,
    samples: draft.samplesValid,
    options: draft.optionsValid && draft.state.optionsFormIssue === null,
    review: draft.reviewReady,
  };

  return (
    <section className="input-workbench min-w-0 flex-1">
      <Panel title="Input authoring" className="min-w-0 overflow-hidden">
        <header className="flex min-w-0 flex-col gap-2 border-b border-[var(--color-border)] pb-3 sm:flex-row sm:items-start sm:justify-between">
          <div className="min-w-0">
            <Link
              className="inline-flex min-h-9 items-center gap-1 text-xs font-medium text-[var(--color-link)] hover:underline"
              to={`/workflows/${encodeURIComponent(workflowId)}`}
            >
              <ArrowLeft aria-hidden="true" size={14} />
              Workflow details
            </Link>
            <h2 className="text-xl font-semibold leading-7">Input workbench</h2>
            {workflowName && (
              <p className="mt-1 text-sm font-medium text-[var(--color-text)]">
                {workflowName}
              </p>
            )}
            <p className="mt-0.5 truncate font-mono text-xs text-[var(--color-text-faint)]" title={workflowId}>
              {workflowId}
            </p>
          </div>
          <div className="flex flex-wrap gap-2 text-xs">
            <span className="rounded border border-[var(--color-border)] bg-[var(--color-bg)] px-2 py-1">
              Schema {schema.contract.schema_version}
            </span>
            <span className="rounded-[4px] border border-[var(--color-warning-border)] bg-[var(--color-warning-bg)] px-2 py-1 text-[var(--color-warning)]">
              Draft only · not scientifically validated
            </span>
          </div>
        </header>
        <p className="border-b border-[var(--color-border)] py-2 text-xs text-[var(--color-text-muted)]">
          This draft stays only in this page session; a browser refresh clears this draft.
        </p>
        <ReferenceProfileSelector
          profiles={referenceProfiles}
          selectedRevisionId={draft.state.referenceProfileRevisionId}
          loading={referenceProfilesLoading}
          refreshing={referenceProfilesRefreshing}
          errorMessage={referenceProfilesError}
          selectionMessage={referenceSelectionMessage}
          onSelect={(revisionId) => {
            draft.setReferenceProfileRevisionId(revisionId);
            setReferenceSelectionMessage(null);
          }}
          onRefresh={onRefreshReferenceProfiles}
        />

        <Tabs.Root value={step} onValueChange={setStep} className="min-w-0 pt-3">
          <Tabs.List
            className="grid grid-cols-2 gap-1 border-b border-[var(--color-border)] sm:flex"
            aria-label="Input workbench sections"
          >
            {STEPS.map(({ value, label, icon: Icon }) => (
              <Tabs.Trigger
                key={value}
                value={value}
                aria-label={label}
                disabled={value === 'review' && draft.state.yamlIssue !== null}
                aria-describedby={`workbench-step-${value}-status`}
                className="inline-flex min-h-11 min-w-0 items-center justify-center gap-1.5 border-b-2 border-transparent px-3 py-2 text-sm font-medium text-[var(--color-text-muted)] outline-none transition-colors hover:text-[var(--color-text)] data-[state=active]:border-[var(--color-accent)] data-[state=active]:text-[var(--color-accent)] disabled:cursor-not-allowed disabled:opacity-50"
              >
                <Icon aria-hidden="true" size={15} />
                {label}
                <span
                  aria-hidden="true"
                  className={`h-2 w-2 shrink-0 rounded-full border ${
                    stepReady[value]
                      ? 'border-[var(--color-success)] bg-[var(--color-success)]'
                      : 'border-[var(--color-border-strong)] bg-[var(--color-surface)]'
                  }`}
                  title={stepReady[value] ? 'Ready' : 'Needs attention'}
                />
                <span id={`workbench-step-${value}-status`} className="sr-only">
                  {stepReady[value] ? 'Ready' : 'Needs attention'}
                </span>
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
              referenceSelectionAvailable={
                selectedReference !== null && referenceProfilesError === null
              }
            />
          </Tabs.Content>
        </Tabs.Root>
      </Panel>
    </section>
  );
}
