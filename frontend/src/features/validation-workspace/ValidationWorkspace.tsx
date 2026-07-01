import { Button } from '../../components/Button';

interface ValidationWorkspaceProps {
  workflowId: string;
  configText: string;
  samplesText: string;
  optionsText: string;
  loading: boolean;
  onConfigChange: (value: string) => void;
  onSamplesChange: (value: string) => void;
  onOptionsChange: (value: string) => void;
  onValidate: () => void;
}

export function ValidationWorkspace({
  workflowId,
  configText,
  samplesText,
  optionsText,
  loading,
  onConfigChange,
  onSamplesChange,
  onOptionsChange,
  onValidate,
}: ValidationWorkspaceProps) {
  return (
    <div className="space-y-4">
      <div className="text-sm text-[var(--color-text-muted)]">
        Workflow: {workflowId}
      </div>
      <InputSection
        label="Config (JSON)"
        value={configText}
        onChange={onConfigChange}
      />
      <InputSection
        label="Samples (path string)"
        value={samplesText}
        onChange={onSamplesChange}
      />
      <InputSection
        label="Options (JSON)"
        value={optionsText}
        onChange={onOptionsChange}
      />
      <Button
        variant="primary"
        onClick={onValidate}
        disabled={loading}
        data-testid="validate-button"
      >
        {loading ? 'Validating…' : 'Validate'}
      </Button>
    </div>
  );
}

function InputSection({
  label,
  value,
  onChange,
}: {
  label: string;
  value: string;
  onChange: (value: string) => void;
}) {
  const id = label.toLowerCase().replace(/[^a-z]+/g, '-');
  return (
    <div>
      <label
        htmlFor={id}
        className="mb-1 block text-xs font-semibold uppercase tracking-wide text-[var(--color-text-muted)]"
      >
        {label}
      </label>
      <textarea
        id={id}
        value={value}
        onChange={(e) => onChange(e.target.value)}
        rows={5}
        className="w-full rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-2 font-mono text-sm text-[var(--color-text)] focus:border-[var(--color-accent)] focus:outline-none"
      />
    </div>
  );
}
