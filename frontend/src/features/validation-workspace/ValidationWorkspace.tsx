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

const inputHelp: Record<string, string> = {
  'Config (JSON)': 'JSON configuration values for this validation request.',
  'Samples (path string)': 'Path to a sample sheet on the compute environment.',
  'Options (JSON)': 'Optional validation settings such as strict_inputs.',
};

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
    <div className="space-y-3">
      <div className="text-xs text-[var(--color-text-muted)]">
        Workflow: {workflowId}
      </div>
      <div className="grid grid-cols-1 gap-3 md:grid-cols-3">
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
      </div>
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
    <div className="flex flex-col">
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
        rows={4}
        className="w-full flex-1 rounded border border-[var(--color-border)] bg-[var(--color-bg)] p-2 font-mono text-sm text-[var(--color-text)] focus:border-[var(--color-accent)] focus:outline-none focus:ring-2 focus:ring-[var(--color-accent)]"
      />
      <p className="mt-1 text-xs text-[var(--color-text-muted)]">{inputHelp[label]}</p>
    </div>
  );
}
