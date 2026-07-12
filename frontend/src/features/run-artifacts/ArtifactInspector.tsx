import { useEffect, useState } from 'react';
import { Copy } from 'lucide-react';
import type { ArtifactReferenceResponse } from '../../api/generated/models';
import { Button } from '../../components/Button';
import { formatBytes, formatProducedTime } from './artifactState';

interface ArtifactInspectorProps {
  artifact: ArtifactReferenceResponse | null;
  selectedArtifactId: string | null;
  isLoading: boolean;
  isError: boolean;
  invalidSelection: boolean;
  onRetry: () => void;
}

const metadataLabels: Array<[
  keyof ArtifactReferenceResponse['metadata'],
  string,
]> = [
  ['catalog_id', 'Catalog ID'],
  ['scope', 'Scope'],
  ['level', 'Level'],
  ['sample_id', 'Sample'],
  ['experiment_id', 'Experiment'],
  ['assay', 'Assay'],
  ['target', 'Target'],
  ['genome', 'Genome'],
  ['method', 'Method'],
  ['qc_flag', 'QC flag'],
];

export function ArtifactInspector({
  artifact,
  selectedArtifactId,
  isLoading,
  isError,
  invalidSelection,
  onRetry,
}: ArtifactInspectorProps) {
  const [copyStatus, setCopyStatus] = useState<string | null>(null);

  useEffect(() => setCopyStatus(null), [selectedArtifactId]);

  async function copyValue(label: string, value: string) {
    try {
      if (!navigator.clipboard?.writeText) throw new Error('clipboard unavailable');
      await navigator.clipboard.writeText(value);
      setCopyStatus(`${label} copied.`);
    } catch {
      setCopyStatus(`Could not copy ${label.toLowerCase()}.`);
    }
  }

  return (
    <aside
      className="min-w-0 border-t border-[var(--color-border)] pt-3 xl:border-l xl:border-t-0 xl:pl-4 xl:pt-0"
      aria-labelledby="artifact-inspector-heading"
      data-testid="artifact-inspector"
    >
      <h3 id="artifact-inspector-heading" className="text-sm font-semibold">
        Artifact details
      </h3>

      {!selectedArtifactId && !invalidSelection && (
        <p className="mt-2 text-sm text-[var(--color-text-muted)]">
          Select an artifact to inspect its persisted metadata.
        </p>
      )}

      {invalidSelection && (
        <div className="mt-2" role="status">
          <p className="text-sm text-[var(--color-error)]">
            This artifact link is invalid or unavailable.
          </p>
          <p className="mt-1 text-xs text-[var(--color-text-muted)]">
            Return to the artifact list and choose a persisted item.
          </p>
        </div>
      )}

      {isLoading && (
        <div className="mt-3 min-h-48 animate-pulse space-y-3" role="status" aria-label="Loading artifact details">
          <div className="h-4 w-2/3 rounded bg-gray-200" />
          <div className="h-12 rounded bg-gray-100" />
          <div className="h-12 rounded bg-gray-100" />
        </div>
      )}

      {isError && !isLoading && (
        <div className="mt-3" role="status">
          <p className="text-sm text-[var(--color-error)]">
            Artifact details could not be loaded.
          </p>
          <Button className="mt-2" variant="secondary" onClick={onRetry}>
            Retry details
          </Button>
        </div>
      )}

      {artifact && !isLoading && !isError && (
        <div className="mt-3 min-w-0 space-y-4">
          <div className="min-w-0">
            <p className="break-words text-sm font-medium">{artifact.output_type}</p>
            <p className="break-all text-xs text-[var(--color-text-muted)]" title={artifact.name}>
              {artifact.name}
            </p>
          </div>

          <dl className="grid min-w-0 grid-cols-[6rem_minmax(0,1fr)] gap-x-3 gap-y-2 text-xs">
            <dt className="text-[var(--color-text-muted)]">Size</dt>
            <dd>{formatBytes(artifact.size_bytes)}</dd>
            <dt className="text-[var(--color-text-muted)]">Produced</dt>
            <dd><time dateTime={artifact.produced_at}>{formatProducedTime(artifact.produced_at)}</time></dd>
            <dt className="text-[var(--color-text-muted)]">MIME type</dt>
            <dd className="break-all">{artifact.mime_type ?? '—'}</dd>
            {metadataLabels.map(([key, label]) => {
              const value = artifact.metadata[key];
              if (!value) return null;
              return (
                <div className="contents" key={key}>
                  <dt className="text-[var(--color-text-muted)]">{label}</dt>
                  <dd className="break-all">{value}</dd>
                </div>
              );
            })}
          </dl>

          <CopyField
            label="Relative path"
            value={artifact.relative_path}
            onCopy={copyValue}
          />
          <CopyField label="Opaque URI" value={artifact.uri} onCopy={copyValue} />
        </div>
      )}

      <p className="sr-only" role="status" aria-live="polite">
        {copyStatus ?? ''}
      </p>
      {copyStatus && (
        <p className="mt-3 text-xs text-[var(--color-text-muted)]" aria-hidden="true">
          {copyStatus}
        </p>
      )}
    </aside>
  );
}

function CopyField({
  label,
  value,
  onCopy,
}: {
  label: string;
  value: string;
  onCopy: (label: string, value: string) => Promise<void>;
}) {
  return (
    <div className="min-w-0 border-t border-[var(--color-border)] pt-3">
      <div className="flex min-w-0 items-start justify-between gap-2">
        <div className="min-w-0">
          <p className="text-xs font-medium text-[var(--color-text-muted)]">{label}</p>
          <code className="mt-1 block break-all text-xs" title={value}>{value}</code>
        </div>
        <Button
          variant="secondary"
          className="shrink-0 px-2"
          onClick={() => void onCopy(label, value)}
          aria-label={`Copy ${label.toLowerCase()}`}
          title={`Copy ${label.toLowerCase()}`}
        >
          <Copy className="h-4 w-4" aria-hidden="true" />
        </Button>
      </div>
    </div>
  );
}

