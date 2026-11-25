import type { ArchivedTrack } from '../../api/archivedTracks';
import { archivedSourceText, archivedVersionText } from './archiveView';

export function ArchivedTrackFacts({ track }: { track: ArchivedTrack }) {
  return (
    <dl className="grid gap-1 text-sm">
      <div>
        <dt className="text-[var(--color-text-muted)]">名称</dt>
        <dd data-testid="archived-track-label">{track.display_name}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">版本</dt>
        <dd data-testid="archived-track-version-label">{archivedVersionText(track.version_number)}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">来源</dt>
        <dd data-testid="archived-track-source">{archivedSourceText(track)}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">组装</dt>
        <dd data-testid="archived-track-assembly">{track.assembly}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">实验修订</dt>
        <dd data-testid="archived-track-experiment-revision">{track.experiment_revision_id}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">轨道</dt>
        <dd data-testid="archived-track-key">{track.track_key}</dd>
      </div>
    </dl>
  );
}
