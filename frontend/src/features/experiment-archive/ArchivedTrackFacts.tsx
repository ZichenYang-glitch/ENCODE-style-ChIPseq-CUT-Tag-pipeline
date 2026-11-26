import type { UiLanguage } from '../../lib/uiLanguage';
import { archiveText } from './archiveText';
import type { ArchivedTrack } from '../../api/archivedTracks';
import { archivedSourceText, archivedVersionText } from './archiveView';

export function ArchivedTrackFacts({ track, language }: { track: ArchivedTrack; language: UiLanguage }) {
  const text = archiveText[language];
  return (
    <dl className="grid gap-1 text-sm">
      <div>
        <dt className="text-[var(--color-text-muted)]">{text.name}</dt>
        <dd data-testid="archived-track-label">{track.display_name}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">{text.version}</dt>
        <dd data-testid="archived-track-version-label">{archivedVersionText(track.version_number, language)}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">{text.source}</dt>
        <dd data-testid="archived-track-source">{archivedSourceText(track, language)}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">{text.assembly}</dt>
        <dd data-testid="archived-track-assembly">{track.assembly}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">{text.experimentRevision}</dt>
        <dd data-testid="archived-track-experiment-revision">{track.experiment_revision_id}</dd>
      </div>
      <div>
        <dt className="text-[var(--color-text-muted)]">{text.key}</dt>
        <dd data-testid="archived-track-key">{track.track_key}</dd>
      </div>
    </dl>
  );
}
