import type { ArchivedTrack } from '../../api/archivedTracks';

export const PROJECT_ID = `prj_${'1'.repeat(32)}`;
export const EXPERIMENT_ID = `exp_${'4'.repeat(32)}`;
export const REVISION_V1 = `expr_${'5'.repeat(32)}`;
export const REVISION_V2 = `expr_${'6'.repeat(32)}`;
export const ARCHIVE_V1 = `arv_${'a'.repeat(32)}`;
export const ARCHIVE_V2 = `arv_${'b'.repeat(32)}`;

export function sealedTrack(
  version: 1 | 2,
  overrides: Partial<ArchivedTrack> = {},
): ArchivedTrack {
  const archiveVersionId = version === 1 ? ARCHIVE_V1 : ARCHIVE_V2;
  const track: ArchivedTrack = {
    archive_version_id: archiveVersionId,
    project_id: PROJECT_ID,
    experiment_id: EXPERIMENT_ID,
    experiment_revision_id: version === 1 ? REVISION_V1 : REVISION_V2,
    track_key: 'signal.bw',
    version_number: version,
    display_name: 'signal.bw',
    assembly: 'hg38',
    track_type: 'bigwig',
    has_index: false,
    source: {
      kind: 'run',
      run_id: version === 1 ? 'run-v1' : 'run-v2',
      artifact_id: 'signal',
      generation: `artifactgen-${version === 1 ? 'a' : 'c'}${'a'.repeat(63)}`,
      revision: `artifactrev-${version === 1 ? 'b' : 'd'}${'b'.repeat(63)}`,
    },
  };
  return {
    ...track,
    ...overrides,
    archive_version_id: overrides.archive_version_id ?? track.archive_version_id,
  };
}
