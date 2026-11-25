export function archivedTrackCollectionPath(projectId: string, experimentId: string): string {
  return `/api/v1/catalog/projects/${projectId}/experiments/${experimentId}/archived-tracks`;
}

export function archivedTrackResourcePath(
  projectId: string,
  experimentId: string,
  archiveVersionId: string,
): string {
  return `${archivedTrackCollectionPath(projectId, experimentId)}/${archiveVersionId}`;
}
