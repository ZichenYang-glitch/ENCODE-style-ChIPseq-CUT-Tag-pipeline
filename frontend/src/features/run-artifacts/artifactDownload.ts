const UNSAFE_FILENAME = /[\\/\u0000-\u001f\u007f]/;

export function safeArtifactDownloadFilename(
  candidate: string,
  artifactId: string,
): string {
  if (
    candidate.length > 0 &&
    candidate.length <= 255 &&
    candidate !== '.' &&
    candidate !== '..' &&
    !UNSAFE_FILENAME.test(candidate)
  ) {
    return candidate;
  }
  return `${artifactId}.download`;
}

export function saveArtifactBlob(blob: Blob, filename: string): void {
  const objectUrl = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = objectUrl;
  anchor.download = filename;
  anchor.rel = 'noopener';
  anchor.hidden = true;
  document.body.appendChild(anchor);
  try {
    anchor.click();
  } finally {
    anchor.remove();
    URL.revokeObjectURL(objectUrl);
  }
}
