import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  safeArtifactDownloadFilename,
  saveArtifactBlob,
} from './artifactDownload';

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('artifact download browser boundary', () => {
  it('keeps safe persisted names and replaces header-injection or path names', () => {
    expect(
      safeArtifactDownloadFilename('result manifest.tsv', 'artifact-a'),
    ).toBe('result manifest.tsv');
    for (const candidate of [
      '../private.tsv',
      'private\\source.tsv',
      'source\r\nX-Private: yes',
      '.',
      '..',
    ]) {
      expect(safeArtifactDownloadFilename(candidate, 'artifact-a')).toBe(
        'artifact-a.download',
      );
    }
  });

  it('always revokes the object URL and removes its temporary anchor', () => {
    const createObjectURL = vi.fn().mockReturnValue('blob:artifact');
    const revokeObjectURL = vi.fn();
    vi.stubGlobal('URL', { ...URL, createObjectURL, revokeObjectURL });
    const click = vi
      .spyOn(HTMLAnchorElement.prototype, 'click')
      .mockImplementation(() => undefined);

    saveArtifactBlob(new Blob(['artifact']), 'artifact.tsv');

    expect(click).toHaveBeenCalledTimes(1);
    expect(revokeObjectURL).toHaveBeenCalledWith('blob:artifact');
    expect(document.querySelector('a[download="artifact.tsv"]')).toBeNull();
  });
});
