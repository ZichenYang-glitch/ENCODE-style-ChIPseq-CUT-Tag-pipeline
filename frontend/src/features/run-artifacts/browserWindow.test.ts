import { describe, expect, it } from 'vitest';
import {
  declaredAssembly,
  localTrackHasForeignUrl,
  placeLocalTrack,
  placePlatformRun,
  seenWindow,
  userSpecifiedReference,
  viewerScaffold,
  windowAuthority,
  type ReferenceBinding,
} from './browserWindow';

function binding(overrides: Partial<ReferenceBinding> = {}): ReferenceBinding {
  return {
    status: 'bound',
    identity: 'a'.repeat(64),
    scaffold: 'hg38',
    assemblyLabel: 'GRCh38',
    ...overrides,
  };
}

describe('browser window authority', () => {
  it('maps assembly labels onto WashU scaffolds without treating the names as identity', () => {
    expect(viewerScaffold('GRCh38')).toBe('hg38');
    expect(viewerScaffold('hg38')).toBe('hg38');
    expect(viewerScaffold('mm10')).toBe('mm10');
    expect(viewerScaffold('not-a-genome')).toBeNull();
    expect(userSpecifiedReference('hg38')).toBe('用户指定 GRCh38');
    expect(declaredAssembly('hg38')).toBe('GRCh38');
    expect(seenWindow('chr7:100-200')).toEqual({ chrom: 'chr7', start: 99, end: 200 });
    expect(seenWindow('chr1:1-2')).toEqual({ chrom: 'chr1', start: 0, end: 2 });
  });

  it('keeps two runs on one assembly when their public identities differ', () => {
    const authority = windowAuthority({ runId: 'run-1', reference: binding() }, 'mm10');
    expect(authority).toMatchObject({ scaffold: 'hg38', identity: 'a'.repeat(64) });
    const differentIdentity = placePlatformRun(
      authority,
      'run-2',
      binding({ identity: 'b'.repeat(64), assemblyLabel: 'GRCh38' }),
    );
    expect(differentIdentity.ok).toBe(true);
    if (differentIdentity.ok) {
      expect(differentIdentity.warning).toMatch(/different reference revisions/i);
      expect(differentIdentity.warning).toMatch(/do not prove/i);
    }
    expect(placePlatformRun(authority, 'run-1', binding())).toMatchObject({ ok: true, warning: null });
  });

  it('keeps an unbound run off a guessed assembly and rejects a different scaffold', () => {
    const authority = windowAuthority({ runId: 'run-1', reference: binding() }, 'mm10');
    expect(
      placePlatformRun(authority, 'run-2', binding({ status: 'unbound', identity: null, scaffold: null })).ok,
    ).toBe(false);
    const otherScaffold = placePlatformRun(
      authority,
      'run-2',
      binding({ scaffold: 'mm10', assemblyLabel: 'mm10' }),
    );
    expect(otherScaffold.ok).toBe(false);
    if (!otherScaffold.ok) expect(otherScaffold.reason).toMatch(/not this window's hg38/i);
  });

  it('allows every same-assembly run on a user-declared window and blocks a different assembly', () => {
    const authority = windowAuthority(null, 'hg38');
    expect(authority).toMatchObject({ source: 'user', scaffold: 'hg38', identity: null });
    expect(placePlatformRun(authority, 'run-a', binding({ identity: 'c'.repeat(64) })).ok).toBe(true);
    expect(placePlatformRun(authority, 'run-b', binding({ identity: 'd'.repeat(64) })).ok).toBe(true);
    expect(
      placePlatformRun(authority, 'run-c', binding({ scaffold: 'mm10', assemblyLabel: 'mm10' })).ok,
    ).toBe(false);
  });

  it('does not relabel a local file onto another scaffold', () => {
    const authority = windowAuthority({ runId: 'run-1', reference: binding() }, 'hg38');
    const mismatch = placeLocalTrack(authority, { label: 'Mouse signal', assembly: 'mm10' });
    expect(mismatch.ok).toBe(false);
    if (!mismatch.ok) expect(mismatch.reason).toMatch(/not relabeled/i);
    const placed = placeLocalTrack(authority, { label: 'Local signal', assembly: 'GRCh38' });
    expect(placed.ok).toBe(true);
    if (placed.ok) {
      expect(placed.warning ?? '').toMatch(/用户指定 GRCh38/);
      expect(placed.warning ?? '').toMatch(/do not prove/i);
    }
  });

  it('rejects a local model that points at another website', () => {
    expect(
      localTrackHasForeignUrl(
        { type: 'bigwig', name: 'Foreign', url: 'https://epigenome.wustl.edu/proxy/signal.bw' },
        'https://helix.example',
      ),
    ).toBe(true);
    expect(
      localTrackHasForeignUrl(
        { type: 'bigwig', name: 'Local', url: 'blob:https://helix.example/1111' },
        'https://helix.example',
      ),
    ).toBe(false);
    expect(
      localTrackHasForeignUrl(
        {
          type: 'bigwig',
          name: 'Api',
          options: { url: 'https://helix.example/api/v1/runs/run-1/browser-tracks/a/bytes/g/r' },
        },
        'https://helix.example',
      ),
    ).toBe(false);
  });
});
