/**
 * One window can show several runs and local files. A shared chromosome
 * name or WashU scaffold is only a navigation list. It does not prove the
 * reference files are the same build.
 *
 * Scaffold names mirror encode_pipeline.adapters.browser.scaffolds.
 */

const SCAFFOLDS: Record<string, string> = {
  hg19: 'hg19',
  grch37: 'hg19',
  hg38: 'hg38',
  grch38: 'hg38',
  mm9: 'mm9',
  mm10: 'mm10',
  grcm38: 'mm10',
  mm39: 'mm39',
  grcm39: 'mm39',
};

export const USER_ASSEMBLIES = [
  { scaffold: 'hg38', declared: 'GRCh38' },
  { scaffold: 'hg19', declared: 'GRCh37' },
  { scaffold: 'mm10', declared: 'GRCm38' },
  { scaffold: 'mm39', declared: 'GRCm39' },
  { scaffold: 'mm9', declared: 'mm9' },
] as const;

const REGION = /^([A-Za-z0-9_.:-]+):([0-9]+)-([0-9]+)$/;

const ABSOLUTE_HTTP = /^https?:\/\//i;

export type ReferenceBinding = {
  status: string;
  identity: string | null;
  scaffold: string | null;
  assemblyLabel: string | null;
};

export type WindowAuthority = {
  scaffold: string;
  source: 'run' | 'user';
  identity: string | null;
  runId: string | null;
};

export type TrackPlacement =
  | { ok: true; warning: string | null }
  | { ok: false; reason: string };

export function viewerScaffold(assembly: string | null | undefined): string | null {
  if (typeof assembly !== 'string') return null;
  const key = assembly.trim().toLowerCase();
  if (key === '') return null;
  return SCAFFOLDS[key] ?? null;
}

export function userSpecifiedReference(scaffold: string | null | undefined): string | null {
  if (scaffold === 'hg38' || scaffold === 'mm10') return `用户指定 ${declaredAssembly(scaffold)}`;
  if (scaffold) return `User-specified ${scaffold}`;
  return null;
}

export function declaredAssembly(scaffold: string | null | undefined): string {
  return USER_ASSEMBLIES.find((item) => item.scaffold === scaffold)?.declared ?? scaffold ?? '';
}

/**
 * Colon text is 1-based and closed. WashU's ChromosomeInterval.parse keeps
 * the end and turns the start into max(0, start - 1). Region summaries use
 * that same half-open window.
 */
export function seenWindow(region: string): { chrom: string; start: number; end: number } | null {
  const parsed = REGION.exec(region);
  if (!parsed) return null;
  return {
    chrom: parsed[1],
    start: Math.max(0, Number(parsed[2]) - 1),
    end: Number(parsed[3]),
  };
}

export function referenceBinding(reference: {
  status: string;
  identity_sha256?: string | null;
  viewer_scaffold?: string | null;
  assembly_label?: string | null;
}): ReferenceBinding {
  return {
    status: reference.status,
    identity: reference.identity_sha256 ?? null,
    scaffold: reference.viewer_scaffold ?? null,
    assemblyLabel: reference.assembly_label ?? null,
  };
}

export function windowAuthority(
  primary: { runId: string; reference: ReferenceBinding } | null,
  userAssembly: string,
): WindowAuthority {
  if (primary && primary.reference.status === 'bound' && primary.reference.scaffold) {
    return {
      scaffold: primary.reference.scaffold,
      source: 'run',
      identity: primary.reference.identity,
      runId: primary.runId,
    };
  }
  return {
    scaffold: viewerScaffold(userAssembly) ?? 'hg38',
    source: 'user',
    identity: null,
    runId: null,
  };
}

export function placePlatformRun(
  authority: WindowAuthority,
  runId: string,
  reference: ReferenceBinding,
): TrackPlacement {
  if (reference.status === 'unbound') {
    return {
      ok: false,
      reason:
        'This run has no reference profile binding, so its tracks are not placed on a guessed assembly.',
    };
  }
  if (reference.status !== 'bound' || !reference.scaffold) {
    return {
      ok: false,
      reason: `Run ${runId} has reference status ${reference.status} and no WashU chromosome scaffold for ${reference.assemblyLabel ?? 'unknown'}. Its tracks stay out of this window.`,
    };
  }
  if (reference.scaffold !== authority.scaffold) {
    return {
      ok: false,
      reason: `Run ${runId} uses scaffold ${reference.scaffold}, not this window's ${authority.scaffold}. Chromosome names do not make those scaffolds the same reference.`,
    };
  }
  // identity_sha256 mixes the profile id, display name, and workflow bindings.
  // A different digest does not prove a different genome. The same assembly
  // label, including GRCh38/hg38, can share the window; the note only says
  // the profiles may be different revisions.
  if (authority.runId === runId) {
    return { ok: true, warning: null };
  }
  if (authority.identity && reference.identity && authority.identity !== reference.identity) {
    return {
      ok: true,
      warning: `Run ${runId} has its own reference profile on ${authority.scaffold}. The profiles may use different reference revisions. Chromosome names do not prove the files are the same build.`,
    };
  }
  if (authority.source === 'user') {
    return {
      ok: true,
      warning: `Run ${runId} keeps its own reference profile on ${declaredAssembly(authority.scaffold)}. This window was declared by the user; chromosome names do not prove it is that reference revision.`,
    };
  }
  return { ok: true, warning: null };
}

export function placeLocalTrack(
  authority: WindowAuthority,
  track: { label: string; assembly: string },
): TrackPlacement {
  const scaffold = viewerScaffold(track.assembly);
  const specified = userSpecifiedReference(scaffold);
  if (!scaffold) {
    return {
      ok: false,
      reason: `${track.label} declares ${track.assembly}, which this WashU build has no chromosome scaffold for. The file was not relabeled.`,
    };
  }
  if (scaffold !== authority.scaffold) {
    return {
      ok: false,
      reason: `${track.label} declares ${track.assembly}, which is not this window's ${authority.scaffold}. The file was not relabeled.`,
    };
  }
  const label = specified ?? `User-specified ${scaffold}`;
  if (authority.identity) {
    return {
      ok: true,
      warning: `${label}. Matching chromosome names do not prove this file is the run reference.`,
    };
  }
  return { ok: true, warning: label };
}

function containsForeignUrl(value: unknown, origin: string, depth: number): boolean {
  if (depth > 4) return false;
  if (typeof value === 'string') {
    if (!ABSOLUTE_HTTP.test(value)) return false;
    try {
      const url = new URL(value);
      return !(url.origin === origin && url.pathname.startsWith('/api/'));
    } catch {
      return true;
    }
  }
  if (Array.isArray(value)) {
    return value.some((item) => containsForeignUrl(item, origin, depth + 1));
  }
  if (value && typeof value === 'object') {
    return Object.values(value as Record<string, unknown>).some((item) =>
      containsForeignUrl(item, origin, depth + 1),
    );
  }
  return false;
}

export function localTrackHasForeignUrl(
  model: { [key: string]: unknown },
  origin: string,
): boolean {
  return containsForeignUrl(model, origin, 0);
}
