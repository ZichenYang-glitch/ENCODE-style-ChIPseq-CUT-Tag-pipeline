import { useCallback, useEffect, useId, useRef, useState } from 'react';
import type { ChangeEvent, KeyboardEvent } from 'react';
import {
  DEFAULT_HG38_ASSEMBLIES,
  assemblyUsesHg38,
  classifyQuery,
  formatRegion,
  loadHg38GeneIndex,
  matchLabel,
  searchGenes,
} from './geneIndex';
import type { GeneSearch, ParsedGeneIndex, QueryPlan } from './geneIndex';

export interface LocusSearchProps {
  assembly: string;
  region: string;
  onNavigate: (region: string) => void;
}

type Phase = 'loading' | 'ready' | 'failed' | 'unsupported';

function statusMessage(
  phase: Phase,
  assembly: string,
  plan: QueryPlan,
  result: GeneSearch | null,
  release: string,
  selectable: number,
  unavailable: string[],
): string | null {
  if (phase === 'failed') {
    return 'Gene annotation did not load. Coordinates still work.';
  }
  if (phase === 'unsupported') {
    const name = assembly.trim() || 'this assembly';
    return `No bundled gene annotation for ${name}. Coordinates still work.`;
  }
  if (plan.kind === 'empty' || plan.kind === 'coordinate') return null;
  if (plan.kind === 'coordinate-draft') return 'Coordinates use chrom:start-end.';
  if (phase === 'loading') return 'Loading hg38 gene annotation…';
  if (!result || result.relation === 'none') {
    return `No genes match "${plan.query}" in ${release}. Coordinates still work.`;
  }
  const reason = unavailable.length ? `WashU hg38 cannot display ${unavailable.join(', ')}.` : '';
  if (!selectable && reason) return `${plan.query} has no displayable location. ${reason}`;
  if (result.relation === 'prefix') {
    return `No exact gene, alias, or gene id for "${plan.query}". Showing ${selectable} of ${result.total} matching locations. ${reason}`.trimEnd();
  }
  const noun = result.total === 1 ? 'location' : 'locations';
  return `${plan.query} matches ${result.total} ${noun}. Choose one. ${reason}`.trimEnd();
}

export default function LocusSearch({ assembly, region, onNavigate }: LocusSearchProps) {
  const listId = useId();
  const helpId = useId();
  const pendingEnter = useRef<{ assembly: string; query: string } | null>(null);
  const [draft, setDraft] = useState(() => region.trim());
  const seenRegion = useRef(region);
  if (seenRegion.current !== region) {
    seenRegion.current = region;
    pendingEnter.current = null;
    setDraft(region.trim());
  }
  const [pinnedDraft, setPinnedDraft] = useState(draft);
  const [activeIndex, setActiveIndex] = useState(-1);
  const [dismissedQuery, setDismissedQuery] = useState<string | null>(null);
  if (pinnedDraft !== draft) {
    setPinnedDraft(draft);
    setActiveIndex(-1);
    setDismissedQuery(null);
  }
  const [phase, setPhase] = useState<Phase>(() =>
    assemblyUsesHg38(assembly) ? 'loading' : 'unsupported',
  );
  const [index, setIndex] = useState<ParsedGeneIndex | null>(null);

  useEffect(() => {
    if (!assemblyUsesHg38(assembly, DEFAULT_HG38_ASSEMBLIES)) {
      setIndex(null);
      setPhase('unsupported');
      return;
    }
    let cancelled = false;
    setPhase((current) => (current === 'ready' ? current : 'loading'));
    loadHg38GeneIndex()
      .then((loaded) => {
        if (cancelled) return;
        if (!assemblyUsesHg38(assembly, loaded.assemblyMatch)) {
          setIndex(null);
          setPhase('unsupported');
          return;
        }
        setIndex(loaded);
        setPhase('ready');
      })
      .catch(() => {
        if (cancelled) return;
        setIndex(null);
        setPhase('failed');
      });
    return () => {
      cancelled = true;
    };
  }, [assembly]);

  const plan = classifyQuery(draft);
  const result =
    phase === 'ready' &&
    index &&
    assemblyUsesHg38(assembly, index.assemblyMatch) &&
    plan.kind === 'gene'
      ? searchGenes(index, plan.query)
      : null;
  const candidates = result?.matches ?? [];
  const matches = candidates.filter(({ locus }) => index?.viewerChromosomes.has(locus.chrom));
  const excludedChromosomes = candidates
    .map(({ locus }) => locus.chrom)
    .filter((chrom) => !index?.viewerChromosomes.has(chrom));
  const open = matches.length > 0 && dismissedQuery !== draft;
  const active = activeIndex < matches.length ? activeIndex : -1;
  const status = statusMessage(
    phase, assembly, plan, result, index?.release ?? 'this hg38 annotation',
    matches.length, [...new Set(excludedChromosomes)],
  );
  const soleMatch = result?.relation === 'exact' && result.total === 1 ? matches[0] : undefined;
  const soleRegion = soleMatch && formatRegion(soleMatch.locus);

  const navigate = useCallback((next: string) => {
    pendingEnter.current = null;
    setDraft(next);
    onNavigate(next);
  }, [onNavigate]);

  useEffect(() => {
    const pending = pendingEnter.current;
    if (!pending) return;
    if (pending.assembly !== assembly || pending.query !== draft) {
      pendingEnter.current = null;
      return;
    }
    if (phase === 'loading') return;
    pendingEnter.current = null;
    if (phase === 'ready' && soleRegion) navigate(soleRegion);
  }, [assembly, draft, phase, soleRegion, navigate]);

  function onChange(event: ChangeEvent<HTMLInputElement>) {
    pendingEnter.current = null;
    const next = event.target.value.trim();
    setDraft(next);
    const nextPlan = classifyQuery(next);
    if (nextPlan.kind === 'coordinate') onNavigate(nextPlan.region);
  }

  function onKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    if (event.key === 'Escape') {
      pendingEnter.current = null;
      setDismissedQuery(draft);
      setActiveIndex(-1);
      return;
    }
    if (event.key !== 'Enter' && event.key !== 'ArrowDown' && event.key !== 'ArrowUp') return;
    if (plan.kind !== 'gene') return;
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      pendingEnter.current = null;
      event.preventDefault();
      if (matches.length === 0) return;
      setDismissedQuery(null);
      const last = matches.length - 1;
      setActiveIndex((current) => {
        const bounded = Math.min(current, last);
        if (event.key === 'ArrowDown') return bounded >= last ? 0 : bounded + 1;
        return bounded <= 0 ? last : bounded - 1;
      });
      return;
    }
    event.preventDefault();
    if (dismissedQuery === draft) {
      setDismissedQuery(null);
      return;
    }
    if (phase === 'loading') {
      pendingEnter.current = { assembly, query: draft };
      return;
    }
    const match = matches[active] ?? soleMatch;
    if (match) navigate(formatRegion(match.locus));
  }

  return (
    <div className="block text-sm">
      <label className="block">
        Genomic region
        <input
          className="mt-1 w-full max-w-md border border-[var(--color-border)] px-2 py-1"
          value={draft}
          onChange={onChange}
          onKeyDown={onKeyDown}
          aria-label="Genomic region"
          aria-describedby={helpId}
          aria-autocomplete="list"
          aria-expanded={open}
          aria-controls={open ? listId : undefined}
          aria-activedescendant={open && active >= 0 ? `${listId}-option-${active}` : undefined}
          spellCheck={false}
          autoCapitalize="off"
          autoCorrect="off"
          autoComplete="off"
        />
      </label>
      <p id={helpId} className="mt-1 max-w-md text-[var(--color-text-muted)]">
        Enter a gene name (for example EGFR) or chrom:start-end. For one exact match, press Enter.
        To choose a location, tap or click it, or use the arrow keys and press Enter.
      </p>
      {status ? (
        <p role="status" className="mt-1 max-w-md text-[var(--color-text-muted)]">
          {status}
        </p>
      ) : null}
      {open ? (
        <ul
          id={listId}
          role="listbox"
          aria-label="Gene matches"
          className="mt-1 max-h-60 w-full max-w-md overflow-auto border border-[var(--color-border)]"
        >
          {matches.map((match, option) => {
            const label = matchLabel(match);
            const selected = option === active;
            return (
              <li
                key={`${match.locus.geneId}:${match.locus.chrom}:${match.locus.start}:${match.locus.end}`}
                id={`${listId}-option-${option}`}
                role="option"
                aria-selected={selected}
                data-region={formatRegion(match.locus)}
                data-symbol={match.locus.symbol}
                data-gene-id={match.locus.geneId}
                className={
                  selected
                    ? 'cursor-pointer bg-[var(--color-border)] px-2 py-1'
                    : 'cursor-pointer px-2 py-1 hover:bg-[var(--color-border)]'
                }
                onMouseDown={(event) => event.preventDefault()}
                onClick={() => navigate(formatRegion(match.locus))}
              >
                {label}
              </li>
            );
          })}
        </ul>
      ) : null}
    </div>
  );
}
