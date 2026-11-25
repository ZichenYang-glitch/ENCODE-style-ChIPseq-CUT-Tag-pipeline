import { useEffect, useId, useRef, useState } from 'react';
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
import type { GeneMatch, GeneSearch, ParsedGeneIndex, QueryPlan } from './geneIndex';

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
  if (result.relation === 'prefix') {
    return `No exact gene, alias, or gene id for "${plan.query}". Showing ${result.matches.length} of ${result.total} names that start with it.`;
  }
  const noun = result.total === 1 ? 'location' : 'locations';
  return `${plan.query} matches ${result.total} ${noun}. Choose one.`;
}

export default function LocusSearch({ assembly, region, onNavigate }: LocusSearchProps) {
  const listId = useId();
  const [draft, setDraft] = useState(() => region.trim());
  const seenRegion = useRef(region);
  if (seenRegion.current !== region) {
    seenRegion.current = region;
    setDraft(region.trim());
  }
  const [pinnedDraft, setPinnedDraft] = useState(draft);
  const [activeIndex, setActiveIndex] = useState(0);
  const [dismissedQuery, setDismissedQuery] = useState<string | null>(null);
  if (pinnedDraft !== draft) {
    setPinnedDraft(draft);
    setActiveIndex(0);
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
    phase === 'ready' && index && plan.kind === 'gene' ? searchGenes(index, plan.query) : null;
  const matches = result?.matches ?? [];
  const open = matches.length > 0 && dismissedQuery !== draft;
  const active = matches.length > 0 ? Math.min(activeIndex, matches.length - 1) : 0;
  const status = statusMessage(phase, assembly, plan, result, index?.release ?? 'this hg38 annotation');

  function choose(match: GeneMatch) {
    const next = formatRegion(match.locus);
    setDraft(next);
    onNavigate(next);
  }

  function onChange(event: ChangeEvent<HTMLInputElement>) {
    const next = event.target.value.trim();
    setDraft(next);
    const nextPlan = classifyQuery(next);
    if (nextPlan.kind === 'coordinate') onNavigate(nextPlan.region);
  }

  function onKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    if (event.key === 'Escape') {
      setDismissedQuery(draft);
      return;
    }
    if (event.key !== 'Enter' && event.key !== 'ArrowDown' && event.key !== 'ArrowUp') return;
    if (classifyQuery(draft).kind !== 'gene') return;
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
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
    const match = matches[active];
    if (match) choose(match);
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
          aria-autocomplete="list"
          aria-expanded={open}
          aria-controls={open ? listId : undefined}
          aria-activedescendant={open ? `${listId}-option-${active}` : undefined}
          spellCheck={false}
          autoCapitalize="off"
          autoCorrect="off"
          autoComplete="off"
        />
      </label>
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
                    : 'cursor-pointer px-2 py-1'
                }
                onMouseEnter={() => setActiveIndex(option)}
                onMouseDown={(event) => {
                  event.preventDefault();
                  choose(match);
                }}
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
