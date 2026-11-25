import { useEffect, useId, useRef, useState } from 'react';
import type { ChangeEvent, KeyboardEvent } from 'react';
import {
  DEFAULT_HG38_ASSEMBLIES,
  assemblyUsesHg38,
  classifyQuery,
  coordinateSequence,
  formatRegion,
  loadHg38GeneIndex,
  matchLabel,
  missingPinnedSequence,
  searchGenes,
} from './geneIndex';
import type { GeneMatch, GeneSearch, ParsedGeneIndex, QueryPlan } from './geneIndex';

export interface LocusSearchProps {
  assembly: string;
  region: string;
  onNavigate: (region: string) => void;
}

type Phase = 'loading' | 'ready' | 'failed' | 'unsupported';

const HINT =
  '基因名（如 EGFR）或 chr:start-end。多个位置请点选，或用方向键选中后再按 Enter；唯一结果可直接按 Enter。';

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
  if (result.total === 1) return `${plan.query} matches 1 location. Press Enter or choose it.`;
  return `${plan.query} matches ${result.total} locations. Choose one. Arrow keys, then Enter.`;
}

export default function LocusSearch({ assembly, region, onNavigate }: LocusSearchProps) {
  const listId = useId();
  const hintId = useId();
  const pendingEnter = useRef<string | null>(null);
  const [draft, setDraft] = useState(() => region.trim());
  const [blocked, setBlocked] = useState<string | null>(null);
  const seenRegion = useRef(region);
  if (seenRegion.current !== region) {
    seenRegion.current = region;
    const next = region.trim();
    setDraft(next);
    if (pendingEnter.current && pendingEnter.current !== next) pendingEnter.current = null;
    setBlocked(null);
  }
  const [pinnedDraft, setPinnedDraft] = useState(draft);
  const [armedIndex, setArmedIndex] = useState<number | null>(null);
  const [dismissedQuery, setDismissedQuery] = useState<string | null>(null);
  if (pinnedDraft !== draft) {
    setPinnedDraft(draft);
    setArmedIndex(null);
    setDismissedQuery(null);
  }
  const [phase, setPhase] = useState<Phase>(() =>
    assemblyUsesHg38(assembly) ? 'loading' : 'unsupported',
  );
  const [index, setIndex] = useState<ParsedGeneIndex | null>(null);

  useEffect(() => {
    pendingEnter.current = null;
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

  function deliver(match: GeneMatch) {
    const absent = assemblyUsesHg38(assembly) ? missingPinnedSequence(match.locus.chrom) : null;
    if (absent) {
      setBlocked(absent);
      return;
    }
    const next = formatRegion(match.locus);
    setBlocked(null);
    setDraft(next);
    onNavigate(next);
  }

  useEffect(() => {
    if (phase !== 'ready' || !index) return;
    const query = pendingEnter.current;
    if (!query || query !== draft) return;
    pendingEnter.current = null;
    const plan = classifyQuery(query);
    if (plan.kind !== 'gene') return;
    const found = searchGenes(index, plan.query);
    if (found.relation !== 'exact' || found.matches.length !== 1) return;
    const match = found.matches[0];
    if (!match) return;
    const absent = assemblyUsesHg38(assembly) ? missingPinnedSequence(match.locus.chrom) : null;
    if (absent) {
      setBlocked(absent);
      return;
    }
    const next = formatRegion(match.locus);
    setBlocked(null);
    setDraft(next);
    onNavigate(next);
  }, [phase, index, draft, assembly, onNavigate]);

  const plan = classifyQuery(draft);
  const result =
    phase === 'ready' && index && plan.kind === 'gene' ? searchGenes(index, plan.query) : null;
  const matches = result?.matches ?? [];
  const open = matches.length > 0 && dismissedQuery !== draft;
  const active =
    armedIndex !== null && matches.length > 0 ? Math.min(armedIndex, matches.length - 1) : null;
  const status =
    blocked ??
    statusMessage(phase, assembly, plan, result, index?.release ?? 'this hg38 annotation');

  function onChange(event: ChangeEvent<HTMLInputElement>) {
    const next = event.target.value.trim();
    setDraft(next);
    if (pendingEnter.current && pendingEnter.current !== next) pendingEnter.current = null;
    const nextPlan = classifyQuery(next);
    if (nextPlan.kind !== 'coordinate') {
      setBlocked(null);
      return;
    }
    pendingEnter.current = null;
    const sequence = coordinateSequence(nextPlan.region);
    const absent =
      assemblyUsesHg38(assembly) && sequence ? missingPinnedSequence(sequence) : null;
    if (absent) {
      setBlocked(absent);
      return;
    }
    setBlocked(null);
    onNavigate(nextPlan.region);
  }

  function onKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    if (event.key === 'Escape') {
      pendingEnter.current = null;
      setDismissedQuery(draft);
      setArmedIndex(null);
      return;
    }
    if (event.key !== 'Enter' && event.key !== 'ArrowDown' && event.key !== 'ArrowUp') return;
    if (classifyQuery(draft).kind !== 'gene') return;
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      pendingEnter.current = null;
      event.preventDefault();
      if (matches.length === 0) return;
      setDismissedQuery(null);
      const last = matches.length - 1;
      setArmedIndex((current) => {
        if (current === null) return event.key === 'ArrowDown' ? 0 : last;
        const bounded = Math.min(current, last);
        if (event.key === 'ArrowDown') return bounded >= last ? 0 : bounded + 1;
        return bounded <= 0 ? last : bounded - 1;
      });
      return;
    }
    event.preventDefault();
    if (phase === 'loading') {
      pendingEnter.current = draft;
      return;
    }
    if (phase !== 'ready' || !index) return;
    const found = searchGenes(index, draft);
    if (found.relation === 'exact' && found.matches.length === 1) {
      const match = found.matches[0];
      if (match) deliver(match);
      return;
    }
    if (dismissedQuery === draft) {
      setDismissedQuery(null);
      return;
    }
    if (armedIndex === null) return;
    const match = found.matches[Math.min(armedIndex, found.matches.length - 1)];
    if (match) deliver(match);
  }

  return (
    <div className="block text-sm">
      <label className="block">
        <span id={hintId} className="block max-w-md">{HINT}</span>
        <input
          className="mt-1 w-full max-w-md border border-[var(--color-border)] px-2 py-1"
          value={draft}
          onChange={onChange}
          onKeyDown={onKeyDown}
          aria-label="Genomic region"
          aria-describedby={hintId}
          aria-autocomplete="list"
          aria-expanded={open}
          aria-controls={open ? listId : undefined}
          aria-activedescendant={open && active !== null ? `${listId}-option-${active}` : undefined}
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
        <div
          id={listId}
          role="listbox"
          aria-label="Gene matches"
          className="mt-1 max-h-60 w-full max-w-md overflow-auto border border-[var(--color-border)]"
        >
          {matches.map((match, option) => {
            const absent = assemblyUsesHg38(assembly)
              ? missingPinnedSequence(match.locus.chrom)
              : null;
            const label = absent
              ? `${matchLabel(match)} · not in pinned hg38`
              : matchLabel(match);
            const selected = option === active;
            return (
              <button
                key={`${match.locus.geneId}:${match.locus.chrom}:${match.locus.start}:${match.locus.end}`}
                id={`${listId}-option-${option}`}
                type="button"
                role="option"
                aria-selected={selected}
                data-region={formatRegion(match.locus)}
                data-symbol={match.locus.symbol}
                data-gene-id={match.locus.geneId}
                data-navigable={absent ? 'false' : 'true'}
                className={
                  selected
                    ? 'block w-full cursor-pointer bg-[var(--color-border)] px-2 py-2 text-left'
                    : 'block w-full cursor-pointer px-2 py-2 text-left hover:bg-[var(--color-border)]'
                }
                onClick={() => deliver(match)}
              >
                {label}
              </button>
            );
          })}
        </div>
      ) : null}
    </div>
  );
}
