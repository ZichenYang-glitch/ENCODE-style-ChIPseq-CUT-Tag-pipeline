import { useEffect, useId, useRef, useState } from 'react';
import { useUiLanguage, type UiLanguage } from '../../lib/uiLanguage';
import type { ChangeEvent, KeyboardEvent, PointerEvent } from 'react';
import {
  annotationAssembly,
  classifyQuery,
  coordinateSequence,
  formatRegion,
  loadGeneIndex,
  matchLabel,
  missingPinnedSequence,
  searchGenes,
} from './geneIndex';
import type { AnnotationAssembly, GeneMatch, GeneSearch, ParsedGeneIndex, QueryPlan } from './geneIndex';

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
  language: UiLanguage,
): string | null {
  const english = language === 'en';
  if (phase === 'failed') {
    return english ? 'Gene annotation did not load. Coordinates still work.' : '基因注释加载失败，仍可输入坐标。';
  }
  if (phase === 'unsupported') {
    const name = assembly.trim() || (english ? 'this assembly' : '当前组装');
    return english
      ? `No bundled gene annotation for ${name}. Coordinates still work.`
      : `${name} 没有内置基因注释，仍可输入坐标。`;
  }
  if (plan.kind === 'empty' || plan.kind === 'coordinate') return null;
  if (plan.kind === 'coordinate-draft') return english ? 'Coordinates use chrom:start-end.' : '坐标格式为 chr:start-end。';
  if (phase === 'loading') return english ? `Loading ${assembly} gene annotation…` : `正在加载 ${assembly} 的基因注释…`;
  if (!result || result.relation === 'none') {
    return english
      ? `No genes match "${plan.query}" in ${release}. Coordinates still work.`
      : `在 ${release} 中未找到“${plan.query}”，仍可输入坐标。`;
  }
  if (result.relation === 'prefix') {
    return english
      ? `No exact gene, alias, or gene id for "${plan.query}". Showing ${result.matches.length} of ${result.total} names that start with it.`
      : `未找到与“${plan.query}”完全匹配的基因、别名或基因 ID。以此开头的名称共 ${result.total} 个，显示 ${result.matches.length} 个。`;
  }
  if (result.total === 1) {
    return english
      ? `${plan.query} matches 1 location. Press Enter or choose it.`
      : `“${plan.query}”匹配一个位置。按 Enter 或点选即可跳转。`;
  }
  return english
    ? `${plan.query} matches ${result.total} locations. Choose one. Arrow keys, then Enter.`
    : `“${plan.query}”匹配 ${result.total} 个位置。请点选，或用方向键选择后按 Enter。`;
}

export default function LocusSearch(props: LocusSearchProps) {
  const annotation = annotationAssembly(props.assembly);
  return <AssemblySearch key={annotation ?? props.assembly} {...props} annotation={annotation} />;
}

function AssemblySearch({ assembly, region, onNavigate, annotation }: LocusSearchProps & {
  annotation: AnnotationAssembly | null;
}) {
  const language = useUiLanguage();
  const listId = useId();
  const pendingEnter = useRef<string | null>(null);
  const [draft, setDraft] = useState(() => region.trim());
  const [blockedSequence, setBlockedSequence] = useState<string | null>(null);
  const seenRegion = useRef(region);
  if (seenRegion.current !== region) {
    seenRegion.current = region;
    const next = region.trim();
    setDraft(next);
    if (pendingEnter.current && pendingEnter.current !== next) pendingEnter.current = null;
    setBlockedSequence(null);
  }
  const [pinnedDraft, setPinnedDraft] = useState(draft);
  const [armedIndex, setArmedIndex] = useState<number | null>(null);
  const [dismissedQuery, setDismissedQuery] = useState<string | null>(null);
  if (pinnedDraft !== draft) {
    setPinnedDraft(draft);
    setArmedIndex(null);
    setDismissedQuery(null);
  }
  const [requested, setRequested] = useState(false);
  const [phase, setPhase] = useState<Phase>(annotation ? 'loading' : 'unsupported');
  const [index, setIndex] = useState<ParsedGeneIndex | null>(null);

  useEffect(() => {
    if (!annotation || !requested) return;
    let cancelled = false;
    loadGeneIndex(annotation)
      .then((loaded) => {
        if (cancelled) return;
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
  }, [annotation, requested]);

  function deliver(match: GeneMatch) {
    const absent = annotation ? missingPinnedSequence(match.locus.chrom, annotation) : null;
    if (absent) {
      setBlockedSequence(match.locus.chrom);
      return;
    }
    const next = formatRegion(match.locus);
    setBlockedSequence(null);
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
    const absent = annotation ? missingPinnedSequence(match.locus.chrom, annotation) : null;
    if (absent) {
      setBlockedSequence(match.locus.chrom);
      return;
    }
    const next = formatRegion(match.locus);
    setBlockedSequence(null);
    setDraft(next);
    onNavigate(next);
  }, [phase, index, draft, annotation, onNavigate]);

  const plan = classifyQuery(draft);
  const result =
    phase === 'ready' && index && plan.kind === 'gene' ? searchGenes(index, plan.query) : null;
  const matches = result?.matches ?? [];
  const open = matches.length > 0 && dismissedQuery !== draft;
  const active =
    armedIndex !== null && matches.length > 0 ? Math.min(armedIndex, matches.length - 1) : null;
  const status =
    (blockedSequence && annotation ? missingPinnedSequence(blockedSequence, annotation, language) : null) ??
    statusMessage(phase, annotation ?? assembly, plan, result, index?.release ?? assembly, language);

  function onChange(event: ChangeEvent<HTMLInputElement>) {
    const next = event.target.value.trim();
    setDraft(next);
    if (pendingEnter.current && pendingEnter.current !== next) pendingEnter.current = null;
    const nextPlan = classifyQuery(next);
    if (annotation && nextPlan.kind === 'gene') setRequested(true);
    if (nextPlan.kind !== 'coordinate') {
      setBlockedSequence(null);
      return;
    }
    pendingEnter.current = null;
    const sequence = coordinateSequence(nextPlan.region);
    const absent =
      annotation && sequence ? missingPinnedSequence(sequence, annotation) : null;
    if (absent) {
      setBlockedSequence(sequence);
      return;
    }
    setBlockedSequence(null);
    onNavigate(nextPlan.region);
  }

  function onKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    if (event.key === 'Escape') {
      setDismissedQuery(draft);
      setArmedIndex(null);
      return;
    }
    if (event.key !== 'Enter' && event.key !== 'ArrowDown' && event.key !== 'ArrowUp') return;
    if (classifyQuery(draft).kind !== 'gene') return;
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
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
    if (phase !== 'ready' || !index) {
      pendingEnter.current = draft;
      return;
    }
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

  function choose(match: GeneMatch, event: PointerEvent<HTMLButtonElement> | KeyboardEvent<HTMLButtonElement>) {
    event.preventDefault();
    deliver(match);
  }

  const example = annotation === 'mm10' ? 'Egfr' : 'EGFR';
  const hint = language === 'en'
    ? (annotation
      ? `Enter a gene name (e.g. ${example}) or coordinates (chr:start-end). Press Enter for one match. For multiple matches, click one or use the arrow keys and Enter.`
      : 'Enter coordinates (chr:start-end).')
    : (annotation
      ? `基因名（如 ${example}）或 chr:start-end。多个位置请点选，或用方向键选中后再按 Enter。唯一结果可直接按 Enter。`
      : '输入坐标 chr:start-end。');

  return (
    <div className="block text-sm" lang={language}>
      <label className="block">
        <span className="block max-w-md">
          {hint}
        </span>
        <input
          className="mt-1 w-full max-w-md border border-[var(--color-border)] px-2 py-1"
          value={draft}
          onChange={onChange}
          onKeyDown={onKeyDown}
          aria-label={language === 'en' ? 'Genomic region' : '基因组位置'}
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
          aria-label={language === 'en' ? 'Gene matches' : '基因匹配结果'}
          className="mt-1 max-h-60 w-full max-w-md overflow-auto border border-[var(--color-border)]"
        >
          {matches.map((match, option) => {
            const absent = annotation
              ? missingPinnedSequence(match.locus.chrom, annotation)
              : null;
            let label = matchLabel(match, language);
            if (absent) {
              label += language === 'en' ? ` · not in pinned ${annotation}` : ` · 超出 ${annotation} 的染色体范围`;
            }
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
                onPointerDown={(event) => choose(match, event)}
                onKeyDown={(event) => {
                  if (event.key === 'Enter' || event.key === ' ') choose(match, event);
                }}
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
