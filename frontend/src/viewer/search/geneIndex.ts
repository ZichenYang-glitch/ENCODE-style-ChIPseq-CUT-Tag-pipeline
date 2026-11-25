import bundledIndexUrl from './hg38-gene-index.json?url';

/**
 * Gene spans in the bundled index are 1-based and closed, the GTF's own
 * coordinates. A typed chrom:start-end is returned unchanged: the same
 * pattern the genome panel already treats as a legal region. The provenance
 * URL is a citation, not a request.
 */

const LEGAL_COORDINATE = /^([A-Za-z0-9_.:-]+):([0-9]+)-([0-9]+)$/;
const GENE_ID = /^ENSG\d+\.\d+$/;
const SHA256 = /^[0-9a-f]{64}$/;

export const DEFAULT_HG38_ASSEMBLIES = [
  'hg38',
  'GRCh38',
  'GRCh38.p14',
  'GCA_000001405.15',
] as const;

const PREFIX_MIN = 2;
const PREFIX_LIMIT = 12;
const VIA_RANK = { symbol: 0, alias: 1, 'gene-id': 2, prefix: 3 } as const;

type GeneRow = [
  symbol: string,
  geneId: string,
  chromIndex: number,
  start: number,
  end: number,
  strand: string,
  geneType: string,
];

export interface Locus {
  symbol: string;
  geneId: string;
  chrom: string;
  start: number;
  end: number;
  strand: string;
  geneType: string;
}

export interface GeneMatch {
  locus: Locus;
  via: 'symbol' | 'alias' | 'gene-id' | 'prefix';
  matchedAlias?: string;
}

export interface GeneSearch {
  matches: GeneMatch[];
  relation: 'exact' | 'prefix' | 'none';
  total: number;
}

export interface GeneProvenance {
  annotation: {
    release: string;
    source_url: string;
    filename: string;
    bytes: number;
    sha256: string;
    feature: string;
    assembly_accession: string;
    license: string;
  };
  coordinate_system: string;
  coordinate_note: string;
  counts: {
    genes: number;
    aliases: number;
    chromosomes: number;
    dropped_incomplete: number;
    dropped_unknown_sequence: number;
  };
  aliases: {
    sha256: string;
    license: string;
    source: string;
  };
  assembly_report: {
    sha256: string;
    filename: string;
  };
}

export interface ParsedGeneIndex {
  release: string;
  assemblyMatch: readonly string[];
  chromosomes: readonly string[];
  provenance: GeneProvenance;
  rows: readonly GeneRow[];
  bySymbol: ReadonlyMap<string, readonly number[]>;
  byAlias: ReadonlyMap<string, readonly number[]>;
  byGeneId: ReadonlyMap<string, readonly number[]>;
  aliasSpelling: ReadonlyMap<string, string>;
  symbolsSorted: readonly { key: string; row: number }[];
}

export type QueryPlan =
  | { kind: 'empty' }
  | { kind: 'coordinate'; region: string }
  | { kind: 'coordinate-draft' }
  | { kind: 'gene'; query: string };

export function legalCoordinate(value: string): string | null {
  const trimmed = value.trim();
  if (!LEGAL_COORDINATE.test(trimmed)) return null;
  return trimmed;
}

export function classifyQuery(value: string): QueryPlan {
  const trimmed = value.trim();
  if (!trimmed) return { kind: 'empty' };
  const region = legalCoordinate(trimmed);
  if (region) return { kind: 'coordinate', region };
  if (trimmed.includes(':')) return { kind: 'coordinate-draft' };
  return { kind: 'gene', query: trimmed };
}

export function assemblyUsesHg38(
  assembly: string,
  accepted: readonly string[] = DEFAULT_HG38_ASSEMBLIES,
): boolean {
  const key = assembly.trim().toLowerCase();
  if (!key) return false;
  return accepted.some((item) => item.toLowerCase() === key);
}

export function formatRegion(locus: Locus): string {
  return `${locus.chrom}:${locus.start}-${locus.end}`;
}

export function matchLabel(match: GeneMatch): string {
  const { locus } = match;
  const name =
    match.via === 'alias' && match.matchedAlias
      ? `${match.matchedAlias}, alias of ${locus.symbol}`
      : locus.symbol;
  return `${name} · ${formatRegion(locus)} · ${locus.strand} · ${locus.geneType} · ${locus.geneId}`;
}

export function geneIndexRequestUrl(indexUrl: string, origin: string): string {
  let url: URL;
  try {
    url = new URL(indexUrl, origin);
  } catch {
    throw new Error('gene index URL is not same-origin');
  }
  if (url.origin !== new URL(origin).origin || url.username || url.password) {
    throw new Error('gene index URL is not same-origin');
  }
  return url.href;
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function requiredString(record: Record<string, unknown>, key: string): string {
  const value = record[key];
  if (typeof value !== 'string' || value.trim() === '') {
    throw new Error(`gene index ${key} is missing`);
  }
  return value;
}

function requiredSha256(record: Record<string, unknown>, key: string): string {
  const value = requiredString(record, key);
  if (!SHA256.test(value)) throw new Error(`gene index ${key} is not a sha256`);
  return value;
}

function requiredCount(record: Record<string, unknown>, key: string): number {
  const value = record[key];
  if (typeof value !== 'number' || !Number.isSafeInteger(value) || value < 0) {
    throw new Error(`gene index ${key} is not a count`);
  }
  return value;
}

function parseProvenance(value: unknown): GeneProvenance {
  if (!isRecord(value)) throw new Error('gene index provenance is missing');
  const annotation = value.annotation;
  const counts = value.counts;
  const aliases = value.aliases;
  const report = value.assembly_report;
  if (!isRecord(annotation) || !isRecord(counts) || !isRecord(aliases) || !isRecord(report)) {
    throw new Error('gene index provenance is incomplete');
  }
  const bytes = annotation.bytes;
  if (typeof bytes !== 'number' || !Number.isSafeInteger(bytes) || bytes < 1) {
    throw new Error('gene index annotation length is missing');
  }
  requiredString(annotation, 'release');
  requiredString(annotation, 'source_url');
  requiredString(annotation, 'filename');
  requiredSha256(annotation, 'sha256');
  requiredString(annotation, 'feature');
  requiredString(annotation, 'assembly_accession');
  requiredString(annotation, 'license');
  requiredString(value, 'coordinate_system');
  requiredString(value, 'coordinate_note');
  requiredCount(counts, 'genes');
  requiredCount(counts, 'aliases');
  requiredCount(counts, 'chromosomes');
  requiredCount(counts, 'dropped_incomplete');
  requiredCount(counts, 'dropped_unknown_sequence');
  requiredSha256(aliases, 'sha256');
  requiredString(aliases, 'license');
  requiredString(aliases, 'source');
  requiredSha256(report, 'sha256');
  requiredString(report, 'filename');
  return value as unknown as GeneProvenance;
}

function parseRow(value: unknown, chromosomes: readonly string[]): GeneRow {
  if (!Array.isArray(value) || value.length !== 7) {
    throw new Error('gene index row is malformed');
  }
  const [symbol, geneId, chromIndex, start, end, strand, geneType] = value;
  if (
    typeof symbol !== 'string' ||
    symbol.trim() === '' ||
    typeof geneId !== 'string' ||
    !GENE_ID.test(geneId) ||
    typeof chromIndex !== 'number' ||
    !Number.isSafeInteger(chromIndex) ||
    chromIndex < 0 ||
    chromIndex >= chromosomes.length ||
    typeof start !== 'number' ||
    typeof end !== 'number' ||
    !Number.isSafeInteger(start) ||
    !Number.isSafeInteger(end) ||
    start < 1 ||
    end < start ||
    (strand !== '+' && strand !== '-' && strand !== '.') ||
    typeof geneType !== 'string' ||
    geneType.trim() === ''
  ) {
    throw new Error('gene index row is malformed');
  }
  return [symbol, geneId, chromIndex, start, end, strand, geneType];
}

function pushRow(map: Map<string, number[]>, key: string, row: number): void {
  const list = map.get(key);
  if (!list) {
    map.set(key, [row]);
    return;
  }
  if (!list.includes(row)) list.push(row);
}

export function parseGeneIndex(value: unknown): ParsedGeneIndex {
  if (!isRecord(value)) throw new Error('gene index is not an object');
  if (value.format !== 'helixweave.gene-index.v1') {
    throw new Error('unsupported gene index');
  }
  if (value.assembly !== 'hg38' || value.coordinate_system !== '1-based-closed') {
    throw new Error('gene index is not the hg38 1-based annotation');
  }
  const release = requiredString(value, 'release');
  const provenance = parseProvenance(value.provenance);
  if (provenance.annotation.release !== release) {
    throw new Error('gene index release does not match its provenance');
  }
  if (provenance.coordinate_system !== '1-based-closed') {
    throw new Error('gene index provenance is not 1-based');
  }
  if (!Array.isArray(value.assembly_match) || value.assembly_match.length === 0) {
    throw new Error('gene index assembly list is missing');
  }
  const assemblyMatch = value.assembly_match.map((item) => {
    if (typeof item !== 'string' || item.trim() === '') {
      throw new Error('gene index assembly list is malformed');
    }
    return item;
  });
  if (!assemblyUsesHg38('hg38', assemblyMatch)) {
    throw new Error('gene index does not apply to hg38');
  }
  if (!Array.isArray(value.chromosomes) || value.chromosomes.length === 0) {
    throw new Error('gene index chromosomes are missing');
  }
  const chromosomes = value.chromosomes.map((item) => {
    if (typeof item !== 'string' || item.trim() === '') {
      throw new Error('gene index chromosome name is missing');
    }
    return item;
  });
  if (!Array.isArray(value.rows)) throw new Error('gene index rows are missing');
  const rows = value.rows.map((row) => parseRow(row, chromosomes));
  if (!Array.isArray(value.aliases)) throw new Error('gene index aliases are missing');
  const bySymbol = new Map<string, number[]>();
  const byAlias = new Map<string, number[]>();
  const byGeneId = new Map<string, number[]>();
  const aliasSpelling = new Map<string, string>();
  rows.forEach((row, index) => {
    pushRow(bySymbol, row[0].toUpperCase(), index);
    pushRow(byGeneId, row[1].toUpperCase(), index);
    const versionless = row[1].toUpperCase().split('.')[0] ?? row[1].toUpperCase();
    if (versionless !== row[1].toUpperCase()) pushRow(byGeneId, versionless, index);
  });
  for (const entry of value.aliases) {
    if (!Array.isArray(entry) || entry.length !== 2) {
      throw new Error('gene index alias is malformed');
    }
    const [alias, row] = entry;
    if (
      typeof alias !== 'string' ||
      alias.trim() === '' ||
      typeof row !== 'number' ||
      !Number.isSafeInteger(row) ||
      row < 0 ||
      row >= rows.length
    ) {
      throw new Error('gene index alias is malformed');
    }
    const key = alias.toUpperCase();
    if (!aliasSpelling.has(key)) aliasSpelling.set(key, alias);
    pushRow(byAlias, key, row);
  }
  if (
    provenance.counts.genes !== rows.length ||
    provenance.counts.aliases !== value.aliases.length ||
    provenance.counts.chromosomes !== chromosomes.length
  ) {
    throw new Error('gene index counts do not match its rows');
  }
  const symbolsSorted = rows
    .map((row, index) => ({ key: row[0].toUpperCase(), row: index }))
    .sort((left, right) => {
      if (left.key < right.key) return -1;
      if (left.key > right.key) return 1;
      return left.row - right.row;
    });
  return {
    release,
    assemblyMatch,
    chromosomes,
    provenance,
    rows,
    bySymbol,
    byAlias,
    byGeneId,
    aliasSpelling,
    symbolsSorted,
  };
}

function locusAt(index: ParsedGeneIndex, row: number): Locus {
  const item = index.rows[row];
  if (!item) throw new Error('gene index row is missing');
  const chrom = index.chromosomes[item[2]];
  if (!chrom) throw new Error('gene index chromosome is missing');
  return {
    symbol: item[0],
    geneId: item[1],
    chrom,
    start: item[3],
    end: item[4],
    strand: item[5],
    geneType: item[6],
  };
}

function chromRank(chrom: string): number {
  const body = chrom.replace(/^chr/i, '');
  if (/^\d+$/.test(body)) return Number(body);
  if (body === 'X') return 23;
  if (body === 'Y') return 24;
  if (body === 'M' || body === 'MT') return 25;
  return 1000;
}

function compareLoci(left: Locus, right: Locus): number {
  const chrom = chromRank(left.chrom) - chromRank(right.chrom);
  if (chrom !== 0) return chrom;
  if (left.chrom !== right.chrom) return left.chrom < right.chrom ? -1 : 1;
  if (left.start !== right.start) return left.start - right.start;
  if (left.end !== right.end) return left.end - right.end;
  if (left.symbol !== right.symbol) return left.symbol < right.symbol ? -1 : 1;
  if (left.geneId !== right.geneId) return left.geneId < right.geneId ? -1 : 1;
  return 0;
}

function addRows(
  index: ParsedGeneIndex,
  chosen: Map<number, GeneMatch>,
  rows: readonly number[] | undefined,
  via: GeneMatch['via'],
  matchedAlias?: string,
): void {
  if (!rows) return;
  for (const row of rows) {
    const current = chosen.get(row);
    if (current && VIA_RANK[current.via] <= VIA_RANK[via]) continue;
    chosen.set(row, {
      locus: locusAt(index, row),
      via,
      matchedAlias: via === 'alias' ? matchedAlias : undefined,
    });
  }
}

function prefixRows(index: ParsedGeneIndex, key: string): number[] {
  const sorted = index.symbolsSorted;
  let low = 0;
  let high = sorted.length;
  while (low < high) {
    const mid = (low + high) >> 1;
    const item = sorted[mid];
    if (item && item.key < key) low = mid + 1;
    else high = mid;
  }
  const found: number[] = [];
  const seen = new Set<number>();
  for (let i = low; i < sorted.length; i += 1) {
    const item = sorted[i];
    if (!item || !item.key.startsWith(key)) break;
    if (seen.has(item.row)) continue;
    seen.add(item.row);
    found.push(item.row);
  }
  found.sort((left, right) => compareLoci(locusAt(index, left), locusAt(index, right)));
  return found;
}

export function searchGenes(index: ParsedGeneIndex, query: string): GeneSearch {
  const plan = classifyQuery(query);
  if (plan.kind !== 'gene' || plan.query.length > 80) {
    return { matches: [], relation: 'none', total: 0 };
  }
  const key = plan.query.toUpperCase();
  const chosen = new Map<number, GeneMatch>();
  addRows(index, chosen, index.bySymbol.get(key), 'symbol');
  addRows(index, chosen, index.byAlias.get(key), 'alias', index.aliasSpelling.get(key));
  addRows(index, chosen, index.byGeneId.get(key), 'gene-id');
  if (chosen.size > 0) {
    const matches = [...chosen.values()].sort((left, right) => {
      const via = VIA_RANK[left.via] - VIA_RANK[right.via];
      if (via !== 0) return via;
      return compareLoci(left.locus, right.locus);
    });
    return { matches, relation: 'exact', total: matches.length };
  }
  if (key.length < PREFIX_MIN) return { matches: [], relation: 'none', total: 0 };
  const rows = prefixRows(index, key);
  return {
    matches: rows.slice(0, PREFIX_LIMIT).map((row) => ({
      locus: locusAt(index, row),
      via: 'prefix' as const,
    })),
    relation: rows.length > 0 ? 'prefix' : 'none',
    total: rows.length,
  };
}

let cached: Promise<ParsedGeneIndex> | null = null;

export function resetGeneIndexCacheForTests(): void {
  cached = null;
}

async function requestGeneIndex(fetcher: typeof fetch): Promise<ParsedGeneIndex> {
  const url = geneIndexRequestUrl(bundledIndexUrl, window.location.origin);
  const response = await fetcher(url, { credentials: 'same-origin' });
  if (!response.ok) {
    throw new Error(`gene index request failed (${response.status})`);
  }
  return parseGeneIndex(await response.json());
}

export function loadHg38GeneIndex(fetcher: typeof fetch = fetch): Promise<ParsedGeneIndex> {
  if (fetcher !== fetch) return requestGeneIndex(fetcher);
  if (!cached) {
    cached = requestGeneIndex(fetch).catch((error: unknown) => {
      cached = null;
      throw error;
    });
  }
  return cached;
}
