import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import {
  DEFAULT_HG38_ASSEMBLIES,
  assemblyUsesHg38,
  classifyQuery,
  PINNED_HG38_SEQUENCES,
  formatRegion,
  geneIndexRequestUrl,
  loadHg38GeneIndex,
  missingPinnedSequence,
  parseGeneIndex,
  pinnedHg38Has,
  resetGeneIndexCacheForTests,
  searchGenes,
} from './geneIndex';

const directory = dirname(fileURLToPath(import.meta.url));
const SHA = 'a'.repeat(64);

function fixtureIndex() {
  const rows = [
    ['BHLHE40', 'ENSG00000000005.1', 0, 1, 3, '+', 'protein_coding'],
    ['DELEC1', 'ENSG00000000006.1', 2, 4, 7, '+', 'lncRNA'],
    ['EGFR', 'ENSG00000000001.1', 1, 10, 40, '+', 'protein_coding'],
    ['IL9R', 'ENSG00000000003.1', 4, 5, 8, '+', 'protein_coding'],
    ['IL9R', 'ENSG00000000004.1', 5, 9, 12, '+', 'protein_coding'],
    ['TP53', 'ENSG00000000002.1', 3, 100, 110, '-', 'protein_coding'],
  ];
  const aliases = [
    ['DEC1', 0],
    ['DEC1', 1],
    ['ERBB1', 2],
    ['p53', 5],
  ];
  const chromosomes = ['chr3', 'chr7', 'chr9', 'chr17', 'chrX', 'chrY'];
  return {
    format: 'helixweave.gene-index.v1',
    assembly: 'hg38',
    assembly_match: [...DEFAULT_HG38_ASSEMBLIES],
    coordinate_system: '1-based-closed',
    release: 'GENCODE human 49',
    chromosomes,
    rows,
    aliases,
    provenance: {
      annotation: {
        release: 'GENCODE human 49',
        source_url:
          'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz',
        filename: 'annotation.gtf',
        bytes: 12,
        sha256: SHA,
        feature: 'gene',
        assembly_accession: 'GCA_000001405.15',
        license: 'open access',
      },
      coordinate_system: '1-based-closed',
      coordinate_note: '1-based closed',
      counts: {
        genes: rows.length,
        aliases: aliases.length,
        chromosomes: chromosomes.length,
        dropped_incomplete: 0,
        dropped_unknown_sequence: 0,
      },
      aliases: { sha256: SHA, license: 'CC0', source: 'HGNC complete set' },
      assembly_report: { sha256: SHA, filename: 'assembly_report.txt' },
    },
  };
}

describe('coordinate and assembly classification', () => {
  it('keeps the panel pattern for a legal region, including a reversed span', () => {
    expect(classifyQuery('  chr1:500-2700  ')).toEqual({
      kind: 'coordinate',
      region: 'chr1:500-2700',
    });
    expect(classifyQuery('chr1:10-1')).toEqual({ kind: 'coordinate', region: 'chr1:10-1' });
    expect(classifyQuery('chr1:500')).toEqual({ kind: 'coordinate-draft' });
    expect(classifyQuery('EGFR')).toEqual({ kind: 'gene', query: 'EGFR' });
    expect(classifyQuery('  ')).toEqual({ kind: 'empty' });
  });

  it('recognizes hg38 names and refuses another assembly', () => {
    expect(assemblyUsesHg38(' hg38 ')).toBe(true);
    expect(assemblyUsesHg38('GRCh38.p14')).toBe(true);
    expect(assemblyUsesHg38('GCA_000001405.15')).toBe(true);
    expect(assemblyUsesHg38('mm10')).toBe(false);
    expect(assemblyUsesHg38('hg19')).toBe(false);
    expect(assemblyUsesHg38('')).toBe(false);
  });

  it('names only the sequences present in the pinned hg38 genome', () => {
    expect(PINNED_HG38_SEQUENCES).toHaveLength(25);
    expect(pinnedHg38Has('chr7')).toBe(true);
    expect(pinnedHg38Has('chrM')).toBe(true);
    expect(missingPinnedSequence('chrY')).toBeNull();
    expect(missingPinnedSequence('chr1_KI270706v1_random')).toMatch(/viewer was not moved/);
    expect(missingPinnedSequence('CHR7')).toMatch(/viewer was not moved/);
  });

  it('refuses a gene index URL on another origin', () => {
    const origin = 'http://localhost:3000';
    expect(geneIndexRequestUrl('/gene-index/hg38.json', origin)).toBe(
      'http://localhost:3000/gene-index/hg38.json',
    );
    expect(() => geneIndexRequestUrl('https://example.com/genes.json?q=EGFR', origin)).toThrow(
      /same-origin/,
    );
    expect(() => geneIndexRequestUrl('//genes.example/hg38.json', origin)).toThrow(/same-origin/);
  });
});

describe('bundled index request', () => {
  beforeEach(() => {
    resetGeneIndexCacheForTests();
  });

  it('requests only the same-origin index and does not put the query in the URL', async () => {
    const fetcher = vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
      void input;
      void init;
      return new Response(JSON.stringify(fixtureIndex()), { status: 200 });
    });
    const loaded = await loadHg38GeneIndex(fetcher);
    expect(fetcher).toHaveBeenCalledTimes(1);
    const call = fetcher.mock.calls[0];
    if (!call) throw new Error('missing gene index request');
    const [url, init] = call;
    const requested = String(url);
    expect(requested.toLowerCase()).not.toContain('egfr');
    expect(new URL(requested).origin).toBe(window.location.origin);
    expect(init).toEqual({ mode: 'same-origin', credentials: 'same-origin' });
    expect(searchGenes(loaded, 'ERBB1').matches.map((match) => match.locus.symbol)).toEqual([
      'EGFR',
    ]);
  });

  it('leaves a failed index request for the caller to keep coordinates usable', async () => {
    const fetcher = vi.fn(async () => new Response('missing', { status: 404 }));
    await expect(loadHg38GeneIndex(fetcher)).rejects.toThrow(/failed \(404\)/);
  });
});

describe('fixture search', () => {
  const index = parseGeneIndex(fixtureIndex());

  it('resolves one symbol, an alias, a versionless id, and two same-name loci', () => {
    expect(formatRegion(searchGenes(index, 'egfr').matches[0]?.locus ?? fail())).toBe('chr7:10-40');
    expect(searchGenes(index, 'p53').matches.map((match) => [match.via, match.locus.symbol])).toEqual(
      [['alias', 'TP53']],
    );
    expect(
      searchGenes(index, 'ENSG00000000001').matches.map((match) => match.locus.symbol),
    ).toEqual(['EGFR']);
    expect(searchGenes(index, 'IL9R').matches.map((match) => formatRegion(match.locus))).toEqual([
      'chrX:5-8',
      'chrY:9-12',
    ]);
    expect(searchGenes(index, 'DEC1').matches.map((match) => match.locus.symbol)).toEqual([
      'BHLHE40',
      'DELEC1',
    ]);
    expect(searchGenes(index, 'chr7:10-40').matches).toEqual([]);
    expect(searchGenes(index, 'ZZZNOGENE').relation).toBe('none');
  });
});

describe('shipped GRCh38 index', () => {
  const index = parseGeneIndex(
    JSON.parse(readFileSync(join(directory, 'hg38-gene-index.json'), 'utf8')),
  );
  const provenance = JSON.parse(
    readFileSync(join(directory, 'hg38-gene-index.provenance.json'), 'utf8'),
  );

  it('is the hashed GENCODE 49 extract, not a handwritten locus table', () => {
    expect(index.provenance).toEqual(provenance);
    expect(index.release).toBe('GENCODE human 49');
    expect(index.assemblyMatch).toEqual([...DEFAULT_HG38_ASSEMBLIES]);
    expect(provenance.annotation.sha256).toBe(
      '5cabe88e93064851dc6cb826efe0e8b7babf7985e73df280cfb05910d05cd023',
    );
    expect(provenance.annotation.bytes).toBe(93496634);
    expect(provenance.annotation.source_url).toContain(
      'gencode.v49.primary_assembly.annotation.gtf',
    );
    expect(provenance.aliases.source_url).toContain('hgnc_complete_set_2025-10-07.tsv');
    expect(provenance.viewer_scaffold.sequences).toEqual([...PINNED_HG38_SEQUENCES]);
    expect(index.chromosomes.filter(pinnedHg38Has).sort()).toEqual([...PINNED_HG38_SEQUENCES].sort());
    expect(index.provenance.annotation.assembly_accession).toBe('GCA_000001405.29');
    expect(index.provenance.assembly_report.sha256).toBe(
      'aa733ce92719f6339b3a4540137d4322894434363314c4b7704ba6691e9eab66',
    );
    expect(index.provenance.aliases.sha256).toBe(
      '30c2a235934e0f2ec579cb5a7158b8df0091a6674793341d08b8ea9e42e797e5',
    );
    expect(index.provenance.counts.genes).toBe(index.rows.length);
    expect(index.provenance.counts.genes).toBeGreaterThan(70000);
    expect(index.provenance.counts.dropped_unknown_sequence).toBe(0);
    expect(provenance.counts.outside_viewer_scaffold).toBe(208);
    for (const name of ['LocusSearch.tsx', 'geneIndex.ts', 'build_grch38_gene_index.py']) {
      const source = readFileSync(join(directory, name), 'utf8');
      expect(source).not.toContain('55018820');
      expect(source).not.toContain('55019017');
    }
  });

  it('resolves EGFR, TP53, same-name IL9R, and alias DEC1 from the annotation', () => {
    const egfr = searchGenes(index, 'EGFR');
    expect(egfr.relation).toBe('exact');
    expect(egfr.matches.map((match) => [match.via, match.locus.geneId, formatRegion(match.locus)])).toEqual(
      [['symbol', 'ENSG00000146648.22', 'chr7:55018820-55211628']],
    );
    expect(egfr.matches[0]?.locus.strand).toBe('+');
    expect(egfr.matches[0]?.locus.geneType).toBe('protein_coding');

    expect(
      searchGenes(index, 'tp53').matches.map((match) => [
        match.locus.geneId,
        formatRegion(match.locus),
        match.locus.strand,
      ]),
    ).toEqual([['ENSG00000141510.20', 'chr17:7661779-7687546', '-']]);

    expect(
      searchGenes(index, 'IL9R').matches.map((match) => [
        formatRegion(match.locus),
        match.locus.geneId,
      ]),
    ).toEqual([
      ['chrX:155997696-156013017', 'ENSG00000124334.20'],
      ['chrY:57184216-57199537', 'ENSG00000292373.3'],
    ]);

    expect(
      searchGenes(index, 'dec1').matches.map((match) => [
        match.via,
        match.matchedAlias,
        match.locus.symbol,
        match.locus.geneId,
        formatRegion(match.locus),
        match.locus.geneType,
      ]),
    ).toEqual([
      ['alias', 'DEC1', 'BHLHE40', 'ENSG00000134107.6', 'chr3:4978715-4985323', 'protein_coding'],
      ['alias', 'DEC1', 'DELEC1', 'ENSG00000173077.17', 'chr9:114850968-115402644', 'lncRNA'],
    ]);

    expect(searchGenes(index, 'ERBB1').matches.map((match) => match.locus.symbol)).toEqual(['EGFR']);
    expect(searchGenes(index, 'P53').matches.map((match) => [match.via, match.locus.symbol])).toEqual([
      ['alias', 'TP53'],
    ]);
    expect(searchGenes(index, 'ensg00000146648').matches.map((match) => match.via)).toEqual([
      'gene-id',
    ]);
    expect(searchGenes(index, 'chr7:55018820-55211628').matches).toEqual([]);
  });

  it('preserves both pseudoautosomal locations for HGNC synonyms and names unavailable sequences', () => {
    for (const term of ['SHOX', 'PHOG', 'SHOXY']) {
      expect(searchGenes(index, term).matches.map(({locus}) => formatRegion(locus))).toEqual([
        'chrX:624344-659411', 'chrY:624344-659411',
      ]);
    }
    const [offScaffold] = searchGenes(index, 'ENSG00000278457.1').matches;
    expect(offScaffold?.locus.chrom).toBe('chrUn_KI270442v1');
    expect(pinnedHg38Has(offScaffold!.locus.chrom)).toBe(false);
    expect(searchGenes(index, 'IMAGE:4942737').matches.map(({locus}) => locus.symbol)).toEqual(['FAM20C']);
  });
});

function fail(): never {
  throw new Error('missing locus');
}
