import { act, fireEvent, render, screen } from '@testing-library/react';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import LocusSearch from './LocusSearch';
import { DEFAULT_HG38_ASSEMBLIES, parseGeneIndex } from './geneIndex';

const SHA = 'b'.repeat(64);

vi.mock('./geneIndex', async () => {
  const actual = await vi.importActual<typeof import('./geneIndex')>('./geneIndex');
  return {
    ...actual,
    loadHg38GeneIndex: vi.fn(),
  };
});

const { loadHg38GeneIndex } = await import('./geneIndex');
const loadIndex = vi.mocked(loadHg38GeneIndex);

function parsedFixture() {
  const rows = [
    ['BHLHE40', 'ENSG00000000005.1', 0, 1, 3, '+', 'protein_coding'],
    ['DELEC1', 'ENSG00000000006.1', 2, 4, 7, '+', 'lncRNA'],
    ['EGFR', 'ENSG00000000001.1', 1, 10, 40, '+', 'protein_coding'],
    ['IL9R', 'ENSG00000000003.1', 4, 5, 8, '+', 'protein_coding'],
    ['IL9R', 'ENSG00000000004.1', 5, 9, 12, '+', 'protein_coding'],
    ['TP53', 'ENSG00000000002.1', 3, 100, 110, '-', 'protein_coding'],
    ['OFFSCAFFOLD', 'ENSG00000000007.1', 6, 5, 9, '+', 'lncRNA'],
  ];
  const aliases = [
    ['DEC1', 0],
    ['DEC1', 1],
    ['ERBB1', 2],
    ['p53', 5],
    ['MIXED', 2], ['MIXED', 6],
  ];
  const viewerChromosomes = ['chr3', 'chr7', 'chr9', 'chr17', 'chrX', 'chrY'];
  const chromosomes = [...viewerChromosomes, 'chrUn_KI270442v1'];
  return parseGeneIndex({
    format: 'helixweave.gene-index.v1',
    assembly: 'hg38',
    assembly_match: [...DEFAULT_HG38_ASSEMBLIES],
    coordinate_system: '1-based-closed',
    release: 'GENCODE human 49',
    chromosomes,
    viewer_chromosomes: viewerChromosomes,
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
  });
}

function input() {
  return screen.getByRole('textbox', { name: 'Genomic region' });
}

async function settleIndex() {
  await act(async () => {
    await Promise.resolve();
  });
}

describe('LocusSearch', () => {
  const fetchSpy = vi.fn();

  beforeEach(() => {
    loadIndex.mockReset();
    loadIndex.mockResolvedValue(parsedFixture());
    fetchSpy.mockReset();
    vi.stubGlobal('fetch', fetchSpy);
  });

  it('shows the current region and does not navigate until the user edits it', async () => {
    const onNavigate = vi.fn();
    const view = render(
      <LocusSearch assembly="hg38" region="chr1:500-2700" onNavigate={onNavigate} />,
    );
    await settleIndex();
    expect(screen.getByText('Genomic region')).toBeInTheDocument();
    expect(input()).toHaveValue('chr1:500-2700');
    expect(onNavigate).not.toHaveBeenCalled();
    view.rerender(<LocusSearch assembly="hg38" region="chr2:3-9" onNavigate={onNavigate} />);
    expect(input()).toHaveValue('chr2:3-9');
    expect(onNavigate).not.toHaveBeenCalled();
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it('navigates as soon as a legal coordinate is typed, and trims it', async () => {
    const onNavigate = vi.fn();
    render(<LocusSearch assembly="hg38" region="" onNavigate={onNavigate} />);
    await settleIndex();
    fireEvent.change(input(), { target: { value: '  chr1:500-2700  ' } });
    expect(input()).toHaveValue('chr1:500-2700');
    expect(onNavigate).toHaveBeenCalledTimes(1);
    expect(onNavigate).toHaveBeenCalledWith('chr1:500-2700');
    expect(screen.queryByRole('option')).not.toBeInTheDocument();

    fireEvent.change(input(), { target: { value: 'chr1:500' } });
    expect(onNavigate).toHaveBeenCalledTimes(1);
    expect(screen.getByRole('status')).toHaveTextContent('Coordinates use chrom:start-end.');

    fireEvent.change(input(), { target: { value: 'chr1:10-1' } });
    expect(onNavigate).toHaveBeenLastCalledWith('chr1:10-1');
  });

  it.each([
    ['EGFR', undefined, ['chr7:10-40']],
    ['EGFR', 'TP53', []],
    ['EGFR', 'chr7:3-8', ['chr7:3-8']],
    ['IL9R', undefined, []],
  ] as const)('honors a pending Enter for %s when the input becomes %s', async (query, edit, expected) => {
    let release!: (value: ReturnType<typeof parsedFixture>) => void;
    loadIndex.mockReturnValueOnce(new Promise(resolve => { release = resolve; }));
    const onNavigate = vi.fn();
    render(<LocusSearch assembly="hg38" region="" onNavigate={onNavigate} />);
    fireEvent.change(input(), { target: { value: query } });
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).not.toHaveBeenCalled();
    if (edit) fireEvent.change(input(), { target: { value: edit } });
    await act(async () => release(parsedFixture()));
    expect(onNavigate.mock.calls.map(([region]) => region)).toEqual(expected);
  });

  it('explains out-of-scaffold loci and allows the user to choose an available alias hit', async () => {
    const onNavigate = vi.fn();
    render(<LocusSearch assembly="hg38" region="chr7:3-8" onNavigate={onNavigate} />);
    fireEvent.change(input(), { target: { value: 'OFFSCAFFOLD' } });
    await screen.findByText(/WashU hg38 cannot display chrUn_KI270442v1/);
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(screen.queryByRole('option')).not.toBeInTheDocument();
    fireEvent.change(input(), { target: { value: 'MIXED' } });
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).not.toHaveBeenCalled();
    expect(screen.getByRole('status')).toHaveTextContent('cannot display chrUn_KI270442v1');
    fireEvent.click(screen.getByRole('option', { name: /EGFR/ }));
    expect(onNavigate).toHaveBeenCalledWith('chr7:10-40');
  });

  it('chooses one gene from the index and keeps an ambiguous name until a location is chosen', async () => {
    const onNavigate = vi.fn();
    render(<LocusSearch assembly="hg38" region="" onNavigate={onNavigate} />);
    fireEvent.change(input(), { target: { value: 'EGFR' } });
    const egfr = await screen.findByRole('option', { name: /EGFR/ });
    expect(egfr).toHaveAttribute('data-region', 'chr7:10-40');
    expect(onNavigate).not.toHaveBeenCalled();
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).toHaveBeenCalledWith('chr7:10-40');

    onNavigate.mockClear();
    fireEvent.change(input(), { target: { value: 'IL9R' } });
    const loci = await screen.findAllByRole('option');
    expect(loci.map((option) => option.getAttribute('data-region'))).toEqual([
      'chrX:5-8',
      'chrY:9-12',
    ]);
    expect(screen.getByRole('status')).toHaveTextContent('IL9R matches 2 locations. Choose one.');
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).not.toHaveBeenCalled();
    expect(input()).not.toHaveAttribute('aria-activedescendant');
    fireEvent.keyDown(input(), { key: 'ArrowDown' });
    expect(input()).toHaveAttribute('aria-activedescendant', loci[0]!.id);
    fireEvent.keyDown(input(), { key: 'ArrowDown' });
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).toHaveBeenCalledWith('chrY:9-12');

    onNavigate.mockClear();
    fireEvent.change(input(), { target: { value: 'dec1' } });
    const aliases = await screen.findAllByRole('option');
    expect(aliases.map((option) => option.textContent)).toEqual([
      'DEC1, alias of BHLHE40 · chr3:1-3 · + · protein_coding · ENSG00000000005.1',
      'DEC1, alias of DELEC1 · chr9:4-7 · + · lncRNA · ENSG00000000006.1',
    ]);
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).not.toHaveBeenCalled();
    fireEvent.click(aliases[0]!);
    expect(onNavigate).toHaveBeenCalledWith('chr3:1-3');
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it('requires an explicit choice even when a prefix has just one suggestion', async () => {
    const onNavigate = vi.fn();
    render(<LocusSearch assembly="hg38" region="" onNavigate={onNavigate} />);
    fireEvent.change(input(), { target: { value: 'EGF' } });
    const option = await screen.findByRole('option', { name: /EGFR/ });
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).not.toHaveBeenCalled();
    fireEvent.keyDown(input(), { key: 'ArrowUp' });
    expect(option).toHaveAttribute('aria-selected', 'true');
    fireEvent.keyDown(input(), { key: 'Escape' });
    expect(screen.queryByRole('option')).not.toBeInTheDocument();
    fireEvent.keyDown(input(), { key: 'Enter' });
    fireEvent.keyDown(input(), { key: 'Enter' });
    expect(onNavigate).not.toHaveBeenCalled();
    fireEvent.click(screen.getByRole('option', { name: /EGFR/ }));
    expect(onNavigate).toHaveBeenCalledWith('chr7:10-40');
  });

  it('says when nothing matches and still accepts coordinates without an annotation', async () => {
    const onNavigate = vi.fn();
    const { rerender } = render(<LocusSearch assembly="hg38" region="" onNavigate={onNavigate} />);
    fireEvent.change(input(), { target: { value: 'ZZZNOGENE' } });
    expect(await screen.findByRole('status')).toHaveTextContent(
      'No genes match "ZZZNOGENE" in GENCODE human 49. Coordinates still work.',
    );
    expect(screen.queryByRole('option')).not.toBeInTheDocument();

    loadIndex.mockRejectedValueOnce(new Error('offline'));
    rerender(<LocusSearch assembly="GRCh38" region="" onNavigate={onNavigate} />);
    expect(await screen.findByRole('status')).toHaveTextContent(
      'Gene annotation did not load. Coordinates still work.',
    );
    fireEvent.change(input(), { target: { value: 'chr1:8-9' } });
    expect(onNavigate).toHaveBeenCalledWith('chr1:8-9');
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it('keeps coordinates usable when this assembly has no bundled annotation', () => {
    const onNavigate = vi.fn();
    render(<LocusSearch assembly="mm10" region="" onNavigate={onNavigate} />);
    expect(screen.getByRole('status')).toHaveTextContent(
      'No bundled gene annotation for mm10. Coordinates still work.',
    );
    expect(loadIndex).not.toHaveBeenCalled();
    fireEvent.change(input(), { target: { value: 'chr1:1-5' } });
    expect(onNavigate).toHaveBeenCalledWith('chr1:1-5');
    fireEvent.change(input(), { target: { value: 'EGFR' } });
    expect(onNavigate).toHaveBeenCalledTimes(1);
    expect(screen.queryByRole('option')).not.toBeInTheDocument();
    expect(fetchSpy).not.toHaveBeenCalled();
  });
});
