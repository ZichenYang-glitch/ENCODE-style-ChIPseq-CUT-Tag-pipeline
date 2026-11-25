import assert from 'node:assert/strict';
import test from 'node:test';
import zlib from 'node:zlib';
import {
  bufferHasFloat,
  inflateBgzf,
  isExternalUrl,
  isLocalUpload,
  isWashuProxy,
  pairLines,
  parseRegion,
  regionCovers,
  sumSliceBytes,
  textHasGene,
} from './washu-local-overlay-observe.mjs';

test('gene words are caught in any case and coordinates are not gene words', () => {
  assert.equal(textHasGene('lookup EGFR on chr7'), true);
  assert.equal(textHasGene('no_such_gene_helixweave'), true);
  assert.equal(textHasGene('chr7:55018900-55019400'), false);
});

test('only off-origin hosts and the WashU proxy count as egress', () => {
  assert.equal(isExternalUrl('http://127.0.0.1:18765/api/v1/runs'), false);
  assert.equal(isExternalUrl('blob:http://127.0.0.1/1234'), false);
  assert.equal(isExternalUrl('https://epigenome.wustl.edu/cors/file'), true);
  assert.equal(isWashuProxy('https://epigenome.wustl.edu/cors/file'), true);
  assert.equal(isWashuProxy('http://127.0.0.1/api/v1/runs'), false);
});

test('local markers and oversized posts are uploads', () => {
  const markers = ['egfr-peak'];
  assert.equal(isLocalUpload('POST', '{"ok":true}', markers), false);
  assert.equal(isLocalUpload('POST', 'track\tegfr-peak\t80', markers), true);
  assert.equal(isLocalUpload('PUT', 'x'.repeat(8193), markers), true);
  assert.equal(isLocalUpload('GET', 'x'.repeat(9000), markers), false);
});

test('slice totals stay with the named file', () => {
  const records = [
    { name: 'local-pet.longrange.gz', bytes: 78 },
    { name: 'local-pet.longrange.gz.tbi', bytes: 5000 },
    { name: 'local-pet.longrange.gz', bytes: 12 },
  ];
  assert.equal(sumSliceBytes(records, 'local-pet.longrange.gz'), 90);
  assert.equal(sumSliceBytes(records, 'local-pet.longrange.gz.tbi'), 5000);
});

test('a compressed bigWig block still exposes its float32 signal', () => {
  const raw = Buffer.alloc(4);
  raw.writeFloatLE(4.5, 0);
  const zipped = zlib.deflateSync(raw);
  assert.equal(bufferHasFloat(zipped, 4.5), true);
  assert.equal(bufferHasFloat(zipped, 2.25), false);
});

test('regions use 0-based half-open coordinates', () => {
  const region = parseRegion('chr7:55018900-55019400');
  assert.equal(regionCovers(region, 'chr7', 55019000, 55019300), true);
  assert.equal(regionCovers(region, 'chr7', 55115200, 55115500), false);
  assert.equal(parseRegion('EGFR'), null);
});

test('a planted pair round-trips to two longrange anchors', () => {
  const [left, right] = pairLines('chr7', { a: [55019000, 55019100], b: [55019200, 55019300] });
  assert.equal(left, 'chr7\t55019000\t55019100\tchr7:55019200-55019300,1');
  assert.equal(right, 'chr7\t55019200\t55019300\tchr7:55019000-55019100,1');
});

test('non-bgzf bytes do not invent records', () => {
  assert.equal(inflateBgzf(Buffer.from('not a track')).length, 0);
});
