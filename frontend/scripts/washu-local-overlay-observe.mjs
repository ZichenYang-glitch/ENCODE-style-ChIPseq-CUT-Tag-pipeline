import zlib from 'node:zlib';

export const GENE_WORDS = ['EGFR', 'TP53', 'DDX11L16', 'NO_SUCH_GENE_HELIXWEAVE'];

export function anchorLine(chrom, start, end, mateStart, mateEnd) {
  return `${chrom}\t${start}\t${end}\t${chrom}:${mateStart}-${mateEnd},1`;
}

export function pairLines(chrom, pair) {
  return [
    anchorLine(chrom, pair.a[0], pair.a[1], pair.b[0], pair.b[1]),
    anchorLine(chrom, pair.b[0], pair.b[1], pair.a[0], pair.a[1]),
  ];
}

export function peakLine(peak) {
  return `${peak.chrom}\t${peak.start}\t${peak.end}\t${peak.name}\t${peak.score}\t+`;
}

export function parseRegion(value) {
  const match = /^([^:]+):(\d+)-(\d+)$/.exec(String(value ?? '').trim());
  if (!match) return null;
  const start = Number(match[2]);
  const end = Number(match[3]);
  if (!Number.isInteger(start) || !Number.isInteger(end) || end <= start) return null;
  return { chrom: match[1], start, end };
}

export function regionCovers(region, chrom, start, end) {
  return Boolean(
    region && region.chrom === chrom && region.start <= start && region.end >= end,
  );
}

export function textHasGene(value) {
  const haystack = String(value ?? '').toLowerCase();
  return GENE_WORDS.some((word) => haystack.includes(word.toLowerCase()));
}

export function isExternalUrl(url) {
  let parsed;
  try {
    parsed = new URL(url);
  } catch {
    return true;
  }
  if (parsed.protocol === 'blob:' || parsed.protocol === 'data:') return false;
  return parsed.hostname !== '127.0.0.1' && parsed.hostname !== 'localhost';
}

export function isWashuProxy(url) {
  return String(url).includes('epigenome.wustl.edu');
}

export function isLocalUpload(method, body, markers) {
  if (body == null || body === '') return false;
  const text = Buffer.isBuffer(body) ? body.toString('utf8') : String(body);
  if (markers.some((marker) => marker && text.includes(marker))) return true;
  const verb = String(method ?? '').toUpperCase();
  return (verb === 'POST' || verb === 'PUT' || verb === 'PATCH') && text.length > 8192;
}

export function sumSliceBytes(records, filename) {
  return records.reduce(
    (total, record) => (record.name === filename ? total + Number(record.bytes) : total),
    0,
  );
}

export function bufferHasFloat(buffer, value) {
  const needle = Buffer.alloc(4);
  needle.writeFloatLE(value, 0);
  if (buffer.includes(needle)) return true;
  for (let start = 0; start < buffer.length - 2; start += 1) {
    if (buffer[start] !== 0x78) continue;
    const flag = buffer[start + 1];
    if (flag !== 0x01 && flag !== 0x5e && flag !== 0x9c && flag !== 0xda) continue;
    for (let end = start + 8; end <= buffer.length; end += 1) {
      try {
        if (zlib.inflateSync(buffer.subarray(start, end)).includes(needle)) return true;
      } catch {
        // This slice is not a complete zlib member yet.
      }
      if (end - start > 4096) break;
    }
  }
  return false;
}

export function inflateBgzf(buffer) {
  const parts = [];
  let offset = 0;
  while (offset + 18 <= buffer.length) {
    if (buffer[offset] !== 0x1f || buffer[offset + 1] !== 0x8b) break;
    const extraLength = buffer.readUInt16LE(offset + 10);
    let blockSize = null;
    let cursor = offset + 12;
    const extraEnd = cursor + extraLength;
    while (cursor + 4 <= extraEnd && cursor + 4 <= buffer.length) {
      const subfieldLength = buffer.readUInt16LE(cursor + 2);
      if (buffer[cursor] === 66 && buffer[cursor + 1] === 67 && subfieldLength === 2) {
        blockSize = buffer.readUInt16LE(cursor + 4);
      }
      cursor += 4 + subfieldLength;
    }
    if (blockSize == null) break;
    const blockLength = blockSize + 1;
    if (offset + blockLength > buffer.length) break;
    parts.push(zlib.gunzipSync(buffer.subarray(offset, offset + blockLength)));
    offset += blockLength;
  }
  return Buffer.concat(parts);
}

export function bytesForSlices(file, slices, readAt) {
  const parts = [];
  const ordered = [...slices].sort((left, right) => left.start - right.start);
  for (const slice of ordered) {
    const length = Number(slice.end) - Number(slice.start);
    if (length <= 0) continue;
    parts.push(readAt(file, Number(slice.start), length));
  }
  return Buffer.concat(parts);
}
