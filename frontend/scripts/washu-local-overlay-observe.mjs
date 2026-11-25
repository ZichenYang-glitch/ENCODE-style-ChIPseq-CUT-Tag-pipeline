// Pinned GenomeHub 62.0.2 renders numerical bins to canvas and places BED
// and PET features in the React drawing components. Read the geometry those
// components received, after clipping and layout. Names alone cannot pass.
export function renderedWindow(expected) {
  const host = document.querySelector('[data-testid="washu-genome-viewer"]');
  if (!host) return null;
  const rootKey = Object.keys(host).find(key => key.startsWith('__reactFiber'));
  let root = host[rootKey];
  if (!root) return null;
  while (root.return) root = root.return;
  function* subtree(node) {
    if (!node) return;
    yield node;
    for (let child = node.child; child; child = child.sibling) yield* subtree(child);
  }
  // The DOM's fiber can point to the previous render after navigation.
  let viewer;
  for (const node of subtree(root.stateNode.current)) {
    if (node.stateNode === host) {
      viewer = node;
      break;
    }
  }
  const hits = {};
  for (const fiber of subtree(viewer)) {
    const props = fiber.memoizedProps;
    if (!props) continue;
    // The visible ruler uses a half-open interval. Search uses a 1-based start.
    if (expected.region && props.trackModel?.type === 'ruler' && props.windowWidth > 0) {
      const intervals = props.selectedRegion?.getGenomeIntervals?.();
      if (intervals?.length === 1) {
        const { chr, start, end } = intervals[0];
        if (`${chr}:${start + 1}-${end}` === expected.region) hits.viewport = { chr, start, end };
      }
    }
    for (const key of ['data', 'features', 'items', 'placements', 'arrangedData']) {
      const entries = props[key];
      if (!Array.isArray(entries) || !entries.length) continue;
      if (props.viewRegion?.getGenomeIntervals && entries.includes(expected.signal.value)) {
        const interval = props.viewRegion.getGenomeIntervals().find(r => r.chr === expected.chrom);
        if (interval) {
          const pixel = (interval.end - interval.start) / entries.length;
          const first = interval.start + entries.indexOf(expected.signal.value) * pixel;
          const last = interval.start + (entries.lastIndexOf(expected.signal.value) + 1) * pixel;
          if (Math.abs(first - expected.signal.start) <= 3 * pixel &&
              Math.abs(last - expected.signal.end) <= 3 * pixel) {
            hits.signal = { value: expected.signal.value, first, last, bpPerPixel: pixel };
          }
        }
      }
      for (const item of entries) {
        if (!item || typeof item !== 'object') continue;
        const a = item.interaction?.locus1;
        const b = item.interaction?.locus2;
        if (item.points?.length && a?.chr === expected.chrom && b?.chr === expected.chrom &&
            a.start === expected.pet[0] && a.end === expected.pet[1] &&
            b.start === expected.pet[2] && b.end === expected.pet[3] &&
            item.xSpan1?.end > item.xSpan1?.start && item.xSpan2?.end > item.xSpan2?.start) {
          hits.pet = { a, b, score: item.interaction.score };
        }
        const feature = item.feature;
        if (item.placedFeatures?.length && feature?.chr === expected.chrom &&
            feature.start === expected.peak[0] && feature.end === expected.peak[1] &&
            item.xSpan?.end > item.xSpan?.start) {
          hits.peak = { chr: feature.chr, start: feature.start, end: feature.end, name: feature[3] };
        }
      }
    }
  }
  if (expected.region && !hits.viewport) return null;
  return hits.pet && hits.signal && hits.peak ? hits : null;
}

export function renderedChromosome(chromosome) {
  const host = document.querySelector('[data-testid="washu-genome-viewer"]');
  if (!host) return false;
  const rootKey = Object.keys(host).find(key => key.startsWith('__reactFiber'));
  const pending = [host[rootKey]];
  while (pending.length) {
    const fiber = pending.pop();
    if (!fiber) continue;
    if (fiber.child) pending.push(fiber.child);
    if (fiber.sibling) pending.push(fiber.sibling);
    if (fiber.memoizedProps?.viewRegion?.getGenomeIntervals?.().some(r => r.chr === chromosome)) {
      return true;
    }
  }
  return false;
}
