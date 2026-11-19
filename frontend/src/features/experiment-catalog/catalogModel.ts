import type { CatalogListExperimentsParams } from '../../api/generated/models';

export const CATALOG_PAGE_SIZE = 20;

const STATUS_LABELS: Record<string, string> = {
  unknown: '未知',
  not_applicable: '不适用',
  not_detected: '未检出',
  ambiguous: '歧义',
};

export const STABLE_KEY = /^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$/;

export const RUN_INPUT_SCOPE_TEXT =
  '运行输入范围。关联样本是该次运行的输入范围，不是文件贡献者。';

export interface ReadableSlot {
  status: string;
  display_name?: string;
  term_id?: string;
  term_revision_id?: string;
  values?: ReadonlyArray<{ display_name: string }>;
}

export type NegativeStatus = 'unknown' | 'not_applicable' | 'not_detected';

export type WriteSlot =
  | { status: 'known'; term_id: string; term_revision_id: string }
  | { status: NegativeStatus };

export interface PinnedSample {
  sample_revision_id: string;
  replicate_group?: string | null;
  replicate_role: ReadableSlot;
}

export interface SampleWrite {
  sample_revision_id: string;
  replicate_group?: string;
  replicate_role: WriteSlot;
}

export function singleSlotLabel(slot: ReadableSlot | undefined): string {
  if (slot === undefined) return '未分类';
  if (slot.status === 'known') {
    return slot.display_name && slot.display_name.length > 0
      ? slot.display_name
      : '已知';
  }
  return STATUS_LABELS[slot.status] ?? slot.status;
}

export function multiSlotLabel(slot: ReadableSlot | undefined): string {
  if (slot === undefined) return '未分类';
  if (slot.status === 'known') {
    const names = (slot.values ?? [])
      .map((value) => value.display_name)
      .filter((name) => name.length > 0);
    return names.length > 0 ? names.join(', ') : '已知';
  }
  return STATUS_LABELS[slot.status] ?? slot.status;
}

export function experimentCellTypeLabel(
  samples: ReadonlyArray<{
    catalog_state?: 'absent';
    cell_type?: ReadableSlot;
  }>,
): string {
  if (samples.length === 0) return '未分类';
  const labels = samples.map((sample) =>
    sample.catalog_state === 'absent' || sample.cell_type === undefined
      ? '未分类'
      : multiSlotLabel(sample.cell_type),
  );
  return Array.from(new Set(labels)).join(', ');
}

export function writeSlot(slot: ReadableSlot): WriteSlot | null {
  if (slot.status === 'known') {
    if (
      slot.term_id === undefined ||
      slot.term_id === '' ||
      slot.term_revision_id === undefined ||
      slot.term_revision_id === ''
    ) {
      return null;
    }
    return {
      status: 'known',
      term_id: slot.term_id,
      term_revision_id: slot.term_revision_id,
    };
  }
  if (
    slot.status === 'unknown' ||
    slot.status === 'not_applicable' ||
    slot.status === 'not_detected'
  ) {
    return { status: slot.status };
  }
  return null;
}

export function samplesForWrite(
  samples: readonly PinnedSample[],
): SampleWrite[] | null {
  const written: SampleWrite[] = [];
  for (const sample of samples) {
    const role = writeSlot(sample.replicate_role);
    if (role === null) return null;
    const item: SampleWrite = {
      sample_revision_id: sample.sample_revision_id,
      replicate_role: role,
    };
    if (
      typeof sample.replicate_group === 'string' &&
      sample.replicate_group.length > 0
    ) {
      item.replicate_group = sample.replicate_group;
    }
    written.push(item);
  }
  return written.length > 0 ? written : null;
}

export function knownTermSlot(
  terms: ReadonlyArray<{ term_id: string; term_revision_id: string }>,
  termId: string,
): Extract<WriteSlot, { status: 'known' }> | null {
  const term = terms.find((item) => item.term_id === termId);
  if (term === undefined) return null;
  return {
    status: 'known',
    term_id: term.term_id,
    term_revision_id: term.term_revision_id,
  };
}

export function experimentListParams(input: {
  cursor?: string;
  assayTermIds: readonly string[];
  cellTypeTermIds: readonly string[];
}): CatalogListExperimentsParams {
  const params: CatalogListExperimentsParams = { limit: CATALOG_PAGE_SIZE };
  if (input.cursor !== undefined && input.cursor !== '') {
    params.cursor = input.cursor;
  }
  if (input.assayTermIds.length > 0) {
    params.assay_term_id = [...input.assayTermIds];
  }
  if (input.cellTypeTermIds.length > 0) {
    params.cell_type_term_id = [...input.cellTypeTermIds];
  }
  return params;
}

export function headRevision<T extends { revision_number: number }>(
  revisions: readonly T[],
): T | null {
  return revisions.reduce<T | null>(
    (best, item) =>
      best === null || item.revision_number > best.revision_number ? item : best,
    null,
  );
}

export function editingHead<T extends { experiment_id: string; revision_number: number }>(
  listRow: T | null,
  pinned: T | null,
): T | null {
  if (
    pinned !== null &&
    (listRow === null ||
      (pinned.experiment_id === listRow.experiment_id &&
        pinned.revision_number >= listRow.revision_number))
  ) {
    return pinned;
  }
  return listRow;
}

export function revisionChoices<T extends { experiment_revision_id: string }>(
  extras: readonly (T | null)[],
  page: readonly T[],
): T[] {
  const seen = new Set(page.map((item) => item.experiment_revision_id));
  const leading: T[] = [];
  for (const extra of extras) {
    if (extra === null || seen.has(extra.experiment_revision_id)) continue;
    seen.add(extra.experiment_revision_id);
    leading.push(extra);
  }
  return [...leading, ...page];
}

export function experimentTraceText(trace: {
  trace_status: string;
  reason?: string | null;
}): { title: string; detail: string } {
  if (trace.trace_status === 'unresolved') {
    const reason =
      trace.reason === 'experiment_snapshot_binding_absent'
        ? 'experiment_snapshot_binding_absent'
        : (trace.reason ?? 'unresolved');
    return {
      title: '未绑定',
      detail: `这条实验修订没有绑定运行（${reason}）。这里不列出产物，也不把运行样本写成文件贡献者。`,
    };
  }
  return {
    title: trace.trace_status,
    detail: '实验追溯状态不是未绑定。若出现产物行，它仍只表示运行输入范围。',
  };
}

export function newIdempotencyHeaders(): { 'Idempotency-Key': string } {
  const token = `${Date.now().toString(36)}${Math.random().toString(36).slice(2, 12)}`;
  return { 'Idempotency-Key': `cat${token}`.slice(0, 128) };
}

export async function collectPages<T>(
  load: (cursor: string | undefined) => Promise<{
    items: readonly T[];
    nextCursor: string | null | undefined;
  }>,
): Promise<{ items: T[]; truncated: boolean }> {
  const items: T[] = [];
  let cursor: string | undefined;
  for (let page = 0; page < 20; page += 1) {
    const result = await load(cursor);
    items.push(...result.items);
    if (result.nextCursor == null || result.nextCursor === '') {
      return { items, truncated: false };
    }
    cursor = result.nextCursor;
  }
  return { items, truncated: true };
}
