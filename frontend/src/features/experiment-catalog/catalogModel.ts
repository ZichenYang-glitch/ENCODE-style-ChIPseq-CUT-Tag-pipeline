import type { CatalogListExperimentsParams } from '../../api/generated/models';

export const CATALOG_PAGE_SIZE = 20;

export type CatalogCopyLanguage = 'en' | 'zh-CN';

const STATUS_LABELS: Record<CatalogCopyLanguage, Record<string, string>> = {
  'zh-CN': {
    unknown: '未知',
    not_applicable: '不适用',
    not_detected: '未检出',
    ambiguous: '歧义',
  },
  en: {
    unknown: 'Unknown',
    not_applicable: 'Not applicable',
    not_detected: 'Not detected',
    ambiguous: 'Ambiguous',
  },
};

export const STABLE_KEY = /^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$/;

export const RUN_INPUT_SCOPE_TEXT =
  '运行输入范围。关联样本是该次运行的输入范围，不是文件贡献者。';

const RUN_INPUT_SCOPE_TEXT_EN =
  'Run input scope. Associated samples are the input scope of that run, not file contributors.';

function copyLanguage(lang?: CatalogCopyLanguage): CatalogCopyLanguage {
  return lang === 'en' ? 'en' : 'zh-CN';
}

function unclassifiedLabel(lang?: CatalogCopyLanguage): string {
  return copyLanguage(lang) === 'en' ? 'Unclassified' : '未分类';
}

function knownFallbackLabel(lang?: CatalogCopyLanguage): string {
  return copyLanguage(lang) === 'en' ? 'Known' : '已知';
}

function statusLabel(status: string, lang?: CatalogCopyLanguage): string {
  return STATUS_LABELS[copyLanguage(lang)][status] ?? status;
}

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

export function singleSlotLabel(
  slot: ReadableSlot | undefined,
  lang?: CatalogCopyLanguage,
): string {
  if (slot === undefined) return unclassifiedLabel(lang);
  if (slot.status === 'known') {
    return slot.display_name && slot.display_name.length > 0
      ? slot.display_name
      : knownFallbackLabel(lang);
  }
  return statusLabel(slot.status, lang);
}

export function multiSlotLabel(
  slot: ReadableSlot | undefined,
  lang?: CatalogCopyLanguage,
): string {
  if (slot === undefined) return unclassifiedLabel(lang);
  if (slot.status === 'known') {
    const names = (slot.values ?? [])
      .map((value) => value.display_name)
      .filter((name) => name.length > 0);
    return names.length > 0 ? names.join(', ') : knownFallbackLabel(lang);
  }
  return statusLabel(slot.status, lang);
}

export function experimentCellTypeLabel(
  samples: ReadonlyArray<{
    catalog_state?: 'absent';
    cell_type?: ReadableSlot;
  }>,
  lang?: CatalogCopyLanguage,
): string {
  if (samples.length === 0) return unclassifiedLabel(lang);
  const labels = samples.map((sample) =>
    sample.catalog_state === 'absent' || sample.cell_type === undefined
      ? unclassifiedLabel(lang)
      : multiSlotLabel(sample.cell_type, lang),
  );
  return Array.from(new Set(labels)).join(', ');
}

export function revisionNumberLabel(
  revisionNumber: number,
  lang?: CatalogCopyLanguage,
): string {
  return copyLanguage(lang) === 'en'
    ? `Revision ${revisionNumber}`
    : `修订 ${revisionNumber}`;
}

export function revisionMarkLabel(
  isCurrent: boolean,
  revisionNumber: number,
  lang?: CatalogCopyLanguage,
): string {
  if (copyLanguage(lang) === 'en') {
    return isCurrent
      ? `Current revision ${revisionNumber}`
      : `Earlier revision ${revisionNumber}`;
  }
  return isCurrent ? `当前 ${revisionNumber}` : `旧版 ${revisionNumber}`;
}

export function runInputScopeHeading(lang?: CatalogCopyLanguage): string {
  return copyLanguage(lang) === 'en' ? 'Run input scope' : '运行输入范围';
}

export function runInputScopeText(lang?: CatalogCopyLanguage): string {
  return copyLanguage(lang) === 'en' ? RUN_INPUT_SCOPE_TEXT_EN : RUN_INPUT_SCOPE_TEXT;
}

export function emptySampleRunsText(lang?: CatalogCopyLanguage): string {
  return copyLanguage(lang) === 'en'
    ? 'This sample revision has no run results.'
    : '这条样本修订没有运行结果。';
}

export function runSampleLine(
  sampleId: string,
  revisionNumber: number,
  lang?: CatalogCopyLanguage,
): string {
  return copyLanguage(lang) === 'en'
    ? `Run sample ${sampleId} revision ${revisionNumber}`
    : `运行样本 ${sampleId} 修订 ${revisionNumber}`;
}

export function catalogConflictText(lang?: CatalogCopyLanguage): string {
  return copyLanguage(lang) === 'en'
    ? 'Conflict. This submission was not written.'
    : '冲突。这次提交没有写入。';
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

export function experimentTraceText(
  trace: {
    trace_status: string;
    reason?: string | null;
  },
  lang?: CatalogCopyLanguage,
): { title: string; detail: string } {
  if (trace.trace_status === 'unresolved') {
    const reason =
      trace.reason === 'experiment_snapshot_binding_absent'
        ? 'experiment_snapshot_binding_absent'
        : (trace.reason ?? 'unresolved');
    if (copyLanguage(lang) === 'en') {
      return {
        title: 'Unbound',
        detail: `This experiment revision has no bound run (${reason}). Artifacts are not listed here, and run samples are not described as file contributors.`,
      };
    }
    return {
      title: '未绑定',
      detail: `这条实验修订没有绑定运行（${reason}）。这里不列出产物，也不把运行样本写成文件贡献者。`,
    };
  }
  if (copyLanguage(lang) === 'en') {
    return {
      title: trace.trace_status,
      detail:
        'The experiment trace status is not unbound. If artifact rows appear, they still only describe the run input scope.',
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
