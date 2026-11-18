import { describe, expect, it } from 'vitest';
import {
  CATALOG_PAGE_SIZE,
  RUN_INPUT_SCOPE_TEXT,
  experimentCellTypeLabel,
  experimentListParams,
  experimentTraceText,
  headRevision,
  knownTermSlot,
  multiSlotLabel,
  samplesForWrite,
  singleSlotLabel,
  writeSlot,
} from './catalogModel';

describe('catalog display labels', () => {
  it('uses the contract labels and keeps known names pinned', () => {
    expect(singleSlotLabel({ status: 'known', display_name: 'CUT&Tag' })).toBe(
      'CUT&Tag',
    );
    expect(singleSlotLabel({ status: 'unknown' })).toBe('未知');
    expect(singleSlotLabel({ status: 'not_applicable' })).toBe('不适用');
    expect(singleSlotLabel({ status: 'not_detected' })).toBe('未检出');
    expect(singleSlotLabel({ status: 'ambiguous' })).toBe('歧义');
    expect(singleSlotLabel(undefined)).toBe('未分类');
    expect(
      multiSlotLabel({
        status: 'known',
        values: [{ display_name: 'regulatory T cell' }],
      }),
    ).toBe('regulatory T cell');
    expect(
      experimentCellTypeLabel([
        { catalog_state: 'absent' },
        {
          cell_type: {
            status: 'known',
            values: [{ display_name: 'keratinocyte' }],
          },
        },
      ]),
    ).toBe('未分类, keratinocyte');
  });

  it('describes an unresolved experiment trace as unbound', () => {
    expect(
      experimentTraceText({
        trace_status: 'unresolved',
        reason: 'experiment_snapshot_binding_absent',
      }),
    ).toEqual({
      title: '未绑定',
      detail:
        '这条实验修订没有绑定运行（experiment_snapshot_binding_absent）。这里不列出产物，也不把运行样本写成文件贡献者。',
    });
    expect(RUN_INPUT_SCOPE_TEXT).toContain('运行输入范围');
    expect(RUN_INPUT_SCOPE_TEXT).toContain('不是文件贡献者');
  });
});

describe('catalog write and query shaping', () => {
  it('copies known and negative slots and refuses ambiguous ones', () => {
    expect(
      writeSlot({
        status: 'known',
        term_id: 'term_a',
        term_revision_id: 'trmv_a',
        display_name: 'CUT&Tag',
      }),
    ).toEqual({
      status: 'known',
      term_id: 'term_a',
      term_revision_id: 'trmv_a',
    });
    expect(writeSlot({ status: 'unknown' })).toEqual({ status: 'unknown' });
    expect(
      writeSlot({
        status: 'ambiguous',
        term_id: 'term_a',
        term_revision_id: 'trmv_a',
      }),
    ).toBeNull();
    expect(
      samplesForWrite([
        {
          sample_revision_id: 'smpr_1',
          replicate_group: 'biorep-1',
          replicate_role: {
            status: 'known',
            term_id: 'term_role',
            term_revision_id: 'trmv_role',
          },
        },
      ]),
    ).toEqual([
      {
        sample_revision_id: 'smpr_1',
        replicate_group: 'biorep-1',
        replicate_role: {
          status: 'known',
          term_id: 'term_role',
          term_revision_id: 'trmv_role',
        },
      },
    ]);
    expect(
      samplesForWrite([
        {
          sample_revision_id: 'smpr_1',
          replicate_role: { status: 'ambiguous' },
        },
      ]),
    ).toBeNull();
  });

  it('asks for one page of project, assay, and cell type filters', () => {
    expect(
      experimentListParams({
        assayTermIds: ['term_assay'],
        cellTypeTermIds: ['term_cell', 'term_other'],
      }),
    ).toEqual({
      limit: CATALOG_PAGE_SIZE,
      assay_term_id: ['term_assay'],
      cell_type_term_id: ['term_cell', 'term_other'],
    });
    expect(
      experimentListParams({
        cursor: 'catcur_next',
        assayTermIds: [],
        cellTypeTermIds: [],
      }).cursor,
    ).toBe('catcur_next');
    expect(CATALOG_PAGE_SIZE).toBe(20);
  });

  it('picks the current head and a known term revision', () => {
    expect(
      headRevision([
        { revision_number: 1, name: 'old' },
        { revision_number: 3, name: 'head' },
        { revision_number: 2, name: 'middle' },
      ]),
    ).toEqual({ revision_number: 3, name: 'head' });
    expect(
      knownTermSlot(
        [{ term_id: 'term_a', term_revision_id: 'trmv_a' }],
        'term_a',
      ),
    ).toEqual({
      status: 'known',
      term_id: 'term_a',
      term_revision_id: 'trmv_a',
    });
    expect(knownTermSlot([], 'missing')).toBeNull();
  });
});
