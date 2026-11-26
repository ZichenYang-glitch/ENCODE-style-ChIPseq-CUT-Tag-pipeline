import { describe, expect, it } from 'vitest';
import {
  CATALOG_PAGE_SIZE,
  RUN_INPUT_SCOPE_TEXT,
  catalogConflictText,
  editingHead,
  emptySampleRunsText,
  experimentCellTypeLabel,
  experimentListParams,
  experimentTraceText,
  headRevision,
  knownTermSlot,
  multiSlotLabel,
  revisionChoices,
  revisionMarkLabel,
  revisionNumberLabel,
  runInputScopeHeading,
  runInputScopeText,
  runSampleLine,
  samplesForWrite,
  singleSlotLabel,
  writeSlot,
} from './catalogModel';

describe('catalog display labels', () => {
  it('uses the contract labels and keeps known names pinned', () => {
    expect(singleSlotLabel({ status: 'known', display_name: 'CUT&Tag' }, 'zh-CN')).toBe(
      'CUT&Tag',
    );
    expect(singleSlotLabel({ status: 'unknown' }, 'zh-CN')).toBe('未知');
    expect(singleSlotLabel({ status: 'not_applicable' }, 'zh-CN')).toBe('不适用');
    expect(singleSlotLabel({ status: 'not_detected' }, 'zh-CN')).toBe('未检出');
    expect(singleSlotLabel({ status: 'ambiguous' }, 'zh-CN')).toBe('歧义');
    expect(singleSlotLabel(undefined, 'zh-CN')).toBe('未分类');
    expect(
      multiSlotLabel(
        {
          status: 'known',
          values: [{ display_name: 'regulatory T cell' }],
        },
        'zh-CN',
      ),
    ).toBe('regulatory T cell');
    expect(
      experimentCellTypeLabel(
        [
          { catalog_state: 'absent' },
          {
            cell_type: {
              status: 'known',
              values: [{ display_name: 'keratinocyte' }],
            },
          },
        ],
        'zh-CN',
      ),
    ).toBe('未分类, keratinocyte');
    expect(revisionNumberLabel(2, 'zh-CN')).toBe('修订 2');
    expect(revisionMarkLabel(true, 24, 'zh-CN')).toBe('当前 24');
    expect(revisionMarkLabel(false, 1, 'zh-CN')).toBe('旧版 1');
    expect(runSampleLine('smp_treg', 2, 'zh-CN')).toBe('运行样本 smp_treg 修订 2');
    expect(catalogConflictText('zh-CN')).toBe('冲突。这次提交没有写入。');
    expect(emptySampleRunsText('zh-CN')).toBe('这条样本修订没有运行结果。');
    expect(runInputScopeHeading('zh-CN')).toBe('运行输入范围');
    expect(runInputScopeText('zh-CN')).toBe(RUN_INPUT_SCOPE_TEXT);
  });

  it('describes an unresolved experiment trace as unbound', () => {
    expect(
      experimentTraceText(
        {
          trace_status: 'unresolved',
          reason: 'experiment_snapshot_binding_absent',
        },
        'zh-CN',
      ),
    ).toEqual({
      title: '未绑定',
      detail:
        '这条实验修订没有绑定运行（experiment_snapshot_binding_absent）。这里不列出产物，也不把运行样本写成文件贡献者。',
    });
    expect(
      experimentTraceText({ trace_status: 'resolved', reason: null }, 'zh-CN'),
    ).toEqual({
      title: 'resolved',
      detail: '实验追溯状态不是未绑定。若出现产物行，它仍只表示运行输入范围。',
    });
    expect(RUN_INPUT_SCOPE_TEXT).toContain('运行输入范围');
    expect(RUN_INPUT_SCOPE_TEXT).toContain('不是文件贡献者');
  });

  it('uses natural English and keeps researcher names and API values', () => {
    expect(singleSlotLabel({ status: 'known', display_name: 'CUT&Tag' }, 'en')).toBe(
      'CUT&Tag',
    );
    expect(
      singleSlotLabel({ status: 'known', display_name: '调节性T细胞-批次A' }, 'en'),
    ).toBe('调节性T细胞-批次A');
    expect(singleSlotLabel({ status: 'known', display_name: '' }, 'en')).toBe('Known');
    expect(singleSlotLabel({ status: 'unknown' }, 'en')).toBe('Unknown');
    expect(singleSlotLabel({ status: 'not_applicable' }, 'en')).toBe('Not applicable');
    expect(singleSlotLabel({ status: 'not_detected' }, 'en')).toBe('Not detected');
    expect(singleSlotLabel({ status: 'ambiguous' }, 'en')).toBe('Ambiguous');
    expect(singleSlotLabel({ status: 'custom_api_status' }, 'en')).toBe(
      'custom_api_status',
    );
    expect(singleSlotLabel(undefined, 'en')).toBe('Unclassified');
    expect(
      multiSlotLabel(
        {
          status: 'known',
          values: [{ display_name: 'regulatory T cell' }, { display_name: '角质形成细胞' }],
        },
        'en',
      ),
    ).toBe('regulatory T cell, 角质形成细胞');
    expect(
      experimentCellTypeLabel(
        [
          { catalog_state: 'absent' },
          {
            cell_type: {
              status: 'known',
              values: [{ display_name: 'keratinocyte' }],
            },
          },
        ],
        'en',
      ),
    ).toBe('Unclassified, keratinocyte');
    expect(
      experimentTraceText(
        {
          trace_status: 'unresolved',
          reason: 'experiment_snapshot_binding_absent',
        },
        'en',
      ),
    ).toEqual({
      title: 'Unbound',
      detail:
        'This experiment revision has no bound run (experiment_snapshot_binding_absent). Its trace contains no artifacts or file-level sample attributions.',
    });
    expect(experimentTraceText({ trace_status: 'resolved' }, 'en')).toEqual({
      title: 'resolved',
      detail:
        'Each artifact lists the samples used as inputs to its run.',
    });
    expect(runInputScopeHeading('en')).toBe('Run input scope');
    expect(runInputScopeText('en')).toBe(
      'These samples were inputs to this run. File-level contributions are not recorded here.',
    );
    expect(revisionMarkLabel(true, 2, 'en')).toBe('Current revision 2');
    expect(revisionMarkLabel(false, 1, 'en')).toBe('Earlier revision 1');
    expect(revisionNumberLabel(2, 'en')).toBe('Revision 2');
    expect(runSampleLine('smp_treg', 2, 'en')).toBe('Run sample smp_treg revision 2');
    expect(emptySampleRunsText('en')).toBe('This sample revision has no run results.');
    expect(catalogConflictText('en')).toBe('Revision conflict. Your changes were not saved.');
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

  it('keeps a newer pinned head and still offers revisions from another page', () => {
    const listRow = { experiment_id: 'exp_1', revision_number: 401 };
    const pinned = { experiment_id: 'exp_1', revision_number: 402 };
    expect(editingHead(listRow, null)).toEqual(listRow);
    expect(editingHead(listRow, pinned)).toEqual(pinned);
    expect(editingHead({ ...listRow, revision_number: 403 }, pinned)).toEqual({
      experiment_id: 'exp_1',
      revision_number: 403,
    });
    const page = [{ experiment_revision_id: 'expr_1' }];
    expect(
      revisionChoices(
        [{ experiment_revision_id: 'expr_401' }, page[0] ?? null],
        page,
      ).map((item) => item.experiment_revision_id),
    ).toEqual(['expr_401', 'expr_1']);
  });
});
