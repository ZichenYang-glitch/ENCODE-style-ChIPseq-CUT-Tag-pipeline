import { useState } from 'react';
import {
  catalogAdoptProposal,
  catalogPreviewProposal,
} from '../../api/generated/catalog/catalog';
import type {
  CatalogAdoptionReceiptResponse,
  CatalogProposalPreviewResponse,
  CatalogProposalResponse,
  CatalogTermSummaryResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import { useUiLanguage } from '../../lib/uiLanguage';
import { newIdempotencyHeaders } from './catalogModel';

type TargetTerm = Pick<CatalogTermSummaryResponse, 'term_revision_id' | 'display_name'>;

function targetLabel(value: unknown, terms: readonly TargetTerm[], english: boolean): string {
  const unavailable = english ? 'Unavailable' : '无可用值';
  if (!value || typeof value !== 'object') return unavailable;
  if ('display_name' in value && typeof value.display_name === 'string') return value.display_name;
  if ('term_revision_id' in value && typeof value.term_revision_id === 'string') {
    const revision = value.term_revision_id;
    return terms.find((term) => term.term_revision_id === revision)?.display_name ??
      (english ? `Vocabulary revision ${revision}` : `词条修订 ${revision}`);
  }
  if ('status' in value && typeof value.status === 'string') {
    if (value.status === 'unknown') return english ? 'unknown' : '未分类';
    return value.status;
  }
  return unavailable;
}

export function PeriProposalReview({ proposal, terms, onAdopted }: {
  proposal: CatalogProposalResponse;
  terms: readonly TargetTerm[];
  onAdopted: (receipt: CatalogAdoptionReceiptResponse) => void;
}) {
  const english = useUiLanguage() === 'en';
  const [selection, setSelection] = useState<string[]>([]);
  const [preview, setPreview] = useState<CatalogProposalPreviewResponse | null>(null);
  const [adoptionHeaders, setAdoptionHeaders] = useState(newIdempotencyHeaders);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  const [receipt, setReceipt] = useState<CatalogAdoptionReceiptResponse | null>(null);

  function toggle(itemId: string, checked: boolean) {
    setSelection((ids) => checked ? [...ids, itemId] : ids.filter((id) => id !== itemId));
    setPreview(null);
    setError('');
  }

  async function previewSelected() {
    setBusy(true);
    setError('');
    try {
      const next = await catalogPreviewProposal(proposal.proposal_id, { selected_item_ids: selection });
      setPreview(next);
      setAdoptionHeaders(newIdempotencyHeaders());
    } catch (cause) {
      setPreview(null);
      setError(cause instanceof ApiError ? cause.message : (english ? 'Could not load the preview. Try again.' : '无法载入预览，请重试。'));
    } finally {
      setBusy(false);
    }
  }

  async function adoptSelected() {
    if (!preview || preview.has_conflicts) return;
    setBusy(true);
    setError('');
    try {
      const result = await catalogAdoptProposal(proposal.proposal_id, {
        preview_digest: preview.preview_digest,
        expected_revisions: preview.expected_revisions,
        selected_item_ids: preview.selected_item_ids,
      }, adoptionHeaders);
      setReceipt(result);
      onAdopted(result);
    } catch (cause) {
      // Keep the preview and its key after a lost response so retry retrieves the receipt.
      if (cause instanceof ApiError && [403, 409, 410].includes(cause.status)) setPreview(null);
      setError(cause instanceof ApiError ? cause.message : (english ? 'Could not confirm adoption. Retry with the same preview.' : '未收到采用回执，可用同一预览重试。'));
    } finally {
      setBusy(false);
    }
  }

  const supported = proposal.items.every((item) =>
    item.target_kind === 'experiment' && item.allowed_fields.length === 1 && item.allowed_fields[0] === 'target',
  );
  if (!supported) return <p role="alert">{english
    ? 'This proposal contains changes outside the experiment target field.'
    : '这份建议含有实验靶标以外的变更。'}</p>;
  if (receipt) return <p role="status">{english ? 'Adopted. Receipt' : '已采用。回执'} {receipt.receipt_id}</p>;
  if (!proposal.adoptable_by_caller) return <p role="status">{english ? 'Proposal status' : '建议状态'} · {proposal.proposal_state}</p>;
  return (
    <div className="space-y-2">
      <fieldset disabled={busy} className="space-y-1">
        <legend className="font-medium">{english ? 'Choose target changes to preview' : '选择要预览的靶标变更'}</legend>
        {proposal.items.map((item) => (
          <label key={item.item_id} className="flex items-center gap-2">
            <input type="checkbox" checked={selection.includes(item.item_id)}
              onChange={(event) => toggle(item.item_id, event.target.checked)} />
            {english ? 'Target' : '靶标'} · {targetLabel(item.before.target, terms, english)} → {targetLabel(item.after.target, terms, english)} ({english ? 'revision' : '修订'} {item.base_revision_number})
          </label>
        ))}
      </fieldset>
      <Button type="button" disabled={busy || selection.length === 0} onClick={previewSelected}>
        {english ? 'Preview selected' : '预览所选建议'}
      </Button>
      {preview && (
        <div className="space-y-1">
          {preview.has_conflicts ? (
            <p role="alert">{english ? 'The experiment or its terms changed. Review the current revision before generating a new suggestion.' : '实验或词条已变更。请查看当前修订后再生成建议。'}</p>
          ) : (
            <>
              <p>{english
                ? `${preview.selected_item_ids.length} selected, ${preview.excluded_item_ids.length} excluded.`
                : `已选择 ${preview.selected_item_ids.length} 项，排除 ${preview.excluded_item_ids.length} 项。`}</p>
              <ul>
                {preview.items.map((item) => (
                  <li key={item.item_id}>
                    {targetLabel(item.before.target, terms, english)} → {targetLabel(item.after.target, terms, english)}
                  </li>
                ))}
              </ul>
              <Button type="button" variant="primary" disabled={busy} onClick={adoptSelected}>
                {english ? 'Adopt previewed changes' : '采用预览中的变更'}
              </Button>
            </>
          )}
        </div>
      )}
      {error && <p role="alert">{error}</p>}
    </div>
  );
}
