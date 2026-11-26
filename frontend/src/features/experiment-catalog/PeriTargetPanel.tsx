import { useState, type FormEvent } from 'react';
import { catalogGenerateClassificationSuggestion } from '../../api/generated/catalog/catalog';
import type {
  CatalogAdoptionReceiptResponse,
  CatalogExperimentResponse,
  CatalogProposalResponse,
  CatalogTermSummaryResponse,
} from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';
import { Button } from '../../components/Button';
import { useUiLanguage } from '../../lib/uiLanguage';
import { newIdempotencyHeaders } from './catalogModel';
import { PeriProposalReview } from './PeriProposalReview';

interface Props {
  projectId: string;
  experiment: CatalogExperimentResponse;
  terms: readonly CatalogTermSummaryResponse[];
  onAdopted: (receipt: CatalogAdoptionReceiptResponse) => void;
}

export function PeriTargetPanel({ projectId, experiment, terms, onAdopted }: Props) {
  const english = useUiLanguage() === 'en';
  const [clue, setClue] = useState('');
  const [termRevision, setTermRevision] = useState('');
  const [headers, setHeaders] = useState(newIdempotencyHeaders);
  const [busy, setBusy] = useState(false);
  const [proposal, setProposal] = useState<CatalogProposalResponse | null>(null);
  const [insufficient, setInsufficient] = useState(false);
  const [error, setError] = useState('');
  const targets = terms.filter((term) => term.term_type === 'target' && term.archived_at == null);

  function resetRequest() {
    setHeaders(newIdempotencyHeaders());
    setProposal(null);
    setInsufficient(false);
    setError('');
  }

  async function generate(event: FormEvent) {
    event.preventDefault();
    if (busy || !clue.trim()) return;
    setError('');
    setBusy(true);
    try {
      const result = await catalogGenerateClassificationSuggestion(projectId, experiment.experiment_id, {
        experiment_revision_id: experiment.experiment_revision_id,
        field: 'target',
        reason_code: 'TARGET_TERM_MATCH',
        clue: clue.trim(),
        ...(termRevision && { vocabulary_term_revision_id: termRevision }),
      }, headers);
      if (result.experiment_revision_id !== experiment.experiment_revision_id || result.experiment_id !== experiment.experiment_id) {
        setError(english ? 'The suggestion belongs to a different revision.' : '建议对应另一份实验修订。');
      } else if (result.outcome === 'insufficient_evidence') {
        setInsufficient(true);
      } else if (result.proposal?.project_id === projectId) {
        setProposal(result.proposal);
      } else {
        setError(english ? 'Could not load the proposal for this project.' : '无法读取当前课题的建议。');
      }
    } catch (cause) {
      setError(cause instanceof ApiError
        ? cause.message
        : (english ? 'Could not retrieve the suggestion. You can retry the same request.'
          : '未收到建议，可用同一请求重试。'));
    } finally {
      setBusy(false);
    }
  }

  return (
    <section aria-label={english ? 'Peri target suggestions' : 'Peri 靶标建议'}
      className="space-y-3 rounded-[6px] border border-[var(--color-border)] bg-[var(--color-surface)] p-4">
      <h3 className="text-base font-semibold">{english ? 'Suggest an experiment target' : '建议实验靶标'}</h3>
      <form onSubmit={generate} className="space-y-2">
        <fieldset disabled={busy} className="space-y-2">
          <label className="block text-sm">
            {english ? 'Target clue' : '靶标线索'}
            <input className="mt-1 block rounded border px-2 py-1" value={clue} minLength={2} maxLength={64}
              placeholder="H3K27ac" required autoComplete="off" onChange={(event) => {
                setClue(event.target.value);
                resetRequest();
              }} />
          </label>
          <label className="block text-sm">
            {english ? 'Candidate target' : '候选靶标'}
            <select className="mt-1 block rounded border px-2 py-1" value={termRevision} onChange={(event) => {
              setTermRevision(event.target.value);
              resetRequest();
            }}>
              <option value="">{english ? 'Use the active target vocabulary' : '使用当前靶标词表'}</option>
              {targets.map((term) => (
                <option key={term.term_revision_id} value={term.term_revision_id}>{term.display_name}</option>
              ))}
            </select>
          </label>
          <p className="text-sm text-[var(--color-text-muted)]">{english
            ? 'Your clue and candidate target labels go to the bound model. Use a target name or alias.'
            : '靶标线索和候选词条名称会发送给已绑定模型。请输入靶标名称或别名。'}</p>
          <Button type="submit" disabled={busy || clue.trim().length < 2 || !!proposal || insufficient}>
            {busy ? (english ? 'Requesting…' : '请求中…') : (english ? 'Ask Peri' : '向 Peri 请求建议')}
          </Button>
        </fieldset>
      </form>
      {error && <p role="alert">{error}</p>}
      {insufficient && <p role="status">{english
        ? 'Insufficient evidence to suggest a target. Refine the clue or candidate.'
        : '证据不足，无法建议靶标。可调整线索或候选词。'}</p>}
      {proposal && <PeriProposalReview key={proposal.proposal_id} proposal={proposal} terms={targets} onAdopted={onAdopted} />}
    </section>
  );
}
