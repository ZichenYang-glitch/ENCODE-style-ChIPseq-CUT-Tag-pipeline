import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import { catalogAdoptProposal, catalogPreviewProposal } from '../../api/generated/catalog/catalog';
import type { CatalogProposalPreviewResponse, CatalogProposalResponse } from '../../api/generated/models';
import { PeriProposalReview } from './PeriProposalReview';

vi.mock('../../api/generated/catalog/catalog', () => ({
  catalogAdoptProposal: vi.fn(),
  catalogPreviewProposal: vi.fn(),
}));

const terms = ['CTCF', 'H3K27ac'].map((name, i) => ({ term_revision_id: 'termrev' + i, display_name: name }));
const proposal: CatalogProposalResponse = {
  ok: true, proposal_id: 'proposal', project_id: 'project', contract_version: '1',
  source: 'peri', proposal_state: 'pending', adoptable_by_caller: true,
  rule_version: '1', vocabulary_version: '1', projection_version: '1', proposal_digest: 'a',
  created_at: '', expires_at: '', retention_until: '',
  items: ['CTCF', 'H3K27ac'].map((name, i) => ({
    item_id: name, target_kind: 'experiment', target_id: 'exp' + i,
    base_revision_id: 'original' + i, current_revision_id: 'original' + i,
    base_revision_number: 2, current_revision_number: 2, allowed_fields: ['target'],
    before: { target: { status: 'unknown' } }, after: { target: { status: 'known', term_id: 'term' + i, term_revision_id: 'termrev' + i } }, citations: {},
  })),
};
const preview: CatalogProposalPreviewResponse = {
  ok: true, proposal_id: 'proposal', project_id: 'project', contract_version: '1',
  proposal_digest: 'a', preview_digest: 'b', selected_item_ids: ['CTCF'], excluded_item_ids: ['H3K27ac'],
  expected_revisions: { exp0: 2 }, has_conflicts: false, conflicts: [], expires_at: '', items: [],
};

describe('Review a Peri classification', () => {
  beforeEach(() => { vi.resetAllMocks(); });

  it('submits only the previewed subset and preserves the adoption key after a lost response', async () => {
    const user = userEvent.setup();
    const onAdopted = vi.fn();
    vi.mocked(catalogPreviewProposal).mockResolvedValue(preview);
    const receipt = {
      ok: true as const, replayed: true, contract_version: '1', receipt_id: 'receipt',
      proposal_id: 'proposal', project_id: 'project', proposal_digest: 'a', preview_digest: 'b',
      actor_user_id: 'user', idempotency_key: 'key', selected_item_ids: ['CTCF'],
      excluded_item_ids: ['H3K27ac'], applied_targets: [], adopted_at: '',
    };
    vi.mocked(catalogAdoptProposal).mockRejectedValueOnce(new Error('lost')).mockResolvedValueOnce(receipt);
    render(<PeriProposalReview terms={terms} proposal={proposal} onAdopted={onAdopted} />);
    expect(screen.getByRole('button', { name: 'Preview selected' })).toBeDisabled();
    await user.click(screen.getByRole('checkbox', { name: /CTCF/ }));
    await user.click(screen.getByRole('button', { name: 'Preview selected' }));
    expect(catalogPreviewProposal).toHaveBeenCalledWith('proposal', { selected_item_ids: ['CTCF'] });
    expect(catalogAdoptProposal).not.toHaveBeenCalled();
    expect(screen.getByText('1 selected, 1 excluded.')).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: 'Adopt previewed changes' }));
    expect(await screen.findByRole('alert')).toHaveTextContent('Retry');
    await user.click(screen.getByRole('button', { name: 'Adopt previewed changes' }));
    const [first, second] = vi.mocked(catalogAdoptProposal).mock.calls;
    expect(first).toEqual(['proposal', {
      preview_digest: 'b', selected_item_ids: ['CTCF'], expected_revisions: { exp0: 2 },
    }, { 'Idempotency-Key': expect.any(String) }]);
    expect(second).toEqual(first);
    expect(onAdopted).toHaveBeenCalledWith(receipt);
  });

  it('prevents adoption when an unexpected field would hide an additional change', () => {
    render(<PeriProposalReview terms={terms} proposal={{
      ...proposal, items: [{ ...proposal.items[0]!, allowed_fields: ['target', 'batch'] }],
    }} onAdopted={vi.fn()} />);
    expect(screen.getByRole('alert')).toHaveTextContent('outside the experiment target field');
    expect(screen.queryByRole('button')).toBeNull();
  });

  it('invalidates the preview on selection changes and blocks a conflicted preview', async () => {
    const user = userEvent.setup();
    vi.mocked(catalogPreviewProposal).mockResolvedValueOnce(preview).mockResolvedValueOnce({
      ...preview, has_conflicts: true,
    });
    render(<PeriProposalReview terms={terms} proposal={proposal} onAdopted={vi.fn()} />);
    await user.click(screen.getByRole('checkbox', { name: /CTCF/ }));
    await user.click(screen.getByRole('button', { name: 'Preview selected' }));
    expect(screen.getByRole('button', { name: 'Adopt previewed changes' })).toBeInTheDocument();
    await user.click(screen.getByRole('checkbox', { name: /H3K27ac/ }));
    expect(screen.queryByRole('button', { name: 'Adopt previewed changes' })).toBeNull();
    await user.click(screen.getByRole('button', { name: 'Preview selected' }));
    expect(screen.getByRole('alert')).toHaveTextContent('changed');
    expect(screen.queryByRole('button', { name: 'Adopt previewed changes' })).toBeNull();
    expect(catalogAdoptProposal).not.toHaveBeenCalled();
  });
});
