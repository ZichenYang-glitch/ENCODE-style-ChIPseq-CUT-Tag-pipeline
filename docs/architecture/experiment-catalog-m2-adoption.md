# M2 catalog adoption over Slice A

An administrator can register an offline suggestion or ask the bound Peri model
to suggest an experiment target. Both sources produce a pending proposal.
The administrator selects the changes, previews them, and explicitly adopts the
previewed subset in one Slice A transaction.

The reviewed preview contract stays in
`docs/architecture/experiment-catalog-m2-preview-contract.md`. HTTP adoption
uses that digest and citation shape unchanged.

## Requests

All routes are under `/api/v1/catalog`. Writes send `Idempotency-Key`.
The raw body is hashed with the Slice A scope `(actor_user_id, idempotency_key)`,
then compared with the operation, the proposal id, and the body hash.

- `POST /projects/{project_id}/experiments/{experiment_id}/classification-suggestions`
  asks Peri about the selected experiment revision. An owner grant for this
  project, model binding, and `experiment.target` authorizes the target clue
  and up to four controlled-vocabulary labels. Automatic candidates come
  from current target revisions recorded on active experiments in that
  project. An administrator may explicitly select one current target from
  the shared vocabulary, including the first target for a new project.
  An accepted answer yields a
  pending `peri_bound_model` proposal pinned to the experiment and vocabulary
  revisions. An inconclusive answer yields `insufficient_evidence`.
  The registration transaction checks the experiment head, target revision
  and curation versions again. Generation and preview preserve the experiment
  revision. After another write advances the experiment, generating against
  the old revision returns a conflict. Its existing proposal and adoption
  receipt remain available at their own endpoints.
  The owner setup and grant command are documented in
  `integrations/peri-bound/README.md`.
- `POST /projects/{project_id}/proposals` registers one suggestion. The
  proposal id is derived from the actor, the idempotency key, and the project,
  so a lost response retries the same receipt. A retry returns the frozen
  registration view.
- `GET /proposals/{proposal_id}` is the live view. Members can read it.
  Ambiguous candidates and admin-only citations are removed. Preview and
  adoption stay administrator-only, including a role check inside the write
  transaction.
- `POST /proposals/{proposal_id}/preview` returns the subset digest and does
  not append a revision.
- `POST /proposals/{proposal_id}/adoptions` appends revisions only for the
  selected items and stores one adoption receipt. The same key and body
  replay that receipt with `replayed: true`. A second administrator, a
  changed body, a stale head, a bad citation, or a retired rule version
  writes nothing.
- `GET /projects/{project_id}/adoption-receipts` lists committed receipts.

Active rule, vocabulary, and projection versions live in the single
`catalog_curation_versions` row. There is no public route to change them.
Tests update that row through the catalog repository.

## Revision graph

The pull request is based on current `main`. The reviewed IRENE-27 contract
matches `8d621ed395f2c266646f038f9a8edc2335eba397`.

Alembic revision `20261119_32` revises the existing `20261119_18` snapshot
binding revision. That revision keeps its id, filename, and history. The
graph has one head. Adoption does not rewrite the Slice B binding tables, so
rows written at `20261119_18` stay readable after the upgrade.
