# M2 catalog adoption over Slice A

An administrator registers an offline suggestion, previews one subset, and
adopts that subset in a single Slice A transaction. Suggestions are local
synthetic proposals. This path does not call a model.

The reviewed preview contract stays in
`docs/architecture/experiment-catalog-m2-preview-contract.md`. HTTP adoption
uses that digest and citation shape unchanged.

## Requests

All routes are under `/api/v1/catalog`. Writes send `Idempotency-Key`.
The raw body is hashed with the Slice A scope `(actor_user_id, idempotency_key)`,
then compared with the operation, the proposal id, and the body hash.

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
