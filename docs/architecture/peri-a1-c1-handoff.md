# A1 C1 handoff

Status: C1 and C2 remain on this branch. C3-C1 receipt discovery remains.
This round adds an explicit ENCODE workflow sample plan so validation can
build sample rows from admitted revisions. The contract commit is
`f6cd097e86455243bfcb461651e807e4da970483`. It is not byte-identical to
`7ee16d706049e4d5f6abea20b01891b48c7b16b4`. The turn response
still includes the server `proposal_digest` and, for an explanation,
`explanation.reason_codes`. Do not recompute the digest in the frontend. This
handoff is not a frozen C3 contract, and it is not C4 or A1 acceptance.

Base: `b2f15868c6f74ab8764a73ae6fe4394f67304475` (`origin/main` on 2026-11-10).
Branch: `agent/peri-a1-authoring`.

## What the frontend can call

All routes are under `/api/v1`, use the existing session cookie and CSRF header,
and reject unknown JSON fields. Errors use `{ok:false, reason_code, issues}`
with fixed sentences. Validator business results keep the existing validation
envelope plus `draft_revision`.

| Method and path | Purpose |
| --- | --- |
| `GET /projects` | Bounded project id and display name. Legacy and archived projects are omitted. |
| `GET /projects/{id}/samples` | Sample id, display name, and current revision id. |
| `GET /agent/authoring-sessions` | The caller's drafts. Optional `project_id`, `workflow_id`, `cursor`, `limit` (default 20, max 50). |
| `POST /agent/authoring-sessions` | Create an empty draft for one project and workflow. |
| `GET /agent/authoring-sessions/{sid}` | Draft or `draft: null` after expiry. No renewal. |
| `GET /agent/authoring-sessions/{sid}/receipts` | Bounded adopted receipts for this draft. Optional `cursor` and `limit` (default 20, max 50). |
| `PUT /agent/authoring-sessions/{sid}/draft` | CAS save. Renewal only when normalized content changes. |
| `POST /agent/authoring-sessions/{sid}/connection` | Explicit assistant connection. |
| `DELETE /agent/authoring-sessions/{sid}/connection` | Close the current connection. Idempotent. |
| `POST /agent/authoring-sessions/{sid}/turns` | Structured question. |
| `GET /agent/authoring-sessions/{sid}/turns/{turn_id}` | Turn, proposal visibility, and receipt. |
| `POST /agent/authoring-sessions/{sid}/proposals/{proposal_id}/adopt` | Adopt or replay the same body. |
| `POST /agent/authoring-sessions/{sid}/proposals/{proposal_id}/reject` | Reject a pending proposal. |
| `POST /agent/authoring-sessions/{sid}/validate` | Load the stored draft and call the existing validator. |

`authoring_capability` is the field and issue catalog. Do not hard-code adapter
fields. Display names from the project list stay in the UI.

Static examples live in `docs/architecture/peri-a1-c1-examples.json`.

## Local fixtures

`reference_ready_app` plus `seed_test_authentication` in `test/conftest.py`
already create a login and an ENCODE reference. Create a user project with
`DataRegistryService` before `POST /agent/authoring-sessions`. The API test
client in `test/api/api_test_client.py` attaches the cookie and CSRF header.

By default `POST .../connection` returns 503 `BRIDGE_UNAVAILABLE`. Saving and
`POST .../validate` still work. One owner process admits the fixed Peri binary
from `HELIXWEAVE_PERI_BINARY` and listens on `HELIXWEAVE_PERI_OWNER_SOCKET`.
API workers only set the socket variable. Each native runs in its own 512MiB,
swap-0, 128-task scope; the API process is not placed in that scope. The
binary path is not built into the product.

A completed proposal response includes `proposal_digest`. Adopt and reject send
that value back. An explanation has `explanation.reason_codes`, no
`proposal_id`, and `adoptable: false`. There is no free-text explanation field.

The capability token is handed to the owner once at admission and is not an
API, worker, or browser field. A later turn from any authorized API worker
sends an empty token; the owner checks the stored token. Loss or timeout of an
in-flight open or turn is fenced by the owner, then the matching native is
closed. A normal end-of-response socket close does not fence the session.

## Finding an adopted receipt again

`GET /agent/authoring-sessions/{sid}/receipts` is how the owner finds an
adopted receipt after a refresh, a draft switch, or a new login. The operation
id is `authoringListReceipts`. The generated client function is
`authoringListReceipts`. Checking this list does not adopt again and does not
open an agent connection.

The response is `draft_state`, `receipts`, and `next_cursor`. Each receipt item
has `turn_id` and `adoption_receipt` (`proposal_id`, `base_revision`,
`result_revision`, `changed`, `adopted_at`, `draft_expires_at`). A static
example is `receipt_list_response` in `docs/architecture/peri-a1-c1-examples.json`.

- An active draft returns `draft_state: active`.
- After `draft_expires_at`, including the exact timestamp, the draft is
  `expired` and the receipts remain until `purge_after`. The session read still
  returns `draft: null`.
- At `purge_after`, and after it, the response is 404 even if the database row
  is still present.
- A draft with no adopted receipt returns 200 and `receipts: []`.
- Order is `adopted_at` descending, then `turn_id` descending. `limit` defaults
  to 20 and cannot exceed 50. Pass `next_cursor` back as `cursor`. A cursor
  this endpoint did not produce is 422 `INVALID_REQUEST`.
- Another user, an unknown sid, or a project that is no longer visible is 404.
- The list does not include the question, the originating login, model text,
  paths, logs, or any proposal that was not adopted.
- Show the adopted field difference with the existing
  `GET .../turns/{turn_id}` using the returned `turn_id`. That read does not
  reconnect the assistant.
- Listing does not change revision, expiry, or the stored receipt. The original
  adopt request can still be repeated and returns `replayed: true`.
- Logout, closing the connection, or retiring a generation leaves the adopted
  receipt in place. The new login sees `requires_new_connection: true` until
  the user explicitly connects again.

## Workflow sample plan

Selecting `sample_revision_ids` does not by itself create an ENCODE sample
row. `samples` stays `null`. The draft also carries `workflow_samples`, an
explicit plan the server turns into rows only while validating.

Each plan item names one `sample_revision_id`, `assay`, `target`, `peak_mode`,
`layout`, `role`, optional `control_sample_revision_id`, and `reads`. A read
is `1` or `2` plus an `input_file_revision_id`. SE has read 1. PE has reads 1
and 2. A treatment that names a control must point at another row whose role
is `control`. The same input revision cannot be used twice. The order of
`workflow_samples` must equal `sample_revision_ids`.

The server checks that every revision belongs to the draft's project, that
the input file is still the admitted bytes, and then the ENCODE adapter builds
the workflow rows. Host paths are not stored in the draft, returned by the
API, or shown to the model. A successful validation returns the existing
snapshot object. `snapshot.sample_revision_ids` is the bound selection.

Errors, all with `snapshot: null` unless noted:

- No sample revision: `SAMPLE_SELECTION_REQUIRED`.
- Revisions selected without a plan: `WORKFLOW_SAMPLE_PLAN_REQUIRED`.
- Plan shape, pairing, role, or project does not match: save returns 422, or
  validate returns `WORKFLOW_SAMPLE_PLAN_INVALID`.
- The admitted file changed or cannot be verified:
  `WORKFLOW_SAMPLE_INPUT_UNAVAILABLE`.
- A later validation of an unchanged draft returns the same snapshot only
  after the admitted bytes still match.

An example plan is `workflow_sample_plan` in
`docs/architecture/peri-a1-c1-examples.json`. Receipt discovery is unchanged.
This is not C4 or A1 acceptance.

## Not in this handoff

- The bridge can return one structured proposal or explanation for the approved
  fields. It does not accept free text, and it still has no preflight, create,
  start, or cancel tool.
- The assistant can read only `config.threads` and `options.cores` on the ENCODE
  workflow. Other workflows report `available: false`.
- There is no free-text question field and no preflight, create, start, or cancel
  route on these endpoints.
- A new login can edit an unexpired draft. It cannot adopt a proposal created by
  the previous login.
- C3 pages and the C4 browser loop are still open. This is not A1 acceptance.

## Frontend-owned paths

Pages, feature components, and frontend unit or browser tests. Do not edit
`frontend/src/api/generated/` or `frontend/openapi.json`. Ask the backend session
to change the contract and regenerate the client.
