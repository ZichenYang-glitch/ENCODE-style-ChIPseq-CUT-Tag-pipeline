# A1 C1 handoff

Status: C1 and C2 remain on this branch. C3-C1 receipt discovery remains.
The sample-plan contract began at
`f6cd097e86455243bfcb461651e807e4da970483`. This repair keeps the accepted
receipt contract at `7ee16d706049e4d5f6abea20b01891b48c7b16b4` and stores the
selected input revisions in the snapshot input-use binding. A managed
snapshot is validation evidence, not an execution qualification. The turn response
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
row. `samples` stays `null`. The draft also carries `workflow_samples`.

Identity moves like this:

1. The browser sends stable revision ids and business fields. It does not send
   host paths.
2. Save keeps those ids in the draft. The platform checks only bounds, sample
   and input revision id shape, uniqueness, and that plan order matches
   `sample_revision_ids`.
3. For the ENCODE workflow, the adapter then checks assay, target, peak mode,
   layout, role, read pairing, and control relationships. A scientific
   mismatch is `422` `INVALID_REQUEST` on save.
4. Validate resolves each input revision in the draft's project, rejects
   archived or cross-project revisions, and re-reads the admitted bytes
   through the existing descriptor walk. The same check runs for the public
   workflow validate route.

What this validation actually establishes:

- The ENCODE strict check, when `options.strict_inputs` is true, only checks
  that the submitted path is a regular file. It does not parse FASTQ bases.
- Immediately before that check, the service opens the selected revision and
  requires the submitted path to be that same leaf with the admitted bytes.
- That directory entry stays watched through the scientific check, the
  input-use declaration, and the pre-save re-read. A rename, replacement, or
  rewrite in that interval fails the validation even if the original inode
  and bytes are put back before the function returns.
- The snapshot then records the revision id and the content hash confirmed by
  those reads, plus the path string that passed them. It is evidence for this
  validation, not a promise that the file stays unchanged afterward.
- The database insert itself is not a filesystem transaction. Managed
  execution remains unavailable; a later run must recheck the file.
5. The same revision ids are passed to the existing input-use binding. The
   ENCODE adapter declares one `regular_file_v1` use per occupied read
   (`encode-sample-reads-v1`). A validation that does not name input
   revisions leaves those reads transitional and stores no member. When the
   plan names every read, the repository authorizes the project, pool, and
   revision digest, then stores the member's revision id, digest, size, and
   content hash.
6. The public snapshot's `input_binding.input_uses[].members` carries that
   revision id. `fully_managed: true` means every declared read has an
   admitted revision. It does not qualify execution. `verify_run` still
   returns `MANAGED_INPUT_EXECUTION_UNAVAILABLE`.

The generic draft model does not enumerate ENCODE assays or peak modes. The
browser reads the workflow's own field rules from `GET /workflows/{workflow_id}/schema`.
For this adapter the plan object is still:

- `sample_revision_id`
- `assay`: `chipseq`, `cuttag`, `atac`, or `mnase`
- `target`
- `peak_mode`: `narrow`, `broad`, or `nucleosome`
- `layout`: `PE` or `SE`
- `role`: `treatment` or `control`
- `control_sample_revision_id`: another plan row, or `null`
- `reads`: `{read: 1|2, input_file_revision_id}`

SE is read 1. PE is reads 1 and 2. ATAC uses `narrow`. MNase uses `PE` and
`nucleosome`. A named control must be a row with `role: control`. The same
input revision cannot be used twice.

Save, read, and digest rules:

- `PUT .../draft` stores the canonical plan, including
  `control_sample_revision_id: null` when the browser omits it.
- `GET` returns that same object on `draft.workflow_samples`.
- The draft digest covers `workflow_samples`. A frontend that compares or
  resubmits the whole draft must include this field with the other canonical
  keys: `config`, `input_selections`, `options`,
  `reference_profile_revision_id`, `sample_revision_ids`, and `samples`.
- The draft revision increases only when the canonical content changes.

Where the ids come from:

- `sample_revision_id` comes from `GET /api/v1/projects/{project_id}/samples`
  (`authoringListProjectSamples`).
- `input_file_revision_id` comes from
  `GET /api/v1/projects/{project_id}/input-revisions`
  (`authoringListProjectInputRevisions`). The browser copies
  `revisions[].input_file_revision_id`. It does not invent an `inpfr_` value
  or send a host path.
- `reference_profile_revision_id` is the enabled profile revision already
  selected for the workflow. `control_sample_revision_id` is the
  `sample_revision_id` of another row in the same plan whose `role` is
  `control`.

Input revision catalog:

- Login and project visibility match the sample list. An unknown, archived,
  or non-user project is `404` `NOT_FOUND`. A bad cursor is `422`
  `INVALID_REQUEST`.
- `200` with `revisions: []` and `next_cursor: null` is an empty selectable
  catalog. That includes a project with no active storage-pool binding.
- Each item is the current revision of one non-archived file:
  `input_file_id`, `stable_key`, `input_file_revision_id`,
  `revision_number`, `size_bytes`, `content_sha256`.
- Archived files and older revisions are omitted. The response has no storage
  root, relative path, or private configuration.
- Order is `input_file_id` ascending. `limit` defaults to 20 and is at most
  50. `next_cursor` is opaque; pass it back unchanged.
- The call does not admit a file, write a draft, or start a native or
  scientific task. Do not send this catalog to the model. The Peri tool
  allowlist is unchanged.

A saved draft's previous snapshot is not reused when `workflow_samples`,
`sample_revision_ids`, `reference_profile_revision_id`, `config`, `options`,
or `samples` change, or when the admitted bytes or stored input-use members
no longer match. The full draft comparison is `config`, `input_selections`,
`options`, `reference_profile_revision_id`, `sample_revision_ids`, `samples`,
and `workflow_samples`.

The verified public sequence is `catalog_validation_loop` in
`docs/architecture/peri-a1-c1-examples.json`. Receipt discovery is unchanged.

Errors, all with `snapshot: null` unless noted:

- No sample revision: `SAMPLE_SELECTION_REQUIRED`.
- Revisions selected without a plan: `WORKFLOW_SAMPLE_PLAN_REQUIRED`.
- Plan shape, pairing, role, archived revision, or project does not match:
  save returns 422, or validate returns `WORKFLOW_SAMPLE_PLAN_INVALID`.
- The admitted file changed, was replaced, or cannot be verified:
  `WORKFLOW_SAMPLE_INPUT_UNAVAILABLE`.
- A later validation of an unchanged draft returns the same snapshot only
  when the admitted bytes still match and the stored input-use members are
  still the selected revisions.

An example plan is `workflow_sample_plan` in
`docs/architecture/peri-a1-c1-examples.json`. Receipt discovery is unchanged.
This is not an execution handoff, not C4, and not A1 acceptance.

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
