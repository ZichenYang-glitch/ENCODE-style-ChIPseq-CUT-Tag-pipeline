# A1 C1 handoff

Status: C1 and C2 are implemented on this branch and are back with Codex for
review. The public turn response now includes the server `proposal_digest` and,
for an explanation, `explanation.reason_codes`. Do not recompute the digest in
the frontend. This handoff is not a frozen C3 contract until that review
finishes, and it is not C4 or A1 acceptance.

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
