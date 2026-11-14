# Peri runtime qualification

This record binds the restricted A1 runtime to one reviewed binary. It is
not acceptance of C4, a real provider, Agent submit or cancel, or the whole
A1. The public authoring API, draft and receipt contract, OpenAPI document,
and generated client are unchanged.

The product pin is `BINARY_SHA256` in
`src/encode_pipeline/services/peri_runtime/identity.py`. The operator supplies
the file through `HELIXWEAVE_PERI_BINARY`. A digest mismatch stays unavailable
and is not started. This tree does not hard-code a machine path.

## Accepted identity

| Field | Value |
| --- | --- |
| Binary SHA256 | `5cc82caf6e04578f19ff5790fbc4eb8c9e115793ff6950a699c9e2ad28c50f47` |
| Previous accepted SHA256 | `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517` |
| Peri commit | `c7d154c9985cbe08ba25e964cb7468d15014e6d0` |
| Patch SHA256 | `06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b` |
| Cargo.lock SHA256 | `0692683359fb75073eca3d14b364071bc02058be7b11e0335ccaf2b0f314f5dd` |
| Patched files | 8 |
| rustc | 1.99.0 (b940084d7 2026-09-28) |
| cargo | 1.99.0 (5f94df478 2026-08-27) |
| Review | `review-decision.json` |

The previous digest stays in the historical A0 acceptance records. Those
records are not rewritten here. The frozen A0 build manifest is not replaced;
this directory only appends the later qualification decision and its indexes.

Tool names stay `get_authoring_schema` and `get_draft_summary`. Default
approval, prediction `off`, and the native limits stay 512MiB, swap 0, one
core, and 128 tasks.

## Evidence layers

`SHA256SUMS` covers the files in this directory except itself.

| File | What it proves |
| --- | --- |
| `native-matrix-report.json` | This binary, 23 of 23 native cases passed. The runner's `helixweave_baseline` field is the A0 closeout constant `622778f64e59f935d822d92a91bade78807dabce`, not the HelixWeave head under test. |
| `native-instance-evidence-index.json` | 25 session instances, each closed, each with an empty `cleanup_errors` list. |
| `candidate-environment-probe.json` | Effective cgroup values `memory.max=536870912`, `memory.swap.max=0`, `pids.max=128`, `cpu.max=100000 100000`. The `peri` process was inside the scope. The owner process was outside it. The same file keeps the earlier matrix failure `AF_UNIX path too long`. |
| `a1-final.log`, `a1-final.xml` | Codex independent A1 run: 1 passed, 0 failed, 0 skipped, 4.69 seconds, on head `e00130e2d5e55b3197624ef6bff9c25fd687cf87` plus only the digest rebind. |
| `host-contract.log` | 7 host-broker contract tests passed. These tests do not execute the native binary. |
| `host-transport.log` | 4 relay and stub transport tests passed inside the reviewer's own scope. |
| `source-verification.json` | The eight patched sources, prediction test source, wrapper, lock, and toolchain match the archived patch. |
| `acp-archive-verification.json` | The 2026-11-09 ACP archive checksums still match. |
| `a1-native.log`, `a1-native-2.log`, `a1-independent.log` | Reviewer attempts that failed before a test ran: bootstrap rejected `-p`, the first environment lacked pytest, and a symlink layout failed provenance. |
| `environment.json`, `independent-tests.json` | The reviewer's checkout and the commands behind those attempts. |

The historical ACP unit test remains a source-level result on the previous
binary: 1 passed, 0 failed, 0 ignored, crate `peri-acp`, 2026-11-09. It was
not rerun for this artifact, and it is not a new execution of that unit test
on `5cc82caf…`. This binary's prediction behavior is the passed
`A_prediction_paths` case in the native matrix.

## Reproduce the pinned check

Use a task-local environment created from `workflow/envs/ci-fast.lock`, install
this exact checkout there, and invoke the canonical bootstrap. A user systemd
session is required. Point `HELIXWEAVE_PERI_BINARY` at a file whose SHA256 is
the accepted digest.

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  verify-checkout
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  pytest test/services/test_peri_authoring_bridge.py::test_unverified_binary_stays_unavailable
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  pytest -m platform_real_execution \
  test/services/test_peri_authoring_bridge.py::test_fixed_native_reads_two_tools_and_rejects_an_old_generation
```

The first pytest is the digest-mismatch path. The second is the matching
binary on a real cgroup and the local model stub. Do not repeat the 23-case
matrix or the ACP compiler job unless the Peri commit, the eight-file patch,
or the isolation configuration changes.

The scientific execution closure is separate. This rebind does not change
Bulk RNA-seq execution identity inputs, so those generated manifests stay as
they are and the Protected Bulk Gate is not dispatched.
