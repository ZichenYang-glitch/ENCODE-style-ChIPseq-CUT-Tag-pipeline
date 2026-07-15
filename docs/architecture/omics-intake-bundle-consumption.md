# Omics Intake Bundle Consumption Boundary

Status: draft implementation boundary

## Decision

HelixWeave consumes Omics Intake output only through a pinned, versioned Bundle
public contract. It does not import Omics Intake Python modules, call its CLI,
inspect its persistence, or infer fields from its canonical `intake.json`.

The first slice is a service-only, read-only inspection boundary for exactly:

| Coordinate | Supported value |
| --- | --- |
| Bundle file | `intake-bundle.json` at the Bundle project root |
| Bundle schema | `0.2` |
| Schema ID | `urn:omics-intake:schema:intake-bundle:0.2` |
| Schema SHA-256 | `6dcda336c9f0ba763383ddd58bec280946f8970af8bd730eb758fab5e3a8dd71` |
| Omics Intake release | [`v0.2.0`](https://github.com/ZichenYang-glitch/omics-intake/releases/tag/v0.2.0) |
| Annotated tag object | `140a454d1313b19b322a825a1feebbb1494297c7` |
| Pinned release commit | `32680c12465f543214ed7e0173c639e0d40c7113` |
| Release tree | `48aba2f48fa88fc37dab19c10f0ce70f2641add2` |
| Public workflow name | `encode-epigenomics` |
| Render contract | `encode-render-v3` |
| HelixWeave adapter | `encode-style-chipseq-cuttag-atac-mnase` |

HelixWeave consumes the Bundle 0.2 public contract formally released with
Omics Intake `v0.2.0`. A release audit on 2026-07-15 resolved the annotated tag
object above to the pinned commit and tree, and verified that the released
schema bytes exactly match the packaged schema digest. Source provenance is
fixed to that immutable commit; runtime acceptance remains offline and pinned
to the schema version, ID, and digest rather than trusting or resolving a tag.
HelixWeave rejects every other Bundle version rather than guessing
compatibility.

The annotated tag is unsigned. Its object identity records release evidence;
it is not a publisher signature or proof of authenticity.

## Ownership

Omics Intake owns acquisition, normalization, its canonical project model,
download evidence, readiness derivation, and Bundle production. The public
Bundle schema, fixed artifact roles, file-reference records, and render
contract are the only cross-repository inputs.

HelixWeave owns local source safety, contract dispatch, consumer-side digest
verification, adapter selection, mapping into `WorkflowInputs`, and a fresh
adapter validation. The ENCODE adapter—not the workflow-neutral service—owns
the meaning of `samples.encode.tsv`, `config.encode.yaml`, and
`encode-render-v3`.

```text
versioned public Bundle
  -> workflow-neutral bounded reader
  -> pinned JSON Schema + semantic checks
  -> registered adapter pure mapper + logical file-set candidates
  -> service-owned safe candidate resolution + first observations
  -> current adapter validation
  -> service-owned second observations + private identity comparison
  -> ephemeral inspection result
```

The generic service does not contain ENCODE sample columns, assay fields,
index naming, or resource keys. A future adapter may implement the same
optional capability without adding workflow branches to the service.

## Read-only result

Inspection returns fresh in-memory values:

- the exact contract coordinate and pinned schema digest;
- a SHA-256 observation of the current Bundle document;
- producer name/version and canonical project identity from the public
  contract;
- the selected HelixWeave workflow ID;
- adapter-mapped `WorkflowInputs`; and
- relative-path, size, SHA-256, and contract-binding observations for every
  project file required by the adapter.

It does not write the Bundle, create a validated snapshot, create or start a
run, enqueue work, persist provenance, or mutate a workspace. A successful
inspection is not an execution authorization.

The current snapshot model has no durable Bundle provenance field. Persisting
this identity with a validated snapshot would require a separate domain and
persistence decision, an Alembic migration, atomicity coverage, and a public
projection decision. Bundle metadata must not be hidden in adapter options or
run tags in the meantime.

## Acceptance policy

The initial consumer is deliberately stricter than schema validity. It accepts
a Bundle only when all of the following are true:

1. The source is a normalized local path to the exact Bundle filename.
2. Every source component is opened without following symlinks; source files
   are regular, single-link files and remain identity-stable while read.
3. The Bundle is bounded, UTF-8 JSON with no duplicate keys, non-finite or
   fractional numbers, oversized integers, or excessive nesting.
4. The packaged schema bytes match the pinned schema SHA-256 and schema ID,
   then Draft 2020-12 validation with an explicit RFC 3339 `date-time` checker
   passes.
5. Collections are deterministically ordered and unique, canonical identity
   matches its record digest, and planned destinations are unambiguous.
6. The producer reports `ready` and `runnable`, both validation modes carry
   current passing evidence, no requirement is unresolved, and no `error`,
   `needs_review`, or transfer-owned issue remains.
7. The canonical record and fixed sample/config artifacts match their declared
   size and SHA-256 when re-read locally.
8. Every Bundle file reference has a project-local, producer-verified SHA-256
   binding; external or unbound handoffs are rejected.
9. The selected adapter explicitly declares and implements Bundle import,
   accepts the public workflow/render coordinate, and maps paths lexically
   beneath the Bundle root without probing the filesystem. It may declare
   ordered logical alternatives, but receives no descriptor or reader.
10. The service resolves every alternative through its descriptor-relative,
    no-follow reader. For Bowtie2, a complete standard `.bt2` set is preferred
    deterministically over a complete large `.bt2l` set; a partial standard
    set may fall back only to a complete large set. Unsafe entries never become
    a reason to try the next alternative.
11. Adapter-required paths with producer bindings match their declared size
    and SHA-256. Other adapter-required local resources are observed and marked
    unbound rather than promoted to Bundle provenance.
12. Current adapter validation succeeds. Candidate state and required files
    are observed again after validation; service-private leaf and intermediate
    directory identities must match, and the root directory binding is
    rechecked. Identity evidence is not included in the public result.

The current scientific validator retains legacy file-existence checks. Those
checks are workflow-semantic compatibility only: they are not source-safety
evidence and do not authorize a path. A successful inspection is decided only
after the service-owned no-follow observations and private identity comparison
also succeed.

The ENCODE render mapper additionally rejects YAML aliases, duplicate or
non-string mapping keys, excessive YAML depth/node counts, unexpected TSV
columns, oversized cells, absolute mapped input/resource paths, and input paths
without the expected producer file binding. Its Bundle mapping call graph has
no file existence, type, metadata, or content probes.

All failures use stable, path-free issue messages. Raw JSON/YAML/TSV values,
absolute paths, environment data, exception text, and symlink targets are not
returned in issues.

## Integrity limits

The Bundle supplies internal digest closure, not authenticity. It has no
signature or self-hash, so HelixWeave records its own SHA-256 observation but
does not claim who produced those bytes. A producer-declared `verified` digest
is also rehashed by the consumer; declared MD5 evidence alone is never treated
as a local execution binding.

The two render artifacts are contract-bound. Bundle 0.2 does not fully bind
Bowtie2 index members, genome resources, or every optional control path.
HelixWeave can check that such files are local, safe, stable, and presently
hashable, but reports them as `contract_bound = false`. Durable execution from
an imported Bundle remains deferred until provenance and snapshot ownership
for these files are explicit.

Resource consumption is bounded before content reads: the Bundle document is
limited to 16 MiB, each render artifact to 2 MiB, each adapter-required file to
256 GiB, and all required files to 1 TiB per observation pass. The second pass
after adapter validation uses the same bound. These limits prevent sparse or
oversized local paths from turning inspection into an unbounded read; they are
consumer policy rather than facts from the producer contract.

The two-pass identity check closes persistent file or directory replacement
across current adapter validation, including same-byte replacement with a new
inode. It cannot prevent a privileged actor from replacing and restoring the
same binding entirely between observations, or mutation immediately after the
ephemeral inspection returns. Inspection therefore remains evidence about a
bounded interval, not a durable execution authorization.

## Compatibility and deferred surfaces

| Surface | Initial policy |
| --- | --- |
| Bundle 0.2 local directory | Inspect through the pinned public contract |
| Bundle 0.1 | Reject as unsupported |
| Unknown future versions | Reject as unsupported |
| Ready but not runnable Bundle | Reject |
| External or unbound file reference | Reject |
| Remote URL, archive, or upload | Not implemented |
| API, CLI, or frontend import flow | Not implemented |
| Snapshot/run creation | Not implemented |
| Network access or producer subprocess | Forbidden |
| Omics Intake package dependency | Forbidden |

Supporting another Bundle version requires a new vendored schema with recorded
source commit and digest, its own validator dispatch entry, compatibility and
negative fixtures, an adapter render-contract decision, and review of semantic
differences. Updating the adjacent repository or its Python package must never
silently change HelixWeave behavior.

## Verification evidence

The boundary is covered at three layers:

- platform value tests cover immutability, defensive copies, path grammar, and
  integrity fields;
- service tests cover exact-version dispatch, schema and semantic rejection,
  digest mismatch, standard/large candidate selection, incomplete sets,
  source types, symlinks, replacement races, no source writes, fresh results,
  and workflow-neutral adapter dispatch; and
- import-boundary and packaging checks ensure no producer/runtime dependency
  leaks in and that the pinned schema ships with the distribution.

An accepted test Bundle is constructed from the public schema and public
render contract only. Tests never import Omics Intake to create or load it.
