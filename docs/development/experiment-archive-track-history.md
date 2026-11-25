# Open a sealed result version from an experiment

A researcher selects an experiment, selects a result version and opens it in the embedded WashU viewer. The link identifies a sealed copy, including after a later run publishes a file with the same name. Sealing and retention are local operator actions. The browser is read-only. This feature builds on the independent viewer in PR #265.

## Identity and publication

SQLite stores an immutable version ID, project, stable experiment, source experiment revision, stable track key, result version number, display name, assembly, format and exact source. Version numbers increase within an experiment and track key, independently of experiment metadata revisions. The experiment page includes results across its revisions and identifies the source revision. An unsealed result has no link. The authorized terminal operator explicitly chooses a registered result's experiment. Native run publications must already be bound to that experiment revision. A native source pins run ID, artifact ID, generation and revision. A registered source pins result ID and registered revision. An explicitly indexed track pins both members and their exact sources. One seal represents one track.

A terminal command resolves an exact source through the existing byte contract, verifies it and copies it once to an operator-managed absolute directory (`ENCODE_PIPELINE_ARCHIVE_ROOT`). SQLite records size, full SHA-256, fixed-size chunk hashes and a sealed file identity for each member. Publish a complete manifest only after every member is durable. Retrying the same seal returns the same identity. Concurrency cannot reuse a version number. The source can subsequently move. A reader pins and verifies the opened descriptor and the chunks it actually reads. A missing/offline copy or changed copy returns an explicit error. All reads use the requested version. Paths, digests and operator instructions stay off the public response.

## Read contract

Use `/api/v1/catalog/projects/{project_id}/experiments/{experiment_id}/archived-tracks` as the collection. GET accepts `limit` (default 50, max 100) and `offset` (default 0), returning `{items, next_offset}`. GET `/{archive_version_id}` returns the same descriptor as a list item. Fields are `archive_version_id`, `project_id`, `experiment_id`, `experiment_revision_id`, `track_key`, `version_number`, `display_name`, `assembly` (`hg38`, `mm10`, or other explicitly supported WashU assembly), `track_type` (`bigwig`, `bed` or `longrange`), `source` (string-valued object with `kind` equal to `run` or `registered` plus exact identifiers) and `has_index`. Versioned same-origin GET/HEAD paths end in `/{archive_version_id}/content` and, for indexed longrange, `/content.tbi`. Use a single standard HTTP byte range, including suffix ranges, with 206, `Content-Range`, 416 for unsatisfiable requests, and private no-store caching. Each request enforces the existing authorized reader/session policy and project/experiment/source binding. A stale run URL retains its current behavior. An offline copy returns HTTP 409 with `archive_offline`. A changed copy returns HTTP 409 with `archive_integrity`. An unknown or out-of-scope version returns 404; a missing or revoked login returns 401. These conditions do not select a different version.

The catalog action opens `/viewer.html?project=…&experiment=…&archive=…` in the existing WashU viewer. The viewer resolves the descriptor before constructing same-origin byte URLs. It shows source, explicit result version and assembly, and uses the existing gene/coordinate navigation. The browser supplies no filename, path, external URL or upload.

## Storage and operation

Configure the same absolute `ENCODE_PIPELINE_ARCHIVE_ROOT` on the operator command and the API serving the experiment page. Create the directory on the long-term volume first. The operator needs create access and the API process needs read/traverse access. The sealed objects use SHA-256 names and mode 0440. The catalog supplies project, experiment, track and version names. If a file already exists with the same digest, sealing verifies its bytes before reuse. A directory that is missing, unmounted or resolves through a symlink fails closed. A changed object is withheld from the viewer. After a byte-identical disk restore, an operator can repeat the original seal to verify the restored copy and rebind its physical identity. The result keeps its version and URL. Include both SQLite and the sealed objects in the archive backup. Keep the existing retention procedure for original inputs and scratch.

Use the same track key for successive publications of one signal. For example, for a registered result, run the installed operator command once per verified revision.

```sh
helixweave archive seal registered \
  --project-id "$project_id" --experiment-id "$experiment_id" \
  --experiment-revision-id "$experiment_revision_id" --track-key signal \
  --display-name signal.bw --assembly hg38 --track-type bigwig \
  --result-id "$result_id" --result-revision-id "$result_revision_id"
```

The terminal receipt identifies the exact version. A researcher signs in, opens **Experiments**, selects the experiment, then chooses **在 WashU 查看第 1 版** on the desired track. The viewer opens in the assembly's default window. Use the gene or coordinate field to move to the locus of interest. Sealed bigWig and indexed longrange use byte ranges. Plain BED supports nonempty UTF-8 intervals up to 8 MiB. The archived copy can outlive the scratch file. Changes to inode, link count or permissions require operator verification of the full restored bytes before the old address can serve them again.

## Acceptance

From the packaged experiment page, an authorized colleague can select either of two identically named sealed bigWigs and read the distinct expected signal at one locus (4.5 for v1 and 9 for v2). Each synthetic bigWig exceeds 16 MiB. Capture HTTP ranges, transferred bytes and expected/rendered value, including v1 after v2 exists. Check missing and same-size changed copies, missing index, unauthorized/wrong-scope requests, exact registered and run identities, and existing registered/run reads. Keep separate evidence for any trial using laboratory data.
